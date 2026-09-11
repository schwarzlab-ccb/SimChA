"""Regression tests for backward dependencies, cancellation and stable-ID replay."""
import csv
import json
import random
import subprocess
import sys
import tempfile
import unittest
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hidden_event_replay import verify_filtered
from lib_contig_versions import backward_dependencies
from hidden_events import verify_retained_index
from lib_event_dependencies import event_dependencies, protect_groups
from lib_hidden_events import Arc, WGD, cancellation_groups, observed_edges, parse_regions, read_profiles, reduce_sample


class History:
    """Independent tiny simulator: each contig is a list of individual reference bases."""
    def __init__(self, chromosomes=1, size=10, sample='s'):
        self.lengths = {f'chr{i + 1}': size for i in range(chromosomes)}
        self.contigs = [[(chrom, hap, p) for p in range(size)]
                        for hap in ('H1', 'H2') for chrom in self.lengths]
        self.records = []
        self.sample = sample

    def encode(self, bases):
        counts = Counter(bases)
        regions = []
        for chrom, length in self.lengths.items():
            for hap in ('H1', 'H2'):
                cn = [counts[chrom, hap, i] for i in range(length)]
                for layer in range(1, max(cn, default=0) + 1):
                    start = None
                    for p in range(length + 1):
                        covered = p < length and cn[p] >= layer
                        if covered and start is None:
                            start = p
                        if not covered and start is not None:
                            regions.append(f'{hap}:{chrom}[{start}:{p})')
                            start = None
        return '[' + ','.join(regions) + ']'

    def event(self, kind, ci=0, start=None, end=None, *, cb=None, pos_b=None, inverted_b=False):
        gained, lost = [], []
        before = len(self.contigs)
        if kind == WGD:
            desc = 'WGD'
            gained = [base for contig in self.contigs for base in contig]
            self.contigs.extend([list(contig) for contig in self.contigs])
        elif kind == 'Skip':
            desc = ''
        elif kind == 'Translocation':
            desc = f'contigA:{ci};contigB:{cb};posA:{start};posB:{pos_b};invertedB:{inverted_b}'
            if inverted_b:
                self.contigs[cb].reverse()
            a, b = self.contigs[ci], self.contigs[cb]
            self.contigs[ci], self.contigs[cb] = a[:start] + b[pos_b:], b[:pos_b] + a[start:]
        else:
            length = len(self.contigs[ci])
            start = 0 if start is None else start
            end = length if end is None else end
            assert 0 <= start < end <= length
            desc = f'contig:{ci};length:{length};'
            if not kind.startswith(('Chrom', 'Contig')):
                desc += f'start:{start};end:{end}'
            touched = list(self.contigs[ci][start:end])
            if kind in {'ChromDuplication', 'ContigDuplication', 'ArmDuplication'}:
                self.contigs.append(touched)
                gained = touched
            elif kind.endswith('Duplication'):
                self.contigs[ci][end:end] = touched[::-1] if kind == 'InvertedDuplication' else touched
                gained = touched
            elif kind.endswith('Deletion'):
                del self.contigs[ci][start:end]
                lost = touched
            elif kind == 'InternalInversion':
                self.contigs[ci][start:end] = touched[::-1]
            else:
                raise ValueError(kind)
        self.records.append(dict(sample_id=self.sample, depth=str(len(self.records) + 1),
                                 event_type=kind, description=desc,
                                 contig_slots_before=str(before), contig_slots_after=str(len(self.contigs)),
                                 regions_gained=self.encode(gained), regions_lost=self.encode(lost)))
        return self

    def profile(self):
        counts = Counter(base for contig in self.contigs for base in contig)
        return {chrom: [(p, p + 1, counts[chrom, 'H1', p], counts[chrom, 'H2', p])
                        for p in range(length)] for chrom, length in self.lengths.items()}


class HiddenEventsTest(unittest.TestCase):
    def reduce(self, history, keep_skip=False):
        expected = observed_edges(history.profile())
        labels, retained, arcs, signals = reduce_sample(history.sample, history.records,
                                                        history.profile(), keep_skip=keep_skip)
        verify_retained_index(history.records, labels, retained)
        verify_filtered(retained, history.lengths, expected)
        hidden = {r['event_id'] for r in labels if r['hidden']}
        if labels and labels[0]['retention_method'] == 'whole_event_cancellation':
            balance = defaultdict(Counter)
            for arc in arcs:
                if arc['event_id'] in hidden:
                    units = arc['sign'] * arc['wgd_multiplier']
                    balance[arc['chrom'], arc['haplotype']][arc['start']] += units
                    balance[arc['chrom'], arc['haplotype']][arc['end']] -= units
            self.assertFalse(any(value for unit in balance.values() for value in unit.values()))
            for label in labels:
                if label['hidden']:
                    self.assertEqual(label['retained_flow_units'], 0)
                    self.assertTrue(set(json.loads(label['cancellation_group_members'])) <= hidden)
        for label in labels:
            if not label['hidden']:
                self.assertFalse(set(json.loads(label['dependency_event_ids'])) & hidden)
        return labels, retained, arcs, signals

    def test_copy_then_delete_original_are_jointly_hidden(self):
        h = History().event('ChromDuplication').event('ChromDeletion')
        labels, retained, _, _ = self.reduce(h)
        self.assertTrue(all(r['hidden'] for r in labels))
        self.assertEqual(retained, [])

    def test_cross_wgd_cancellation_keeps_only_wgd(self):
        h = History().event('InternalDuplication', start=2, end=5).event(WGD)
        h.event('InternalDeletion', 0, 5, 8).event('InternalDeletion', 2, 5, 8)
        labels, retained, _, _ = self.reduce(h)
        self.assertEqual([r['hidden'] for r in labels], [True, False, True, True])
        self.assertEqual([r['event_type'] for r in retained], [WGD])
        # The logged WGD delta contains the now-removed extra copy. It must not
        # prevent replaying the unchanged WGD command on the filtered state.
        self.assertIn('H1:chr1[2:5)', retained[0]['regions_gained'])

    def test_partial_cross_wgd_cancellation_keeps_full_events(self):
        h = History().event('InternalDuplication', start=2, end=5).event(WGD)
        h.event('InternalDeletion', 0, 5, 8)
        labels, retained, _, _ = self.reduce(h)
        self.assertFalse(any(r['hidden'] for r in labels))
        self.assertEqual(json.loads(retained[0]['effective_cn_intervals']), [['chr1', 'H1', 2, 5, 2]])

    def test_visible_tail_deletion_is_retained(self):
        h = History().event('ChromDuplication').event('TelomereDeletion', 0, 5, 10)
        labels, retained, _, _ = self.reduce(h)
        self.assertFalse(any(r['hidden'] for r in labels))
        self.assertTrue(all(r['retained_flow_units'] > 0 for r in labels))
        self.assertEqual(json.loads(retained[0]['effective_cn_intervals']), [['chr1', 'H1', 0, 10, 1]])

    def test_two_visible_tail_boundaries_after_wgd(self):
        h = History().event('ChromDuplication').event(WGD)
        h.event('TelomereDeletion', 0, 8, 10).event('TelomereDeletion', 3, 6, 10)
        labels, _, _, _ = self.reduce(h)
        self.assertFalse(any(r['hidden'] for r in labels))

    def test_visible_internal_dip_is_retained(self):
        h = History().event('ChromDuplication').event('InternalDeletion', 0, 3, 7)
        labels, _, _, _ = self.reduce(h)
        self.assertFalse(any(r['hidden'] for r in labels))

    def test_adjacent_gains_keep_both_original_operations(self):
        h = History().event('InternalDuplication', 0, 0, 5).event('InternalDuplication', 0, 10, 15)
        labels, _, _, _ = self.reduce(h)
        self.assertFalse(any(r['hidden'] for r in labels))

    def test_nonmonotone_path_does_not_hide_its_contributors(self):
        h = History().event('InternalDuplication', 0, 0, 7).event('InternalDuplication', 0, 10, 17)
        h.event('InternalDeletion', 0, 3, 7)
        labels, _, _, _ = self.reduce(h)
        self.assertFalse(any(r['hidden'] for r in labels))

    def test_gain_erased_by_whole_deletion_is_removed(self):
        h = History().event('InternalDuplication', 0, 3, 7).event('ChromDeletion')
        labels, _, _, _ = self.reduce(h)
        self.assertEqual([r['hidden'] for r in labels], [True, False])
        self.assertEqual(labels[0]['reason'], 'erased_contig_history')

    def test_duplication_dependency_keeps_entire_cancelled_group(self):
        h = History().event('InternalDuplication', 0, 2, 5).event('ChromDuplication')
        h.event('InternalDeletion', 0, 5, 8)
        labels, _, _, _ = self.reduce(h)
        self.assertFalse(any(r['hidden'] for r in labels))
        self.assertIn('s:e1', json.loads(labels[1]['dependency_event_ids']))

    def test_later_append_keeps_its_id_after_hidden_birth(self):
        h = History(chromosomes=2).event('ChromDuplication', 0).event('ChromDeletion', 4)
        h.event('ChromDuplication', 1).event('InternalDuplication', 5, 2, 5)
        labels, _, _, _ = self.reduce(h)
        self.assertEqual([r['hidden'] for r in labels], [True, True, False, False])

    def test_wgd_copy_ids_survive_hidden_birth(self):
        h = History(chromosomes=2).event('ChromDuplication', 0).event('ChromDeletion', 4).event(WGD)
        h.event('InternalDuplication', 8, 2, 5)
        labels, _, _, _ = self.reduce(h)
        self.assertEqual([r['hidden'] for r in labels], [True, True, False, False])

    def test_neutral_translocation_without_consumers_is_hidden(self):
        h = History(chromosomes=2).event('Translocation', 0, 4, cb=1, pos_b=6)
        labels, retained, _, _ = self.reduce(h)
        self.assertTrue(labels[0]['hidden'])
        self.assertEqual(retained, [])

    def test_translocation_is_retained_before_derivative_loss(self):
        for inverted in (False, True):
            h = History(chromosomes=2).event('Translocation', 0, 4, cb=1, pos_b=6, inverted_b=inverted)
            h.event('ChromDeletion', 0)
            labels, _, _, _ = self.reduce(h)
            self.assertFalse(any(r['hidden'] for r in labels))

    def test_inversion_is_retained_before_coordinate_based_deletion(self):
        h = History().event('InternalInversion', 0, 2, 8).event('InternalDeletion', 0, 1, 5)
        labels, _, _, _ = self.reduce(h)
        self.assertFalse(any(r['hidden'] for r in labels))
        self.assertEqual(labels[0]['reason'], 'required_dependency')

    def test_ambiguous_cancellation_is_deterministic_and_replayable(self):
        h = History().event('ChromDuplication').event('ChromDuplication').event('ChromDeletion')
        a = self.reduce(h)
        b = self.reduce(h)
        self.assertEqual(a, b)
        self.assertTrue(any(r['flow_choice_ambiguous'] for r in a[0]))

    def test_cycles_on_different_chromosomes_join_through_shared_event(self):
        arcs = [Arc(0, 0, 'chr1', 'H1', 0, 10, 1, 2),
                Arc(1, 0, 'chr2', 'H1', 0, 10, 1, 2),
                Arc(2, 1, 'chr1', 'H1', 0, 10, -1, 2),
                Arc(3, 2, 'chr2', 'H1', 0, 10, -1, 2)]
        self.assertEqual(cancellation_groups(arcs, dict.fromkeys(range(4), 0), 3), {0: [0, 1, 2]})

    def test_dependency_closure_keeps_group_partners_and_their_dependencies(self):
        groups = {0: [0, 3], 1: [1, 4], 2: [2]}
        self.assertEqual(protect_groups(groups, {2}, [set(), set(), {0}, {1}, set()]), {0, 1, 2})

    def test_malformed_deltas_and_unsupported_dependency_commands_fail(self):
        for value in ('', None, '[bad]', '[H1:chr1[2:2)]', '[H1:chr1[1:2)],bad'):
            with self.assertRaises(ValueError):
                parse_regions(value)
        h = History().event('ChromDuplication')
        h.records[0]['description'] = 'bad'
        with self.assertRaises(ValueError):
            self.reduce(h)
        h.records[0]['event_type'] = 'UnimplementedEvent'
        with self.assertRaisesRegex(ValueError, 'Unsupported event'):
            self.reduce(h)

    def test_random_histories_retain_valid_original_commands(self):
        for seed in range(50):
            rng = random.Random(seed)
            h = History(chromosomes=2, size=12)
            for depth in range(1, 31):
                if depth in (10, 20):
                    h.event(WGD)
                    continue
                available = [i for i, c in enumerate(h.contigs) if c]
                if not available:
                    break
                ci = rng.choice(available)
                kind = rng.choice(['InternalDuplication', 'InternalDeletion', 'ChromDuplication',
                                   'ChromDeletion', 'InternalInversion', 'InvertedDuplication'])
                start = rng.randrange(len(h.contigs[ci]))
                end = rng.randrange(start + 1, len(h.contigs[ci]) + 1)
                if kind.startswith('Chrom'):
                    start, end = None, None
                h.event(kind, ci, start, end)
            with self.subTest(seed=seed):
                self.reduce(h)


class BackwardVersionsTest(unittest.TestCase):
    def kept(self, h):
        lengths = [h.lengths[c] for _hap in ('H1', 'H2') for c in h.lengths]
        tracking = backward_dependencies(h.records, lengths)
        rows = [r for i, r in enumerate(h.records) if i in tracking['kept']]
        verify_filtered(rows, h.lengths, observed_edges(h.profile()))
        return tracking['kept']

    def test_edits_on_both_extinct_copy_branches_disappear(self):
        h = History().event('ChromDuplication').event('ChromDuplication', 2)
        h.event('InternalDuplication', 2, 4, 6).event('ChromDeletion', 2)
        h.event('InternalDuplication', 3, 1, 3).event('ChromDeletion', 3)
        self.assertEqual(self.kept(h), set())

    def test_edit_after_a_template_copy_can_disappear(self):
        h = History().event('ChromDuplication').event('InternalDuplication', 0, 2, 5)
        h.event('ChromDeletion', 0)
        self.assertEqual(self.kept(h), {0, 2})

    def test_used_template_keeps_its_eventual_deletion(self):
        h = History().event('ChromDuplication').event('InternalDuplication', 2, 1, 3)
        h.event('ChromDuplication', 2).event('ChromDeletion', 2)
        self.assertEqual(self.kept(h), {0, 1, 2, 3})

    def test_edit_erased_in_both_wgd_copies_disappears(self):
        h = History().event('InternalDuplication', 0, 2, 5).event(WGD)
        h.event('ChromDeletion', 0).event('ChromDeletion', 2)
        self.assertEqual(self.kept(h), {1, 2, 3})

    def test_edit_in_a_surviving_wgd_copy_stays(self):
        h = History().event('InternalDuplication', 0, 2, 5).event(WGD)
        h.event('ChromDeletion', 0)
        self.assertEqual(self.kept(h), {0, 1, 2})

    def test_two_wgds_preserve_ids_after_an_omitted_birth(self):
        h = History().event('ChromDuplication').event('ChromDeletion', 2)
        h.event(WGD).event(WGD).event('InternalDuplication', 6, 2, 5)
        self.assertEqual(self.kept(h), {2, 3, 4})

    def test_translocation_between_extinct_branches_disappears(self):
        h = History(chromosomes=2).event('ChromDuplication', 0).event('ChromDuplication', 1)
        h.event('Translocation', 4, 3, cb=5, pos_b=5, inverted_b=True)
        h.event('ChromDeletion', 4).event('ChromDeletion', 5)
        self.assertEqual(self.kept(h), set())

    def test_translocation_export_keeps_source_creation_and_cleanup(self):
        h = History(chromosomes=2).event('ChromDuplication', 1)
        h.event('Translocation', 0, 3, cb=4, pos_b=5, inverted_b=True)
        h.event('ChromDeletion', 4)
        self.assertEqual(self.kept(h), {0, 1, 2})

    def test_germline_extinction_keeps_both_deletions(self):
        h = History().event('InternalDuplication', 0, 2, 5)
        h.event('ChromDeletion', 0).event('ChromDeletion', 1)
        self.assertEqual(self.kept(h), {1, 2})

    def test_original_slot_counts_are_checked(self):
        h = History().event('ChromDuplication')
        h.records[0]['contig_slots_after'] = '2'
        with self.assertRaisesRegex(ValueError, 'contig slot counts'):
            self.kept(h)

    def test_omitted_birth_is_absent_not_an_empty_existing_contig(self):
        h = History().event('ChromDuplication').event('ChromDeletion', 2)
        with self.assertRaisesRegex(ValueError, 'contig 2 does not exist'):
            verify_filtered(h.records[1:], h.lengths, observed_edges(h.profile()))

    def test_legacy_full_history_gets_stable_allocation_metadata(self):
        h = History().event('ChromDuplication').event('ChromDeletion', 2).event(WGD)
        for row in h.records:
            row.pop('contig_slots_before')
            row.pop('contig_slots_after')
        labels, retained, _, _ = reduce_sample(h.sample, h.records, h.profile())
        verify_retained_index(h.records, labels, retained, initial_contigs=2)
        verify_filtered(retained, h.lengths, observed_edges(h.profile()))
        self.assertEqual(retained[0]['contig_slots_before'], '3')
        retained[0]['contig_slots_before'] = '2'
        with self.assertRaisesRegex(ValueError, 'allocation differs'):
            verify_retained_index(h.records, labels, retained, initial_contigs=2)


class ArchivedExamplesTest(unittest.TestCase):
    def example(self, sample):
        folder = Path(__file__).parent / 'fixtures/hidden_events'
        profiles, lengths = read_profiles(folder / 'copynumbers.tsv')
        with (folder / 'events.tsv').open() as handle:
            rows = [r for r in csv.DictReader(handle, delimiter='\t') if r['sample_id'] == sample]
        labels, retained, _, _ = reduce_sample(sample, rows, profiles[sample])
        verify_retained_index(rows, labels, retained, initial_contigs=2 * len(lengths))
        replayed = verify_filtered(retained, lengths, observed_edges(profiles[sample]))
        return rows, {int(r['depth']): r for r in labels}, replayed

    def test_sample_650_erased_chr16_branches(self):
        rows, labels, _ = self.example('Sample_650')
        local = [int(r['depth']) for r in rows if r['event_type'] != WGD
                 and 'chr16[' in r['regions_gained'] + r['regions_lost']]
        self.assertEqual([d for d in local if not labels[d]['hidden']], [80, 234, 250, 283, 406])
        self.assertTrue(labels[309]['hidden'])
        self.assertTrue(labels[328]['hidden'])

    def test_sample_4479_visible_duplication_and_tail_boundaries_stay(self):
        _, labels, replayed = self.example('Sample_4479')
        self.assertTrue(all(not labels[d]['hidden'] for d in (16, 79, 85)))
        edges = replayed['chr12', 'H1']
        self.assertEqual([sum(v for p, v in edges.items() if p <= x)
                          for x in (107755915, 108187203, 117426945)], [4, 3, 2])


class ExactReplayAndCliTest(unittest.TestCase):
    def test_intermediate_delta_and_effective_metadata_are_not_checked(self):
        h = History().event('ChromDuplication')
        h.records[0].update(regions_gained='old intermediate log', effective_cn_intervals='not JSON')
        verify_filtered(h.records, h.lengths, observed_edges(h.profile()))

    def test_missing_copy_template_fails_final_profile(self):
        h = History().event('InternalDuplication', 0, 8, 10).event('ChromDuplication')
        with self.assertRaisesRegex(ValueError, 'profile mismatch'):
            verify_filtered(h.records[1:], h.lengths, observed_edges(h.profile()))

    def test_manually_dropped_visible_tail_still_fails_final_profile(self):
        h = History().event('ChromDuplication').event('TelomereDeletion', 0, 5, 10)
        with self.assertRaisesRegex(ValueError, 'profile mismatch'):
            verify_filtered(h.records[:1], h.lengths, observed_edges(h.profile()))

    def test_wrong_contig_is_detected_by_final_profile(self):
        h = History().event('ChromDuplication')
        h.records[0]['description'] = 'contig:1;length:10;'
        with self.assertRaisesRegex(ValueError, 'profile mismatch'):
            verify_filtered(h.records, h.lengths, observed_edges(h.profile()))

    def test_retained_order_is_not_repaired(self):
        h = History().event('ChromDuplication').event('TelomereDeletion', 0, 5, 10)
        with self.assertRaisesRegex(ValueError, 'increasing depth order'):
            verify_filtered(list(reversed(h.records)), h.lengths, observed_edges(h.profile()))

    def test_empty_history_must_match_both_haplotypes(self):
        with self.assertRaisesRegex(ValueError, 'profile mismatch'):
            verify_filtered([], {'chr1': 10}, observed_edges({'chr1': [(0, 10, 1, 0)]}))

    def test_annotation_mismatch_and_wgd_removal_fail(self):
        h = History().event(WGD)
        labels, retained, _, _ = reduce_sample(h.sample, h.records, h.profile())
        labels[0]['hidden'] = True
        with self.assertRaisesRegex(ValueError, 'WGD must remain'):
            verify_retained_index(h.records, labels, [])
        with self.assertRaisesRegex(ValueError, 'Retained rows disagree'):
            verify_retained_index(h.records, labels, retained)

    def write_table(self, path, rows, columns=None):
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns or list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    def inputs(self, folder, histories):
        events, cn = folder / 'events.tsv', folder / 'copynumbers.tsv'
        self.write_table(events, [r for h in histories for r in h.records])
        self.write_table(cn, [dict(sample_id=h.sample, chrom=chrom, start=a + 1, end=b, cn_a=x, cn_b=y)
                             for h in histories for chrom, rows in h.profile().items() for a, b, x, y in rows])
        return [sys.executable, str(Path(__file__).resolve().parents[1] / 'hidden_events.py'),
                '--events', str(events), '--copynumbers', str(cn), '--out', str(folder / 'filtered'), '--progress', '0']

    def test_skip_policy_cli_and_verify_only(self):
        for keep in (False, True):
            with self.subTest(keep_skip=keep), tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                h = History().event('Skip').event('ChromDuplication')
                skip = History(sample='skip_only').event('Skip')
                command = self.inputs(folder, [h, skip])
                run = subprocess.run(command + (['--keep-skip'] if keep else []), capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr)
                summary = json.loads((folder / 'filtered/summary.json').read_text())
                self.assertEqual(summary['counts']['hidden_events'], 0 if keep else 2)
                self.assertEqual(summary['verification']['method'], 'filtered_original_operations_v3')
                self.assertEqual(summary['verification']['filtered_replay_pass'], 2)
                self.assertEqual(summary['verification']['filtered_replay_fail'], 0)
                verified = subprocess.run(command + ['--verify-only'], capture_output=True, text=True)
                self.assertEqual(verified.returncode, 0, verified.stderr)

    def test_skip_with_recorded_cn_change_is_rejected_by_detector(self):
        h = History().event('Skip')
        h.records[0]['regions_gained'] = '[H1:chr1[0:10)]'
        with self.assertRaisesRegex(ValueError, 'Skip events must not change'):
            reduce_sample(h.sample, h.records, h.profile())

    def test_invalid_new_output_is_not_published(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            h = History().event('ChromDuplication')
            command = self.inputs(folder, [h])
            # Delta/profile agree on a whole-chrom gain, but the unchanged command
            # actually copies half. Only executing it exposes the invalid input.
            h.records[0].update(event_type='InternalDuplication', description='contig:0;length:10;start:0;end:5',
                                contig_slots_after=h.records[0]['contig_slots_before'])
            self.write_table(folder / 'events.tsv', h.records)
            run = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn('Exact filtered replay failed', run.stderr)
            self.assertFalse((folder / 'filtered/events_filtered.tsv').exists())
            self.assertFalse((folder / 'filtered/summary.json').exists())

    def test_verify_only_rejects_old_representation_and_replaces_stale_success(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            h = History().event('ChromDuplication').event('TelomereDeletion', 0, 5, 10)
            command = self.inputs(folder, [h])
            labels, retained, _, _ = reduce_sample(h.sample, h.records, h.profile())
            labels[1]['hidden'] = True
            output = folder / 'filtered'
            output.mkdir()
            self.write_table(output / 'hidden_events.tsv', labels)
            retained[0]['effective_cn_intervals'] = '[["chr1","H1",0,5,1]]'
            self.write_table(output / 'events_filtered.tsv', retained[:1])
            (output / 'summary.json').write_text(json.dumps({'verification': {'reduced_replay_pass': 1}}))
            run = subprocess.run(command + ['--verify-only'], capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            summary = json.loads((output / 'summary.json').read_text())
            self.assertEqual(summary['verification']['filtered_replay_fail'], 1)
            self.assertNotIn('reduced_replay_pass', summary['verification'])
            with (output / 'verification.tsv').open() as handle:
                checks = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(checks[0]['filtered_replay'], 'fail')
            self.assertIn('profile mismatch', checks[0]['error'])


if __name__ == '__main__':
    unittest.main()

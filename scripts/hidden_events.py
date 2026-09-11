#!/usr/bin/env python
"""Annotate hidden CN events, emit a deterministic reduction, and verify its output.

Run once per cohort. WGD is fixed in the source history, but CN contributions cancel
across it. Acceptance requires exact replay of the retained original operations in
depth order. Effective CN intervals are reduction metadata and are never replay input.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from itertools import groupby

if __package__:
    from .lib_hidden_events import EVENT_KEY, VERSION, WGD, event_samples, observed_edges, read_profiles, reduce_sample
    from .hidden_event_replay import verify_filtered
else:
    from lib_hidden_events import EVENT_KEY, VERSION, WGD, event_samples, observed_edges, read_profiles, reduce_sample
    from hidden_event_replay import verify_filtered

ANNOTATIONS = [*EVENT_KEY, 'event_id', 'hidden', 'reason', 'retention_method',
               'input_version_ids', 'output_version_ids', 'source_units', 'retained_flow_units',
               'partially_cancelled', 'flow_choice_ambiguous', 'cancellation_group_id',
               'cancellation_group_members', 'dependency_event_ids', 'signal_ids', 'representative_signal_ids']
CONTRIBUTIONS = ['sample_id', 'event_id', 'interval_id', 'chrom', 'haplotype', 'start', 'end', 'sign',
                 'wgd_multiplier', 'retained_units', 'cancelled_units', 'flow_choice_ambiguous']
SIGNALS = ['signal_id', 'sample_id', 'chrom', 'haplotype', 'start', 'end', 'signed_units',
           'representative_event_id', 'support_event_ids', 'path_interval_ids', 'attribution_ambiguous']
REPLAY_METHOD = 'filtered_original_operations_v3'
REPLAY_SPACE = 'original_operations_with_stable_contig_ids'
VERIFICATION = ['sample_id', 'original_events', 'retained_events',
                'filtered_replay', 'profile_units', 'error']


class FilteredReplayError(ValueError):
    def __init__(self, verification, first_error):
        self.verification = verification
        super().__init__(f"Exact filtered replay failed for {verification['filtered_replay_fail']}/"
                         f"{verification['samples']} samples; first failure: {first_error}")


def digest(path):
    result = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def table_writer(folder, name, columns, handles):
    handle = (folder / name).open('w', newline='')
    handles.append(handle)
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    return writer


def filtered_samples(path):
    # Filtered output can have zero retained rows for a sample. Grouping is handled
    # here as a lookup so those samples are still verified against the full CN table.
    result = defaultdict(list)
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            result[row['sample_id']].append(row)
    return result



def annotation_samples(path):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for sample, rows in groupby(reader, key=lambda row: row['sample_id']):
            yield sample, list(rows)


def verify_retained_index(records, labels, retained, initial_contigs=None):
    original_keys = [tuple(row[key] for key in EVENT_KEY) for row in records]
    annotation_keys = [tuple(row[key] for key in EVENT_KEY) for row in labels]
    if annotation_keys != original_keys:
        raise ValueError('Annotations do not match the original event sequence')
    ids = [row['event_id'] for row in labels]
    if len(ids) != len(set(ids)) or any(str(row['hidden']) not in {'True', 'False'} for row in labels):
        raise ValueError('Invalid annotation IDs or hidden flags')
    expected = [row['event_id'] for row in labels if str(row['hidden']) == 'False']
    actual = [row['event_id'] for row in retained]
    if len(actual) != len(set(actual)) or actual != expected:
        raise ValueError('Retained rows disagree with hidden-event annotations or original order')
    if any(row['event_type'] == WGD and str(row['hidden']) != 'False' for row in labels):
        raise ValueError('WGD must remain at its original depth')
    source = {key: row for key, row in zip(original_keys, records)}
    allocation = {}
    count = initial_contigs
    if count is None and records and 'contig_slots_before' in records[0]:
        count = int(records[0]['contig_slots_before'])
    if count is not None:
        for original in records:
            kind = original['event_type']
            after = 2 * count if kind == WGD else count + (kind in {'ChromDuplication', 'ContigDuplication', 'ArmDuplication'})
            allocation[tuple(original[k] for k in EVENT_KEY)] = (count, after)
            count = after
    for row in retained:
        original = source[tuple(row[key] for key in EVENT_KEY)]
        if any(row[key] != value for key, value in original.items()):
            raise ValueError('Retained source metadata differs from original input')
        if allocation:
            try:
                actual_slots = tuple(int(row[k]) for k in ('contig_slots_before', 'contig_slots_after'))
            except (KeyError, ValueError, TypeError) as exc:
                raise ValueError('Missing or invalid retained contig slot metadata') from exc
            if actual_slots != allocation[tuple(row[k] for k in EVENT_KEY)]:
                raise ValueError('Retained contig allocation differs from original history')


def verify_outputs(events_path, cn_path, output_dir, *, sample_limit=None, progress=250):
    """Validate saved retained original operations; record every sample and reject failures."""
    profiles, lengths = read_profiles(cn_path)
    filtered = filtered_samples(output_dir / 'events_filtered.tsv')
    seen = set()
    checks = []
    originals = event_samples(events_path)
    annotation_groups = annotation_samples(output_dir / 'hidden_events.tsv')

    def check_sample(sample, records, labels):
        retained = filtered.get(sample, [])
        expected = observed_edges(profiles[sample])
        error = ''
        try:
            verify_retained_index(records, labels, retained, initial_contigs=2 * len(lengths))
            verify_filtered(retained, lengths, expected)
        except ValueError as exc:
            error = f'{sample}: {exc}'
        checks.append(dict(sample_id=sample, original_events=len(records),
                           retained_events=len(retained), filtered_replay='fail' if error else 'pass',
                           profile_units=len(expected), error=error))

    for index, (sample, records) in enumerate(originals):
        if sample_limit is not None and index >= sample_limit:
            break
        if sample not in profiles:
            raise ValueError(f'{sample}: events have no observed CN profile')
        seen.add(sample)
        annotation_sample, labels = next(annotation_groups, (None, []))
        if annotation_sample != sample:
            raise ValueError('Annotation sample groups disagree with original input')
        check_sample(sample, records, labels)
        if progress and (index + 1) % progress == 0:
            failed = sum(row['filtered_replay'] == 'fail' for row in checks)
            print(f'checked {index + 1} samples: exact filtered replay, {failed} failed', flush=True)
    if sample_limit is None:
        for sample in sorted(set(profiles) - seen):
            check_sample(sample, [], [])
            seen.add(sample)
    if next(annotation_groups, None) is not None:
        raise ValueError('Annotations contain unexpected extra samples')
    if set(filtered) - seen:
        raise ValueError('Filtered output contains unexpected samples')
    with (output_dir / 'verification.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=VERIFICATION, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(checks)
    failures = [row for row in checks if row['filtered_replay'] == 'fail']
    verification = {'method': REPLAY_METHOD, 'samples': len(checks),
                    'filtered_replay_pass': len(checks) - len(failures),
                    'filtered_replay_fail': len(failures),
                    'profile_units': sum(row['profile_units'] for row in checks)}
    if failures:
        raise FilteredReplayError(verification, failures[0]['error'])
    return verification


def update_verification_summary(output_dir, verification):
    """Replace historical replay claims when rechecking already published output."""
    path = output_dir / 'summary.json'
    if path.exists():
        summary = json.loads(path.read_text())
        summary['verification'] = verification
        summary['replay_space'] = REPLAY_SPACE
        summary.pop('source_event_epochs_preserved_as_metadata', None)
        summary['source_event_epochs_preserved'] = True
        temporary = path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(summary, indent=2) + '\n')
        os.replace(temporary, path)


def run(args):
    start = time.monotonic()
    if args.verify_only:
        try:
            verification = verify_outputs(args.events, args.copynumbers, args.out,
                                          sample_limit=args.samples, progress=args.progress)
        except FilteredReplayError as exc:
            update_verification_summary(args.out, exc.verification)
            print(json.dumps(exc.verification, indent=2), flush=True)
            raise
        update_verification_summary(args.out, verification)
        print(json.dumps(verification, indent=2))
        return
    profiles, _lengths = read_profiles(args.copynumbers)
    with args.events.open(newline='') as handle:
        raw_columns = next(csv.reader(handle, delimiter='\t'))
    extra = ['event_id', 'reduction_space', 'effective_cn_intervals', 'flow_choice_ambiguous']
    if set(extra) & set(raw_columns):
        raise ValueError('Input is already annotated; provide original SimChA events.tsv')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f'.{args.out.name}.', dir=args.out.parent))
    handles = []
    totals, by_type = Counter(), defaultdict(Counter)
    seen = set()
    try:
        annotations = table_writer(staging, 'hidden_events.tsv', ANNOTATIONS, handles)
        slot_columns = [k for k in ('contig_slots_before', 'contig_slots_after') if k not in raw_columns]
        filtered = table_writer(staging, 'events_filtered.tsv', raw_columns + slot_columns + extra, handles)
        contributions = table_writer(staging, 'cn_contributions.tsv', CONTRIBUTIONS, handles)
        signals = table_writer(staging, 'reduced_signals.tsv', SIGNALS, handles)
        for index, (sample, records) in enumerate(event_samples(args.events)):
            if args.samples is not None and index >= args.samples:
                break
            if sample not in profiles:
                raise ValueError(f'{sample}: events have no observed CN profile')
            seen.add(sample)
            labels, retained, pieces, reduced = reduce_sample(sample, records, profiles[sample], keep_skip=args.keep_skip)
            annotations.writerows(labels)
            filtered.writerows(retained)
            contributions.writerows(pieces)
            signals.writerows(reduced)
            totals['samples'] += 1
            if labels:
                totals[labels[0]['retention_method'] + '_samples'] += 1
            totals['events_including_wgd'] += len(records)
            totals['signal_rows'] += len(reduced)
            totals['signal_units'] += sum(abs(row['signed_units']) for row in reduced)
            for label in labels:
                if label['event_type'] == WGD:
                    totals['wgd_events'] += 1
                    continue
                counts = by_type[label['event_type']]
                for counter in (totals, counts):
                    counter['non_wgd_events'] += 1
                    counter['hidden_events'] += label['hidden']
                    counter[label['reason']] += 1
                    counter['partially_cancelled_events'] += label['partially_cancelled']
                    counter['flow_choice_ambiguous_events'] += label['flow_choice_ambiguous']
                    counter['cn_changing_events'] += bool(label['source_units'])
                    counter['hidden_cn_changing_events'] += bool(label['source_units']) and label['hidden']
            if args.progress and (index + 1) % args.progress == 0:
                print(f'reduced {index + 1} samples: {totals["hidden_events"]:,}/{totals["non_wgd_events"]:,} hidden', flush=True)
        if args.samples is None:
            for sample in sorted(set(profiles) - seen):
                reduce_sample(sample, [], profiles[sample], keep_skip=args.keep_skip)
                totals['samples'] += 1
        for handle in handles:
            handle.close()
        handles.clear()
        verification = verify_outputs(args.events, args.copynumbers, staging,
                                      sample_limit=args.samples, progress=args.progress)
        type_rows = []
        for kind, counts in sorted(by_type.items()):
            type_rows.append({'event_type': kind, 'events': counts['non_wgd_events'],
                              'hidden': counts['hidden_events'],
                              'hidden_fraction': counts['hidden_events'] / counts['non_wgd_events'],
                              'fully_cancelled': counts['fully_cancelled'],
                              'erased_contig_history': counts['erased_contig_history'],
                              'required_final_state': counts['required_final_state'],
                              'surviving_contribution': counts['surviving_contribution'],
                              'cancellation_group_kept': counts['cancellation_group_kept'],
                              'required_dependency': counts['required_dependency'],
                              'no_direct_cn_change': counts['no_direct_cn_change'],
                              'skip_kept': counts['skip_kept'],
                              'partially_cancelled': counts['partially_cancelled_events'],
                              'flow_choice_ambiguous': counts['flow_choice_ambiguous_events']})
        with (staging / 'summary_by_event_type.tsv').open('w', newline='') as handle:
            columns = ['event_type', 'events', 'hidden', 'hidden_fraction', 'fully_cancelled',
                       'erased_contig_history', 'required_final_state',
                       'surviving_contribution', 'cancellation_group_kept', 'required_dependency',
                       'no_direct_cn_change', 'skip_kept', 'partially_cancelled', 'flow_choice_ambiguous']
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(type_rows)
        summary = dict(algorithm=VERSION, dataset=args.dataset or args.events.parent.name,
                       definition='Erased contig histories or certified whole-event cancellation; WGD excluded from denominator',
                       cross_wgd_cancellation=True, whole_event_cancellation=True, dependency_protection=True,
                       backward_dependencies=True, stable_contig_ids=True,
                       candidate_selection="fewer retained events; backward tracking wins ties",
                       keep_skip=args.keep_skip,
                       source_event_epochs_preserved=True, replay_space=REPLAY_SPACE,
                       counts=dict(totals), hidden_fraction=totals['hidden_events'] / max(1, totals['non_wgd_events']),
                       hidden_cn_changing_fraction=totals['hidden_cn_changing_events'] / max(1, totals['cn_changing_events']),
                       verification=verification, input_events=str(args.events.resolve()),
                       input_copynumbers=str(args.copynumbers.resolve()),
                       input_events_sha256=digest(args.events), input_copynumbers_sha256=digest(args.copynumbers),
                       sample_limit=args.samples, elapsed_seconds=time.monotonic() - start)
        (staging / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        args.out.mkdir(parents=True, exist_ok=True)
        # Publish success metadata last. No output becomes current before all replay checks pass.
        for path in sorted(staging.iterdir(), key=lambda path: path.name == 'summary.json'):
            os.replace(path, args.out / path.name)
        print(json.dumps(summary, indent=2), flush=True)
    finally:
        for handle in handles:
            handle.close()
        shutil.rmtree(staging)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--events', '-i', type=Path, required=True)
    parser.add_argument('--copynumbers', type=Path, required=True)
    parser.add_argument('--out', '-o', type=Path, required=True)
    parser.add_argument('--dataset', default=None)
    parser.add_argument('--keep-skip', action='store_true',
                        help='Retain Skip rows without CN signals (default: remove them)')
    parser.add_argument('--samples', type=int, default=None, help='Limit samples for development; default all')
    parser.add_argument('--progress', type=int, default=250)
    parser.add_argument('--verify-only', action='store_true', help='Independently replay already written filtered output')
    args = parser.parse_args()
    if args.samples is not None and args.samples < 1:
        parser.error('--samples must be positive')
    run(args)


if __name__ == '__main__':
    main()

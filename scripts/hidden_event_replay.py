"""Replay retained operations at stable contig IDs and verify final CN."""
from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, replace

if __package__:
    from .lib_hidden_events import WGD, clean
else:
    from lib_hidden_events import WGD, clean

DESC = re.compile(r'contig:(\d+);length:(\d+)(?:;start:(\d+);end:(\d+))?;?')
PAIR = re.compile(r'contigA:(\d+);contigB:(\d+);posA:(\d+);posB:(\d+);invertedB:(True|False)', re.I)
APPEND = {'ChromDuplication', 'ContigDuplication', 'ArmDuplication'}
TAIL = {'TelomereDuplication', 'TailDuplication', 'ArmDuplication',
        'TelomereDeletion', 'TailDeletion', 'ArmDeletion'}
SLOT_FIELDS = ('contig_slots_before', 'contig_slots_after')


@dataclass(frozen=True)
class Block:
    chrom: str
    hap: str
    start: int
    end: int
    reverse: bool = False

    @property
    def length(self):
        return self.end - self.start


def chromosome_order(chrom):
    value = chrom.removeprefix('chr')
    return (0, int(value)) if value.isdigit() else (1, value)


def edges(blocks):
    result = defaultdict(Counter)
    for block in blocks:
        result[block.chrom, block.hap][block.start] += 1
        result[block.chrom, block.hap][block.end] -= 1
    return {key: clean(value) for key, value in result.items()}


def assert_profiles(actual, expected, context):
    for unit in set(actual) | set(expected):
        if clean(actual.get(unit, {})) != clean(expected.get(unit, {})):
            raise ValueError(f'{context}: profile mismatch on {unit}: '
                             f'{actual.get(unit, {})} != {expected.get(unit, {})}')


class MaterialReplay:
    def __init__(self, lengths):
        self.contigs = [[Block(chrom, hap, 0, lengths[chrom])]
                        for hap in ('H1', 'H2') for chrom in sorted(lengths, key=chromosome_order)]

    def blocks(self):
        return (block for contig in self.contigs if contig is not None for block in contig)

    def profile(self):
        return edges(self.blocks())

    def double(self):
        self.contigs.extend([None if contig is None else list(contig) for contig in self.contigs])

    def get(self, ci):
        if not 0 <= ci < len(self.contigs) or self.contigs[ci] is None:
            raise ValueError(f'contig {ci} does not exist')
        return self.contigs[ci]

    def prepare_slots(self, row):
        kind = row['event_type']
        if any(key in row for key in SLOT_FIELDS):
            try:
                before, after = (int(row[key]) for key in SLOT_FIELDS)
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError('Invalid contig slot counts') from exc
        else:
            before = len(self.contigs)
            after = 2 * before if kind == WGD else before + (kind in APPEND)
        expected = 2 * before if kind == WGD else before + (kind in APPEND)
        if before < len(self.contigs) or after != expected:
            raise ValueError(f'Invalid contig slot allocation: {before} -> {after}')
        # Missing births reserve IDs without creating material.
        self.contigs.extend([None] * (before - len(self.contigs)))
        return after

    def cut(self, contig_id, position):
        blocks = self.get(contig_id)
        offset = 0
        for index, block in enumerate(blocks):
            if position == offset:
                return index
            if offset < position < offset + block.length:
                amount = position - offset
                if block.reverse:
                    left = replace(block, start=block.end - amount)
                    right = replace(block, end=block.end - amount)
                else:
                    left = replace(block, end=block.start + amount)
                    right = replace(block, start=block.start + amount)
                blocks[index:index + 1] = [left, right]
                return index + 1
            offset += block.length
        if offset == position:
            return len(blocks)
        raise ValueError(f'Invalid contig position {contig_id}:{position}, length {offset}')

    def interval(self, ci, start, end):
        if not 0 <= start < end <= sum(block.length for block in self.get(ci)):
            raise ValueError('Invalid contig-space interval')
        left = self.cut(ci, start)
        right = self.cut(ci, end)
        return left, right, list(self.contigs[ci][left:right])

    def apply_original(self, row):
        after = self.prepare_slots(row)
        result = self.apply_operation(row)
        if len(self.contigs) != after:
            raise ValueError('Operation disagrees with logged contig allocation')
        return result

    def apply_operation(self, row):
        kind = row['event_type']
        if kind == WGD:
            added = list(self.blocks())
            self.double()
            return added, []
        if kind in {'Skip', 'Pass'}:
            return [], []
        if kind == 'Translocation':
            match = PAIR.fullmatch(row['description'])
            if not match:
                raise ValueError('Unsupported translocation description')
            a, b, pa, pb = map(int, match.groups()[:4])
            self.get(a)
            self.get(b)
            if a == b:
                raise ValueError('Translocation requires different contigs')
            if match.group(5).lower() == 'true':
                self.contigs[b] = [replace(block, reverse=not block.reverse)
                                   for block in reversed(self.contigs[b])]
            ia, ib = self.cut(a, pa), self.cut(b, pb)
            tail_a, tail_b = self.contigs[a][ia:], self.contigs[b][ib:]
            self.contigs[a][ia:], self.contigs[b][ib:] = tail_b, tail_a
            return [], []
        match = DESC.fullmatch(row['description'])
        if not match:
            raise ValueError(f'Original replay does not support {kind}: {row["description"]}')
        ci, recorded_length = int(match.group(1)), int(match.group(2))
        blocks = self.get(ci)
        if kind in {'ChromDeletion', 'ContigDeletion'}:
            touched = list(blocks)
            self.contigs[ci] = []
            return [], touched
        if kind in {'ChromDuplication', 'ContigDuplication'}:
            touched = list(blocks)
            self.contigs.append(touched)
            return touched, []
        if match.group(3) is None:
            raise ValueError(f'Missing coordinates for {kind}')
        start, end = int(match.group(3)), int(match.group(4))
        actual_length = sum(block.length for block in blocks)
        if kind in TAIL:
            if start:
                if end != recorded_length:
                    raise ValueError('Tail does not reach a recorded contig end')
                end = actual_length
            elif end == recorded_length and actual_length != recorded_length:
                raise ValueError('Ambiguous whole-tail direction in legacy description')
        left, right, touched = self.interval(ci, start, end)
        if kind == 'ArmDuplication':
            self.contigs.append(touched)
            return touched, []
        if kind in {'InternalDuplication', 'TelomereDuplication', 'TailDuplication',
                    'CentromereBoundDuplication', 'InvertedDuplication'}:
            copied = touched
            if kind == 'InvertedDuplication':
                copied = [replace(block, reverse=not block.reverse) for block in reversed(touched)]
            self.contigs[ci][right:right] = copied
            return copied, []
        if kind in {'InternalDeletion', 'TelomereDeletion', 'TailDeletion',
                    'ArmDeletion', 'CentromereBoundDeletion'}:
            del self.contigs[ci][left:right]
            return [], touched
        if kind == 'InternalInversion':
            self.contigs[ci][left:right] = [replace(block, reverse=not block.reverse)
                                           for block in reversed(touched)]
            return [], []
        raise ValueError(f'Original replay does not support event type {kind}')


def verify_filtered(records, lengths, expected):
    """Execute retained commands once; compare only their final material coverage."""
    replay = MaterialReplay(lengths)
    previous_depth = None
    for row in records:
        depth = int(row['depth'])
        context = f'{row["sample_id"]} d{depth} {row["event_type"]}'
        if previous_depth is not None and depth <= previous_depth:
            raise ValueError(f'{context}: filtered events are not in increasing depth order')
        previous_depth = depth
        if not all(key in row for key in SLOT_FIELDS):
            raise ValueError(f'{context}: filtered replay requires contig slot metadata; rerun the filter')
        try:
            replay.apply_original(row)
        except (ValueError, IndexError) as exc:
            raise ValueError(f'{context}: exact filtered replay failed: {exc}') from exc
    profile = replay.profile()
    assert_profiles(profile, expected, 'Exact filtered replay')
    return profile

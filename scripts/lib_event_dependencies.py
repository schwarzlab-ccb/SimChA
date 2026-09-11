"""Conservative dependencies of original contig commands, without material replay.

A retained command needs the prior versions of the contigs it reads. Appended
contig IDs also need the preceding allocation history. WGD copies version pointers;
it does not itself require every earlier edit to remain. This allows whole
cancellation groups to cross WGD when no surviving command reads their changes.
"""
from __future__ import annotations

import re

SINGLE = re.compile(r'contig:(\d+);length:(\d+)(?:;start:(\d+);end:(\d+))?;?')
PAIR = re.compile(r'contigA:(\d+);contigB:(\d+);posA:(\d+);posB:(\d+);invertedB:(True|False)', re.I)
APPEND = frozenset({'ChromDuplication', 'ContigDuplication', 'ArmDuplication'})
EDIT = frozenset({'InternalDuplication', 'TelomereDuplication', 'TailDuplication',
                  'CentromereBoundDuplication', 'InvertedDuplication', 'ChromDeletion',
                  'ContigDeletion', 'InternalDeletion', 'TelomereDeletion', 'TailDeletion',
                  'ArmDeletion', 'CentromereBoundDeletion', 'InternalInversion'})


def event_dependencies(records, initial_contigs):
    """Return earlier event indices required by each original command.

    `state[ci]` identifies the last operation that determines contig ci's material.
    `address[ci]` identifies the allocation prefix needed to preserve its numeric ID.
    Following these backward links retains original states conservatively, even
    when an intervening neutral group could in principle have restored a state.
    """
    dependencies = [set() for _ in records]
    state = [None] * initial_contigs
    address = [None] * initial_contigs
    last_append = None

    def read_contig(owner, ci):
        if not 0 <= ci < len(state):
            raise ValueError(f'd{records[owner]["depth"]}: original contig {ci} does not exist')
        dependencies[owner].update(x for x in (state[ci], address[ci]) if x is not None)

    for owner, row in enumerate(records):
        kind = row['event_type']
        if kind in {'Skip', 'Pass'}:
            continue
        if kind == 'WholeGenomeDoubling':
            # These copies start at the original pre-WGD contig count, which depends
            # on all preceding appends. WGD is always retained. Material dependencies
            # are inherited individually so unrelated earlier edits are not forced in.
            address.extend([last_append] * len(state))
            state.extend(list(state))
            continue
        if kind == 'Translocation':
            match = PAIR.fullmatch(row['description'])
            if not match:
                raise ValueError(f'Unsupported translocation description: {row["description"]}')
            a, b = map(int, match.groups()[:2])
            read_contig(owner, a)
            read_contig(owner, b)
            state[a] = state[b] = owner
            continue
        if kind not in APPEND | EDIT:
            raise ValueError(f'Unsupported event for dependency protection: {kind}')
        match = SINGLE.fullmatch(row['description'])
        if not match:
            raise ValueError(f'Unsupported {kind} description: {row["description"]}')
        ci = int(match.group(1))
        read_contig(owner, ci)
        if kind in APPEND:
            # Retaining this append preserves its own slot only if earlier appends
            # survive. WGD timing is fixed, so this prefix also fixes doubled offsets.
            if last_append is not None:
                dependencies[owner].add(last_append)
            state.append(owner)
            address.append(owner)
            last_append = owner
        else:
            state[ci] = owner
    return dependencies


def protect_groups(groups, seeds, dependencies):
    """Close retained cancellation groups over all original-command dependencies."""
    owner_group = {owner: group for group, members in groups.items() for owner in members}
    kept = set(seeds)
    pending = sorted(kept, reverse=True)
    while pending:
        group = pending.pop()
        for owner in groups[group]:
            for prerequisite in sorted(dependencies[owner]):
                needed = owner_group[prerequisite]
                if needed not in kept:
                    kept.add(needed)
                    pending.append(needed)
    return kept

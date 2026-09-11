"""Backward dependencies of contig versions; no material replay."""
from collections import defaultdict
from dataclasses import dataclass

if __package__:
    from .lib_event_dependencies import SINGLE, PAIR, APPEND, EDIT
else:
    from lib_event_dependencies import SINGLE, PAIR, APPEND, EDIT

WGD = 'WholeGenomeDoubling'
CLEAR = {'ChromDeletion', 'ContigDeletion'}
NOOP = {'Skip', 'Pass'}
SLOT_FIELDS = ('contig_slots_before', 'contig_slots_after')


def slot_counts(row, before, after):
    present = [key in row for key in SLOT_FIELDS]
    if any(present):
        if not all(present):
            raise ValueError('Both contig slot counts are required')
        try:
            logged = tuple(int(row[key]) for key in SLOT_FIELDS)
        except (TypeError, ValueError) as exc:
            raise ValueError('Invalid contig slot counts') from exc
        if logged != (before, after):
            raise ValueError(f'd{row["depth"]}: contig slot counts {logged} != {(before, after)}')
    return before, after


@dataclass(frozen=True)
class Version:
    owner: int | None
    contig: int
    parents: tuple
    length: int


def backward_dependencies(records, initial_lengths, *, keep_skip=False):
    nodes = [Version(None, ci, (), length) for ci, length in enumerate(initial_lengths)]
    state = list(range(len(nodes)))
    birth = [None] * len(state)
    reads = [set() for _ in records]
    writes = [[] for _ in records]
    slots = []

    def read(ci):
        if not 0 <= ci < len(state):
            raise ValueError(f'Original contig {ci} does not exist')
        return state[ci]

    def version(owner, ci, parents, length):
        index = len(nodes)
        nodes.append(Version(owner, ci, tuple(parents), length))
        reads[owner].update(parents)
        writes[owner].append(index)
        return index

    for owner, row in enumerate(records):
        kind = row['event_type']
        before = len(state)
        if kind == WGD:
            for ci in range(before):
                parent = state[ci]
                state.append(version(owner, before + ci, (parent,), nodes[parent].length))
                birth.append(birth[ci])
        elif kind in NOOP:
            pass
        elif kind == 'Translocation':
            match = PAIR.fullmatch(row['description'])
            if not match:
                raise ValueError('Unsupported translocation description')
            a, b, pa, pb = map(int, match.groups()[:4])
            parents = (read(a), read(b))
            la, lb = (nodes[p].length for p in parents)
            if a == b or not (0 <= pa <= la and 0 <= pb <= lb):
                raise ValueError('Invalid translocation coordinates')
            state[a] = version(owner, a, parents, pa + lb - pb)
            state[b] = version(owner, b, parents, pb + la - pa)
        elif kind in APPEND | EDIT:
            match = SINGLE.fullmatch(row['description'])
            if not match:
                raise ValueError(f'Unsupported {kind} description')
            ci, logged_length = map(int, match.groups()[:2])
            parent = read(ci)
            length = nodes[parent].length
            if length != logged_length:
                raise ValueError(f'd{row["depth"]}: original length {logged_length} != {length}')
            start, end = ((int(match.group(3)), int(match.group(4)))
                          if match.group(3) is not None else (0, length))
            if not 0 <= start < end <= length:
                raise ValueError('Invalid original interval coordinates')
            if kind in CLEAR:
                state[ci] = version(owner, ci, (), 0)
            elif kind in APPEND:
                new = version(owner, before, (parent,), end - start)
                state.append(new)
                birth.append(new)
            else:
                size = end - start
                after_length = (length + size if kind.endswith('Duplication') else
                                length - size if kind.endswith('Deletion') else length)
                state[ci] = version(owner, ci, (parent,), after_length)
        else:
            raise ValueError(f'Unsupported event for backward tracking: {kind}')
        slots.append(slot_counts(row, before, len(state)))

    finals_by_birth = defaultdict(list)
    roots = []
    for ci, final in enumerate(state):
        if birth[ci] is not None:
            finals_by_birth[birth[ci]].append(final)
        if birth[ci] is None or nodes[final].length:
            roots.append(final)
    pending = roots.copy()
    demanded = set()
    dependencies = [set() for _ in records]
    kept = {i for i, r in enumerate(records)
            if r['event_type'] == WGD or (keep_skip and r['event_type'] == 'Skip')}
    while pending:
        index = pending.pop()
        if index in demanded:
            continue
        demanded.add(index)
        node = nodes[index]
        if node.owner is not None:
            kept.add(node.owner)
        prerequisites = list(node.parents)
        origin = birth[node.contig]
        if origin is not None and origin != index:
            prerequisites.append(origin)
        for parent in prerequisites:
            pending.append(parent)
            owner = nodes[parent].owner
            if node.owner is not None and owner is not None and owner != node.owner:
                dependencies[node.owner].add(owner)
        # A retained template must also reach its recorded final state.
        pending.extend(finals_by_birth.get(index, ()))
    return dict(kept=kept, dependencies=dependencies, slots=slots, nodes=nodes,
                reads=reads, writes=writes,
                final_owners={nodes[i].owner for i in roots if nodes[i].owner is not None})

"""Hidden-event reduction by backward dependencies and certified cancellation.

Backward contig dependencies identify erased histories. Certified whole-event
cancellation remains available; replay is used only for final verification.
"""
from __future__ import annotations

import bisect
import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import groupby

import networkx as nx

if __package__:
    from .lib_event_dependencies import event_dependencies, protect_groups
    from .lib_contig_versions import backward_dependencies
else:
    from lib_event_dependencies import event_dependencies, protect_groups
    from lib_contig_versions import backward_dependencies

VERSION = 'backward-contig-dependencies-v4'
EVENT_KEY = ('sample_id', 'depth', 'event_type', 'description')
REGION = re.compile(r'(H\d+):(chr[\w]+)\[(\d+):(\d+)\)')
WGD = 'WholeGenomeDoubling'


def clean(counter):
    return {key: value for key, value in counter.items() if value}


def parse_regions(value):
    """Parse the full cell; silently accepting malformed/missing deltas is unsafe."""
    if value == '[]':
        return []
    if not isinstance(value, str) or not value.startswith('[') or not value.endswith(']'):
        raise ValueError(f'Missing or malformed reference delta: {value!r}')
    result = []
    for token in value[1:-1].split(','):
        match = REGION.fullmatch(token)
        if match is None:
            raise ValueError(f'Malformed reference interval: {token!r}')
        hap, chrom, start, end = match.groups()
        start, end = int(start), int(end)
        if start >= end:
            raise ValueError(f'Empty or reversed reference interval: {token!r}')
        result.append((chrom, hap, start, end))
    return result


def event_samples(path):
    """Stream contiguous samples, rejecting repeated sample groups or duplicate keys."""
    seen = set()
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = set(EVENT_KEY) | {'regions_gained', 'regions_lost'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: required columns missing: {required - set(reader.fieldnames or [])}')
        for sample, rows in groupby(reader, key=lambda row: row['sample_id']):
            if sample in seen:
                raise ValueError(f'{path}: sample {sample} occurs in multiple groups')
            seen.add(sample)
            records = sorted(rows, key=lambda row: int(row['depth']))
            keys = [tuple(row[k] for k in EVENT_KEY) for row in records]
            if len(set(keys)) != len(keys):
                raise ValueError(f'{sample}: duplicate event key')
            depths = [int(row['depth']) for row in records]
            if len(set(depths)) != len(depths):
                raise ValueError(f'{sample}: repeated depth; input must contain one complete lineage per sample')
            yield sample, records


def read_profiles(path):
    profiles = defaultdict(lambda: defaultdict(list))
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'sample_id', 'chrom', 'start', 'end', 'cn_a', 'cn_b'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: required CN columns missing')
        for row in reader:
            start, end = int(row['start']) - 1, int(row['end'])
            a, b = int(row['cn_a']), int(row['cn_b'])
            if start < 0 or end <= start or min(a, b) < 0:
                raise ValueError(f'Invalid observed CN row: {row}')
            profiles[row['sample_id']][row['chrom']].append((start, end, a, b))
    lengths = {}
    for sample, chromosomes in profiles.items():
        for chrom, rows in chromosomes.items():
            rows.sort()
            previous = 0
            for start, end, _a, _b in rows:
                if start != previous:
                    raise ValueError(f'{sample} {chrom}: observed coverage is not contiguous from zero')
                previous = end
            if lengths.setdefault(chrom, previous) != previous:
                raise ValueError(f'Inconsistent chromosome length: {sample} {chrom}')
    if not profiles:
        raise ValueError('No observed CN profiles')
    return dict(profiles), lengths


def observed_edges(profile, baseline=0):
    result = defaultdict(Counter)
    for chrom, rows in profile.items():
        for start, end, a, b in rows:
            for hap, cn in (('H1', a), ('H2', b)):
                result[chrom, hap][start] += cn - baseline
                result[chrom, hap][end] -= cn - baseline
    return {key: clean(value) for key, value in result.items()}


@dataclass(frozen=True)
class Arc:
    index: int
    owner: int
    chrom: str
    hap: str
    start: int
    end: int
    sign: int
    multiplier: int

    @property
    def source(self):
        return self.start if self.sign > 0 else self.end

    @property
    def target(self):
        return self.end if self.sign > 0 else self.start


def variable_arcs(arcs, flows):
    """Flag arcs that can vary among unit-cost optima (not full event min/max).

    Residual shortest-path potentials expose zero-cost cycles. An arc changes in an
    alternative optimum only if its residual edge lies on such a cycle. Capacities
    are not expanded into units. Multi-arc event omission is a separate joint query.
    """
    residual = []
    nodes = set()
    for arc in arcs:
        value = flows[arc.index]
        u, v = arc.source, arc.target
        nodes.update((u, v))
        if value < arc.multiplier:
            residual.append((u, v, 1, arc.index))
        if value > 0:
            residual.append((v, u, -1, arc.index))
    distance = dict.fromkeys(nodes, 0)
    for _ in range(len(nodes)):
        changed = False
        for u, v, weight, _index in residual:
            if distance[v] > distance[u] + weight:
                distance[v] = distance[u] + weight
                changed = True
        if not changed:
            break
    else:
        raise AssertionError('Negative-cost cycle in purported optimal flow')
    zero = nx.DiGraph()
    zero.add_nodes_from(nodes)
    for u, v, weight, _index in residual:
        if weight + distance[u] - distance[v] == 0:
            zero.add_edge(u, v)
    component = {node: index for index, members in enumerate(nx.strongly_connected_components(zero))
                 for node in members}
    adjacency = defaultdict(list)
    for u, v, weight, index in residual:
        if weight + distance[u] - distance[v] == 0:
            adjacency[u].append((v, index))
    variable = set()
    for u, v, weight, index in residual:
        if weight + distance[u] - distance[v] or component[u] != component[v]:
            continue
        # A partially used arc has both residual directions. Its own forward/back
        # two-cycle changes nothing; require a return path through other arcs.
        pending, visited = [v], {v}
        while pending:
            node = pending.pop()
            if node == u:
                variable.add(index)
                break
            for target, other in adjacency[node]:
                if other != index and target not in visited:
                    visited.add(target)
                    pending.append(target)
    return variable


def solve_unit(arcs, expected):
    graph = nx.MultiDiGraph()
    balance = Counter()
    for arc in arcs:
        graph.add_edge(arc.source, arc.target, key=arc.index, capacity=arc.multiplier, weight=1)
        balance[arc.source] += arc.multiplier
        balance[arc.target] -= arc.multiplier
    if clean(balance) != expected:
        raise ValueError(f'Original deltas disagree with observed CN: {clean(balance)} != {expected}')
    if not arcs:
        return {}, [], set()
    nx.set_node_attributes(graph, {node: -balance[node] for node in graph}, 'demand')
    _cost, flow = nx.network_simplex(graph)
    flows = {arc.index: flow[arc.source][arc.target][arc.index] for arc in arcs}
    variable = variable_arcs(arcs, flows)
    remaining = flows.copy()
    outgoing = defaultdict(list)
    by_id = {arc.index: arc for arc in arcs}
    for arc in arcs:
        if remaining[arc.index]:
            outgoing[arc.source].append(arc.index)
    # Prefer shorter paths at a tie, then a stable source event/interval order.
    # The dictionary and all input arcs are deterministically ordered.
    for ids in outgoing.values():
        ids.sort(key=lambda index: (by_id[index].owner, index))
    supply = {node: value for node, value in balance.items() if value > 0}
    demand = {node: -value for node, value in balance.items() if value < 0}
    paths = []
    for source in sorted(supply):
        while supply[source]:
            node, path, visited = source, [], set()
            while not demand.get(node, 0):
                if node in visited:
                    raise AssertionError('Positive-cost optimum contains a circulation')
                visited.add(node)
                options = outgoing[node]
                while options and not remaining[options[0]]:
                    options.pop(0)
                if not options:
                    raise AssertionError('Flow path ends without demand')
                index = options[0]
                path.append(index)
                node = by_id[index].target
            amount = min(supply[source], demand[node], *(remaining[index] for index in path))
            supply[source] -= amount
            demand[node] -= amount
            for index in path:
                remaining[index] -= amount
            sign = 1 if source < node else -1
            candidates = [by_id[index] for index in path if by_id[index].sign == sign]
            # Assign the contracted signal to the widest same-direction contributor,
            # breaking ties by original event order. Other owners remain in support_ids.
            representative = min(candidates, key=lambda arc: (-(arc.end - arc.start), arc.owner, arc.index)).owner
            paths.append({'start': min(source, node), 'end': max(source, node), 'sign': sign,
                          'units': amount, 'owner': representative, 'path': path,
                          'support': sorted({by_id[index].owner for index in path})})
    if any(remaining.values()) or any(demand.values()):
        raise AssertionError('Incomplete flow decomposition')
    return flows, paths, variable


def cancellation_groups(arcs, flows, event_count):
    """Union owners of cancelled-flow cycles into atomic source-event groups.

    Capacities are processed in batches, never expanded into individual WGD units.
    Joining cycles that share an event ensures a removal never takes only part of
    that event's footprint. A group with any live units will be retained later.
    """
    parent = list(range(event_count))
    def find(owner):
        while parent[owner] != owner:
            parent[owner] = parent[parent[owner]]
            owner = parent[owner]
        return owner
    def join(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[max(a, b)] = min(a, b)

    by_id = {arc.index: arc for arc in arcs}
    remaining = {arc.index: arc.multiplier - flows[arc.index] for arc in arcs}
    outgoing = defaultdict(list)
    cursor = defaultdict(int)
    def source(arc):
        return arc.chrom, arc.hap, arc.source
    def target(arc):
        return arc.chrom, arc.hap, arc.target
    for arc in arcs:
        if remaining[arc.index]:
            outgoing[source(arc)].append(arc.index)
    for arc in arcs:
        while remaining[arc.index]:
            node = source(arc)
            positions, path = {}, []
            while node not in positions:
                positions[node] = len(path)
                options = outgoing[node]
                while cursor[node] < len(options) and not remaining[options[cursor[node]]]:
                    cursor[node] += 1
                if cursor[node] == len(options):
                    raise AssertionError('Cancelled flow is not a circulation')
                index = options[cursor[node]]
                path.append(index)
                node = target(by_id[index])
            cycle = path[positions[node]:]
            amount = min(remaining[index] for index in cycle)
            owner = by_id[cycle[0]].owner
            for index in cycle:
                remaining[index] -= amount
                join(owner, by_id[index].owner)
    groups = defaultdict(list)
    for owner in range(event_count):
        groups[find(owner)].append(owner)
    return dict(groups)


def reduce_sample(sample, records, profile, *, keep_skip=False):
    wgd_depths = [int(row['depth']) for row in records if row['event_type'] == WGD]
    baseline = 2 ** len(wgd_depths)
    expected = observed_edges(profile, baseline)
    arcs, units = [], defaultdict(list)
    source_count = Counter()
    for owner, row in enumerate(records):
        gained, lost = parse_regions(row['regions_gained']), parse_regions(row['regions_lost'])
        if row['event_type'] == WGD:
            continue
        if row['event_type'] in {'Skip', 'Pass'} and (gained or lost):
            raise ValueError(f'{sample}: Skip events must not change copy number (nor Pass)')
        multiplier = 2 ** (len(wgd_depths) - bisect.bisect_right(wgd_depths, int(row['depth'])))
        for sign, regions in ((1, gained), (-1, lost)):
            for chrom, hap, start, end in regions:
                if (chrom, hap) not in expected or end > profile[chrom][-1][1]:
                    raise ValueError(f'{sample}: delta outside observed reference: {(chrom, hap, start, end)}')
                arc = Arc(len(arcs), owner, chrom, hap, start, end, sign, multiplier)
                arcs.append(arc)
                units[chrom, hap].append(arc)
                source_count[owner] += multiplier
    flows, signals, variable = {}, [], set()
    for unit in sorted(expected):
        local_flow, paths, local_variable = solve_unit(units[unit], expected[unit])
        flows.update(local_flow)
        variable.update(local_variable)
        for path in paths:
            signals.append({'signal_id': f'{sample}:s{len(signals) + 1}', 'sample_id': sample,
                            'chrom': unit[0], 'haplotype': unit[1], **path})
    owned = defaultdict(list)
    support = defaultdict(set)
    for signal in signals:
        owned[signal['owner']].append(signal)
        for owner in signal['support']:
            support[owner].add(signal['signal_id'])
    live = Counter()
    varied = set()
    for arc in arcs:
        live[arc.owner] += flows[arc.index]
        if arc.index in variable:
            varied.add(arc.owner)
    groups = cancellation_groups(arcs, flows, len(records))
    owner_group = {owner: group for group, members in groups.items() for owner in members}
    seeds = {group for group, members in groups.items()
             if any(live[owner] or records[owner]['event_type'] == WGD
                    or (keep_skip and records[owner]['event_type'] == 'Skip') for owner in members)}
    # Every eligible group must cancel FULL original event contributions, not merely
    # cancelled pieces of a retained event. Check its signed certificate explicitly.
    certificates = defaultdict(Counter)
    for arc in arcs:
        group = owner_group[arc.owner]
        if group not in seeds:
            certificates[group][arc.chrom, arc.hap, arc.start] += arc.sign * arc.multiplier
            certificates[group][arc.chrom, arc.hap, arc.end] -= arc.sign * arc.multiplier
    if any(clean(value) for value in certificates.values()):
        raise AssertionError('Candidate cancellation group changes final CN')
    dependencies = event_dependencies(records, initial_contigs=2 * len(profile))
    kept_groups = protect_groups(groups, seeds, dependencies)
    flow_kept = {owner for group in kept_groups for owner in groups[group]}
    chromosomes = sorted(profile, key=lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c))
    lengths = [profile[c][-1][1] for _hap in ('H1', 'H2') for c in chromosomes]
    tracking = backward_dependencies(records, lengths, keep_skip=keep_skip)
    use_backward = len(tracking['kept']) <= len(flow_kept)
    kept = tracking['kept'] if use_backward else flow_kept
    method = 'backward_dependencies' if use_backward else 'whole_event_cancellation'
    if use_backward:
        dependencies = tracking['dependencies']

    def version_id(index):
        node = tracking['nodes'][index]
        epoch = 'germline' if node.owner is None else 'd' + str(records[node.owner]['depth'])
        return f'{sample}:c{node.contig}@{epoch}'

    original_intervals = defaultdict(list)
    for arc in arcs:
        original_intervals[arc.owner].append([arc.chrom, arc.hap, arc.start, arc.end,
                                              arc.sign * arc.multiplier])
    annotations, filtered = [], []
    for owner, row in enumerate(records):
        is_wgd = row['event_type'] == WGD
        skip_kept = keep_skip and row['event_type'] == 'Skip'
        group = owner_group[owner]
        hidden = owner not in kept
        reason = ('wgd_kept' if is_wgd else 'skip_kept' if skip_kept
                  else 'no_direct_cn_change' if hidden and not source_count[owner]
                  else 'erased_contig_history' if hidden and use_backward
                  else 'fully_cancelled' if hidden
                  else 'required_final_state' if use_backward and owner in tracking['final_owners']
                  else 'required_dependency' if use_backward
                  else 'surviving_contribution' if live[owner]
                  else 'cancellation_group_kept' if group in seeds
                  else 'required_dependency')
        annotation = {key: row[key] for key in EVENT_KEY}
        annotation.update(event_id=f'{sample}:e{owner + 1}', hidden=hidden, reason=reason,
                          retention_method=method,
                          input_version_ids=json.dumps([version_id(i) for i in sorted(tracking['reads'][owner])], separators=(',', ':')),
                          output_version_ids=json.dumps([version_id(i) for i in tracking['writes'][owner]], separators=(',', ':')),
                          source_units=source_count[owner], retained_flow_units=live[owner],
                          partially_cancelled=0 < live[owner] < source_count[owner],
                          flow_choice_ambiguous=owner in varied,
                          cancellation_group_id=f'{sample}:c{group + 1}',
                          cancellation_group_members=json.dumps([f'{sample}:e{i + 1}' for i in groups[group]], separators=(',', ':')),
                          dependency_event_ids=json.dumps([f'{sample}:e{i + 1}' for i in sorted(dependencies[owner])], separators=(',', ':')),
                          signal_ids=json.dumps(sorted(support[owner]), separators=(',', ':')),
                          representative_signal_ids=json.dumps([s['signal_id'] for s in owned[owner]], separators=(',', ':')))
        annotations.append(annotation)
        if not hidden:
            # Informational original source deltas, untrimmed and WGD-scaled. The
            # actual replay uses only the original commands, never these intervals.
            filtered.append({**row,
                             'contig_slots_before': str(tracking['slots'][owner][0]),
                             'contig_slots_after': str(tracking['slots'][owner][1]),
                             'event_id': annotation['event_id'],
                             'reduction_space': 'final_copy_number',
                             'effective_cn_intervals': json.dumps(original_intervals[owner], separators=(',', ':')),
                             'flow_choice_ambiguous': owner in varied})
    contribution_rows = []
    for arc in arcs:
        contribution_rows.append(dict(sample_id=sample, event_id=f'{sample}:e{arc.owner + 1}',
                                      interval_id=arc.index, chrom=arc.chrom, haplotype=arc.hap,
                                      start=arc.start, end=arc.end, sign=arc.sign,
                                      wgd_multiplier=arc.multiplier, retained_units=flows[arc.index],
                                      cancelled_units=arc.multiplier - flows[arc.index],
                                      flow_choice_ambiguous=arc.index in variable))
    signal_rows = []
    for signal in signals:
        signal_rows.append(dict(signal_id=signal['signal_id'], sample_id=sample,
                                chrom=signal['chrom'], haplotype=signal['haplotype'],
                                start=signal['start'], end=signal['end'],
                                signed_units=signal['sign'] * signal['units'],
                                representative_event_id=f'{sample}:e{signal["owner"] + 1}',
                                support_event_ids=json.dumps([f'{sample}:e{i + 1}' for i in signal['support']], separators=(',', ':')),
                                path_interval_ids=json.dumps(signal['path'], separators=(',', ':')),
                                attribution_ambiguous=len(signal['support']) > 1 or any(i in variable for i in signal['path'])))
    return annotations, filtered, contribution_rows, signal_rows

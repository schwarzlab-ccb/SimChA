"""Library for detecting events hidden by later same-contig events.

The ground-truth event file stores affected reference intervals in strings such
as ``[H1:chr1[10:20),H1:chr1[30:40)]``.  A physical deletion can be split into
several such intervals when material between them was already lost.  For that
reason, containment is evaluated using the outer span of all fragments from the
same event, allele, chromosome, and contig. An earlier deletion is also hidden
when a later deletion directly adjoins it on the same contig, because the two
losses form one uninterrupted final loss. The equivalent check is applied to
two directly adjoining gains. A contained event is marked as rescued when an
intermediate duplication overlaps it before the containing deletion occurs.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass

import pandas as pd


_REGION_PATTERN = re.compile(
    r"(?P<allele>[^:]+):(?P<chrom>chr[^\[]+)"
    r"\[(?P<start>\d+):(?P<end>\d+)\)"
)
_CONTIG_PATTERN = re.compile(r"(?:^|;)contig:(?P<contig>\d+)(?:;|$)")

_REQUIRED_COLUMNS = {
    "sample_id",
    "event_type",
    "depth",
    "description",
    "regions_gained",
    "regions_lost",
}


@dataclass(frozen=True)
class Region:
    """A zero-based, half-open reference interval on one allele."""

    allele: str
    chrom: str
    start: int
    end: int


def parse_regions(value: str) -> tuple[Region, ...]:
    """Parse a GT_events.tsv region list and reject malformed input."""

    value = value.strip()
    if value == "[]":
        return ()
    if not (value.startswith("[") and value.endswith("]")):
        raise ValueError(f"Malformed region list: {value!r}")

    regions = []
    for item in value[1:-1].split(","):
        match = _REGION_PATTERN.fullmatch(item.strip())
        if match is None:
            raise ValueError(f"Malformed region: {item!r}")
        start = int(match.group("start"))
        end = int(match.group("end"))
        if end <= start:
            raise ValueError(f"Region must have positive width: {item!r}")
        regions.append(
            Region(
                allele=match.group("allele"),
                chrom=match.group("chrom"),
                start=start,
                end=end,
            )
        )
    return tuple(regions)


def parse_contig_id(description: str) -> int:
    """Extract the simulated contig identifier from an event description."""

    match = _CONTIG_PATTERN.search(description)
    if match is None:
        raise ValueError(f"Event description has no contig identifier: {description!r}")
    return int(match.group("contig"))


def build_event_footprints(
    events: pd.DataFrame,
    *,
    excluded_event_types: Iterable[str] = ("WholeGenomeDoubling",),
) -> pd.DataFrame:
    """Reduce each event to one outer span per allele and chromosome.

    Repeated regions are retained in ``region_count`` because they describe
    copy multiplicity, while ``unique_region_count`` ignores such repetition.
    Whole-genome doubling is excluded by default because it is not a local
    chromosome event that can be wholly hidden by one later deletion.
    """

    missing = _REQUIRED_COLUMNS.difference(events.columns)
    if missing:
        missing_list = ", ".join(sorted(missing))
        raise ValueError(f"Missing required event columns: {missing_list}")

    excluded = set(excluded_event_types)
    footprints: list[dict[str, object]] = []

    columns = [
        "sample_id",
        "event_type",
        "depth",
        "description",
        "regions_gained",
        "regions_lost",
    ]
    for row in events.loc[:, columns].itertuples(index=False):
        if row.event_type in excluded:
            continue

        regions = parse_regions(row.regions_gained) + parse_regions(row.regions_lost)
        if not regions:
            continue
        contig = parse_contig_id(row.description)

        grouped_regions: dict[tuple[str, str], list[Region]] = {}
        for region in regions:
            grouped_regions.setdefault((region.allele, region.chrom), []).append(region)

        for (allele, chrom), event_regions in grouped_regions.items():
            start = min(region.start for region in event_regions)
            end = max(region.end for region in event_regions)
            footprints.append(
                {
                    "sample_id": row.sample_id,
                    "chrom": chrom,
                    "allele": allele,
                    "contig": contig,
                    "depth": int(row.depth),
                    "event_type": row.event_type,
                    "start": start,
                    "end": end,
                    "width": end - start,
                    "region_count": len(event_regions),
                    "unique_region_count": len(
                        {(region.start, region.end) for region in event_regions}
                    ),
                    "intervals": tuple(
                        sorted(
                            {
                                (region.start, region.end)
                                for region in event_regions
                            }
                        )
                    ),
                    "description": row.description,
                }
            )

    return pd.DataFrame.from_records(
        footprints,
        columns=[
            "sample_id",
            "chrom",
            "allele",
            "contig",
            "depth",
            "event_type",
            "start",
            "end",
            "width",
            "region_count",
            "unique_region_count",
            "intervals",
            "description",
        ],
    )


def find_ducking_events(
    events: pd.DataFrame,
    *,
    include_wgd_crossing: bool = False,
) -> pd.DataFrame:
    """Find events hidden by containment or adjoining same-contig CN changes.

    Results contain one row per hidden event. If several later deletions
    qualify, ``hiding_*`` describes the earliest one in event history. For a
    contained event, ``rescuing_*`` describes the earliest duplication that
    overlaps the event after it occurs and before that deletion. Containment
    is inclusive at both outer boundaries, overlap must have positive width,
    and event order is taken from ``depth`` within each sample.

    By default, a deletion is not allowed to hide an event from before a
    whole-genome doubling.  The doubling copied that event, so loss of one
    descendant copy does not establish loss of the event as a whole.  Set
    ``include_wgd_crossing=True`` for a coordinate-only containment screen.
    """

    footprints = build_event_footprints(events)
    wgd_depths_by_sample = (
        events.loc[
            events["event_type"] == "WholeGenomeDoubling",
            ["sample_id", "depth"],
        ]
        .groupby("sample_id")["depth"]
        .agg(tuple)
        .to_dict()
    )
    ordered = footprints.sort_values(
        ["sample_id", "chrom", "allele", "contig", "depth"],
        ascending=[True, True, True, True, False],
        kind="stable",
    )

    results: list[dict[str, object]] = []
    active_key: tuple[str, str, str, int] | None = None
    later_events: list[object] = []

    for event in ordered.itertuples(index=False):
        key = (event.sample_id, event.chrom, event.allele, event.contig)
        if key != active_key:
            active_key = key
            later_events = []

        # Events are visited newest-to-oldest. Reversing this list tests the
        # nearest (earliest) subsequent event first.
        hiding_event = None
        ducking_mode = None
        for later_event in reversed(later_events):
            crosses_wgd = any(
                event.depth < wgd_depth < later_event.depth
                for wgd_depth in wgd_depths_by_sample.get(event.sample_id, ())
            )
            if crosses_wgd and not include_wgd_crossing:
                continue

            later_is_deletion = "Deletion" in later_event.event_type
            event_is_deletion = "Deletion" in event.event_type
            later_is_duplication = "Duplication" in later_event.event_type
            event_is_duplication = "Duplication" in event.event_type
            boundaries_touch = (
                event.end == later_event.start
                or later_event.end == event.start
            )

            if (
                later_is_deletion
                and later_event.start <= event.start
                and event.end <= later_event.end
            ):
                hiding_event = later_event
                ducking_mode = "contained_in_deletion"
                break

            is_neighboring_loss = (
                event_is_deletion and later_is_deletion and boundaries_touch
            )
            if is_neighboring_loss:
                hiding_event = later_event
                ducking_mode = "neighboring_loss"
                break

            is_neighboring_gain = (
                event_is_duplication and later_is_duplication and boundaries_touch
            )
            if is_neighboring_gain:
                hiding_event = later_event
                ducking_mode = "neighboring_gain"
                break

        if hiding_event is not None:
            rescuing_event = None
            if ducking_mode == "contained_in_deletion":
                for candidate in reversed(later_events):
                    if candidate.depth >= hiding_event.depth:
                        continue
                    candidate_is_duplication = (
                        "Duplication" in candidate.event_type
                    )
                    overlaps_hidden_event = any(
                        hidden_start < candidate_end
                        and candidate_start < hidden_end
                        for hidden_start, hidden_end in event.intervals
                        for candidate_start, candidate_end in candidate.intervals
                    )
                    if candidate_is_duplication and overlaps_hidden_event:
                        rescuing_event = candidate
                        break

            results.append(
                {
                    "sample_id": event.sample_id,
                    "chrom": event.chrom,
                    "allele": event.allele,
                    "contig": event.contig,
                    "ducking_mode": ducking_mode,
                    "hidden_depth": event.depth,
                    "hidden_event_type": event.event_type,
                    "hidden_start": event.start,
                    "hidden_end": event.end,
                    "hidden_width": event.width,
                    "hidden_region_count": event.region_count,
                    "hidden_unique_region_count": event.unique_region_count,
                    "hidden_description": event.description,
                    "is_rescued": rescuing_event is not None,
                    "rescuing_depth": (
                        rescuing_event.depth
                        if rescuing_event is not None
                        else None
                    ),
                    "rescuing_event_type": (
                        rescuing_event.event_type
                        if rescuing_event is not None
                        else None
                    ),
                    "rescuing_start": (
                        rescuing_event.start
                        if rescuing_event is not None
                        else None
                    ),
                    "rescuing_end": (
                        rescuing_event.end
                        if rescuing_event is not None
                        else None
                    ),
                    "rescuing_width": (
                        rescuing_event.width
                        if rescuing_event is not None
                        else None
                    ),
                    "rescuing_region_count": (
                        rescuing_event.region_count
                        if rescuing_event is not None
                        else None
                    ),
                    "rescuing_unique_region_count": (
                        rescuing_event.unique_region_count
                        if rescuing_event is not None
                        else None
                    ),
                    "rescuing_description": (
                        rescuing_event.description
                        if rescuing_event is not None
                        else None
                    ),
                    "hiding_depth": hiding_event.depth,
                    "hiding_event_type": hiding_event.event_type,
                    "hiding_start": hiding_event.start,
                    "hiding_end": hiding_event.end,
                    "hiding_width": hiding_event.width,
                    "hiding_region_count": hiding_event.region_count,
                    "hiding_unique_region_count": hiding_event.unique_region_count,
                    "hiding_description": hiding_event.description,
                }
            )

        later_events.append(event)

    columns = [
        "sample_id",
        "chrom",
        "allele",
        "contig",
        "ducking_mode",
        "hidden_depth",
        "hidden_event_type",
        "hidden_start",
        "hidden_end",
        "hidden_width",
        "hidden_region_count",
        "hidden_unique_region_count",
        "hidden_description",
        "is_rescued",
        "rescuing_depth",
        "rescuing_event_type",
        "rescuing_start",
        "rescuing_end",
        "rescuing_width",
        "rescuing_region_count",
        "rescuing_unique_region_count",
        "rescuing_description",
        "hiding_depth",
        "hiding_event_type",
        "hiding_start",
        "hiding_end",
        "hiding_width",
        "hiding_region_count",
        "hiding_unique_region_count",
        "hiding_description",
    ]
    ducking_events = pd.DataFrame.from_records(results, columns=columns)
    if ducking_events.empty:
        return ducking_events
    nullable_rescuing_columns = [
        "rescuing_depth",
        "rescuing_start",
        "rescuing_end",
        "rescuing_width",
        "rescuing_region_count",
        "rescuing_unique_region_count",
    ]
    ducking_events[nullable_rescuing_columns] = ducking_events[
        nullable_rescuing_columns
    ].astype("Int64")
    return ducking_events.sort_values(
        ["sample_id", "chrom", "allele", "contig", "hidden_depth"],
        kind="stable",
        ignore_index=True,
    )


def rescued_event_mask(ducking_events: pd.DataFrame) -> pd.Series:
    """Return a Boolean mask selecting rescued rows in a detailed table."""

    if "is_rescued" not in ducking_events.columns:
        return pd.Series(False, index=ducking_events.index, dtype=bool)

    values = ducking_events["is_rescued"]
    if pd.api.types.is_bool_dtype(values.dtype):
        rescued = values.fillna(False)
    else:
        normalized = values.astype("string").str.strip().str.lower()
        valid_values = {"true", "false", "1", "0", "yes", "no"}
        invalid = normalized.notna() & ~normalized.isin(valid_values)
        if invalid.any():
            bad_value = values.loc[invalid].iloc[0]
            raise ValueError(f"Invalid is_rescued value: {bad_value!r}")
        rescued = normalized.isin({"true", "1", "yes"})
    return rescued.astype(bool)


def unrescued_ducking_events(ducking_events: pd.DataFrame) -> pd.DataFrame:
    """Return ducking-event rows that have no overlapping rescue duplication.

    Detailed tables created before rescue detection did not contain an
    ``is_rescued`` column. They remain compatible and are treated as entirely
    unrescued.
    """

    return ducking_events.loc[~rescued_event_mask(ducking_events)].copy()

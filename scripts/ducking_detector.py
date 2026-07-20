"""Command-line tool for detecting and filtering ducking events."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

if __package__:
    from .lib_ducking import (
        find_ducking_events,
        unrescued_ducking_events,
    )
else:
    from lib_ducking import (
        find_ducking_events,
        unrescued_ducking_events,
    )


DEFAULT_DATA_DIR = Path.cwd()

_EVENT_ID_COLUMNS = ["sample_id", "depth", "event_type", "description"]
_DUCKING_ID_COLUMNS = [
    "sample_id",
    "hidden_depth",
    "hidden_event_type",
    "hidden_description",
]


def filter_ducking_events(
    events: pd.DataFrame,
    ducking_events: pd.DataFrame,
) -> pd.DataFrame:
    """Return the event table without unrescued ducking events."""

    missing_event_columns = set(_EVENT_ID_COLUMNS).difference(events.columns)
    if missing_event_columns:
        missing = ", ".join(sorted(missing_event_columns))
        raise ValueError(f"Missing required event columns: {missing}")

    missing_ducking_columns = set(_DUCKING_ID_COLUMNS).difference(
        ducking_events.columns
    )
    if missing_ducking_columns:
        missing = ", ".join(sorted(missing_ducking_columns))
        raise ValueError(f"Missing required ducking-event columns: {missing}")

    hidden_event_ids = (
        unrescued_ducking_events(ducking_events)
        .loc[:, _DUCKING_ID_COLUMNS]
        .rename(
            columns={
                "hidden_depth": "depth",
                "hidden_event_type": "event_type",
                "hidden_description": "description",
            }
        )
        .drop_duplicates()
        .assign(_is_ducking_event=True)
    )
    marked_events = events.merge(
        hidden_event_ids,
        on=_EVENT_ID_COLUMNS,
        how="left",
        sort=False,
    )
    return marked_events.loc[
        marked_events["_is_ducking_event"].isna(), events.columns
    ].reset_index(drop=True)


def _add_input_output_arguments(
    parser: argparse.ArgumentParser,
    *,
    default_output: Path,
    output_help: str,
) -> None:
    parser.add_argument(
        "-I",
        "-i",
        "--input",
        type=Path,
        default=DEFAULT_DATA_DIR / "events.tsv",
        help="Input event TSV (default: ./events.tsv)",
    )
    parser.add_argument(
        "-O",
        "-o",
        "--output",
        type=Path,
        default=default_output,
        help=output_help,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""

    parser = argparse.ArgumentParser(
        prog="python scripts/ducking_detector.py",
        description="Detect ducking events or remove them from an event table.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    detect_parser = subparsers.add_parser(
        "detect",
        help="detect ducking events",
        description="Detect events hidden by later same-contig events.",
    )
    _add_input_output_arguments(
        detect_parser,
        default_output=DEFAULT_DATA_DIR / "ducking_events.tsv",
        output_help="Detailed ducking-event TSV (default: ./ducking_events.tsv)",
    )
    filter_parser = subparsers.add_parser(
        "filter",
        help="remove unrescued ducking events",
        description="Remove unrescued ducking events from an event TSV.",
    )
    _add_input_output_arguments(
        filter_parser,
        default_output=DEFAULT_DATA_DIR / "events_filtered.tsv",
        output_help=(
            "Filtered event TSV (default: ./events_filtered.tsv)"
        ),
    )
    for command_parser in (detect_parser, filter_parser):
        command_parser.add_argument(
            "--include-wgd-crossing",
            action="store_true",
            help=(
                "allow ducking matches across an intervening whole-genome "
                "doubling"
            ),
        )
    return parser


def _write_tsv(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)


def main(argv: Sequence[str] | None = None) -> None:
    """Run the selected command."""

    args = build_parser().parse_args(argv)
    events = pd.read_csv(args.input, sep="\t")
    ducking_events = find_ducking_events(
        events,
        include_wgd_crossing=args.include_wgd_crossing,
    )

    if args.command == "detect":
        _write_tsv(ducking_events, args.output)
        print(
            f"Detected {len(ducking_events):,} ducking events "
            f"({int(ducking_events['is_rescued'].sum()):,} rescued); "
            f"wrote {args.output}"
        )
        return

    filtered_events = filter_ducking_events(events, ducking_events)
    _write_tsv(filtered_events, args.output)
    removed_count = len(events) - len(filtered_events)
    rescued_count = len(ducking_events) - len(
        unrescued_ducking_events(ducking_events)
    )
    print(
        f"Removed {removed_count:,} ducking events; "
        f"preserved {rescued_count:,} rescued events; "
        f"wrote {len(filtered_events):,} events to {args.output}"
    )


if __name__ == "__main__":
    main()

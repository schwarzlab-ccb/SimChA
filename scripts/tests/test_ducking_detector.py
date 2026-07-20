import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.ducking_detector import filter_ducking_events, main


def make_events(rows):
    return pd.DataFrame(
        rows,
        columns=[
            "sample_id",
            "event_type",
            "depth",
            "description",
            "regions_gained",
            "regions_lost",
        ],
    )


class FilterDuckingEventsTests(unittest.TestCase):
    def test_removes_exact_hidden_events_and_preserves_input_schema(self):
        events = make_events(
            [
                (
                    "S1",
                    "InternalDuplication",
                    1,
                    "contig:7;gain",
                    "[H1:chr1[20:30)]",
                    "[]",
                ),
                (
                    "S1",
                    "InternalDeletion",
                    2,
                    "contig:7;loss",
                    "[]",
                    "[H1:chr1[10:40)]",
                ),
                (
                    "S2",
                    "InternalDuplication",
                    1,
                    "contig:7;gain",
                    "[H1:chr1[20:30)]",
                    "[]",
                ),
            ]
        )
        ducking_events = pd.DataFrame(
            {
                "sample_id": ["S1", "S1"],
                "hidden_depth": [1, 1],
                "hidden_event_type": [
                    "InternalDuplication",
                    "InternalDuplication",
                ],
                "hidden_description": ["contig:7;gain", "contig:7;gain"],
            }
        )

        result = filter_ducking_events(events, ducking_events)

        self.assertEqual(result.columns.tolist(), events.columns.tolist())
        self.assertEqual(
            result.loc[:, ["sample_id", "depth"]].values.tolist(),
            [["S1", 2], ["S2", 1]],
        )

    def test_preserves_rescued_ducking_events(self):
        events = make_events(
            [
                (
                    "S1",
                    "InternalDuplication",
                    1,
                    "contig:7;rescued",
                    "[H1:chr1[20:30)]",
                    "[]",
                ),
                (
                    "S1",
                    "InternalDuplication",
                    2,
                    "contig:7;unrescued",
                    "[H1:chr1[40:50)]",
                    "[]",
                ),
            ]
        )
        ducking_events = pd.DataFrame(
            {
                "sample_id": ["S1", "S1"],
                "hidden_depth": [1, 2],
                "hidden_event_type": [
                    "InternalDuplication",
                    "InternalDuplication",
                ],
                "hidden_description": [
                    "contig:7;rescued",
                    "contig:7;unrescued",
                ],
                "is_rescued": [True, False],
            }
        )

        result = filter_ducking_events(events, ducking_events)

        self.assertEqual(result["depth"].tolist(), [1])


class CommandLineTests(unittest.TestCase):
    def test_detect_and_filter_commands_support_short_input_output_options(self):
        events = make_events(
            [
                (
                    "S1",
                    "InternalDuplication",
                    1,
                    "contig:7;gain",
                    "[H1:chr1[20:30)]",
                    "[]",
                ),
                (
                    "S1",
                    "InternalDeletion",
                    2,
                    "contig:7;loss",
                    "[]",
                    "[H1:chr1[10:40)]",
                ),
            ]
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            input_path = directory / "events.tsv"
            ducking_path = directory / "ducking.tsv"
            filtered_path = directory / "filtered.tsv"
            events.to_csv(input_path, sep="\t", index=False)

            main(["detect", "-I", str(input_path), "-O", str(ducking_path)])
            main(
                [
                    "filter",
                    "-I",
                    str(input_path),
                    "-O",
                    str(filtered_path),
                ]
            )

            detected = pd.read_csv(ducking_path, sep="\t")
            filtered = pd.read_csv(filtered_path, sep="\t")

        self.assertEqual(len(detected), 1)
        self.assertEqual(len(filtered), 1)
        self.assertEqual(filtered.loc[0, "depth"], 2)

if __name__ == "__main__":
    unittest.main()

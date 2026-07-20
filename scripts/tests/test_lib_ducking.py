"""Tests for the ducking-event detection core."""

import unittest

import pandas as pd

from scripts.lib_ducking import (
    find_ducking_events,
    parse_contig_id,
    parse_regions,
    unrescued_ducking_events,
)


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


class ParseRegionsTests(unittest.TestCase):
    def test_parses_half_open_regions(self):
        regions = parse_regions("[H1:chr1[10:20),H2:chrX[30:50)]")

        self.assertEqual(
            [(region.allele, region.chrom, region.start, region.end) for region in regions],
            [("H1", "chr1", 10, 20), ("H2", "chrX", 30, 50)],
        )

    def test_rejects_non_positive_width(self):
        with self.assertRaises(ValueError):
            parse_regions("[H1:chr1[10:10)]")

    def test_parses_contig_identifier(self):
        self.assertEqual(parse_contig_id("contig:27;length:100;"), 27)


class FindDuckingEventsTests(unittest.TestCase):
    def test_finds_event_inside_later_same_allele_deletion(self):
        events = make_events(
            [
                ("S1", "InternalDuplication", 1, "contig:7;gain", "[H1:chr1[20:30)]", "[]"),
                ("S1", "InternalDeletion", 2, "contig:7;loss", "[]", "[H1:chr1[10:40)]"),
            ]
        )

        result = find_ducking_events(events)

        self.assertEqual(len(result), 1)
        self.assertEqual(result.loc[0, "ducking_mode"], "contained_in_deletion")
        self.assertEqual(result.loc[0, "hidden_depth"], 1)
        self.assertEqual(result.loc[0, "hiding_depth"], 2)

    def test_uses_outer_deletion_span_across_previously_lost_gap(self):
        events = make_events(
            [
                ("S1", "InternalDeletion", 1, "contig:7;first", "[]", "[H1:chr1[20:30)]"),
                (
                    "S1",
                    "InternalDeletion",
                    2,
                    "contig:7;second",
                    "[]",
                    "[H1:chr1[10:20),H1:chr1[30:40)]",
                ),
            ]
        )

        result = find_ducking_events(events)

        self.assertEqual(len(result), 1)
        self.assertEqual(result.loc[0, "hidden_description"], "contig:7;first")
        self.assertEqual(result.loc[0, "hiding_start"], 10)
        self.assertEqual(result.loc[0, "hiding_end"], 40)

    def test_does_not_match_other_allele_or_partial_overlap(self):
        events = make_events(
            [
                ("S1", "InternalDuplication", 1, "contig:7;gain", "[H1:chr1[20:30)]", "[]"),
                ("S1", "InternalDeletion", 2, "contig:7;other allele", "[]", "[H2:chr1[10:40)]"),
                ("S1", "InternalDeletion", 3, "contig:7;partial", "[]", "[H1:chr1[25:40)]"),
            ]
        )

        result = find_ducking_events(events)

        self.assertTrue(result.empty)

    def test_chooses_earliest_subsequent_containing_deletion(self):
        events = make_events(
            [
                ("S1", "InternalDuplication", 1, "contig:7;gain", "[H1:chr1[20:30)]", "[]"),
                ("S1", "InternalDeletion", 2, "contig:7;near", "[]", "[H1:chr1[10:40)]"),
                ("S1", "ChromDeletion", 3, "contig:7;far", "[]", "[H1:chr1[0:100)]"),
            ]
        )

        result = find_ducking_events(events)

        hidden_gain = result.loc[result["hidden_depth"] == 1].iloc[0]
        self.assertEqual(hidden_gain["hiding_depth"], 2)

    def test_ignores_wgd_as_a_local_hidden_event(self):
        events = make_events(
            [
                ("S1", "WholeGenomeDoubling", 1, "wgd", "[H1:chr1[0:100)]", "[]"),
                ("S1", "ChromDeletion", 2, "contig:7;loss", "[]", "[H1:chr1[0:100)]"),
            ]
        )

        result = find_ducking_events(events)

        self.assertTrue(result.empty)

    def test_wgd_blocks_loss_by_one_later_deletion_by_default(self):
        events = make_events(
            [
                ("S1", "InternalDuplication", 1, "contig:7;gain", "[H1:chr1[20:30)]", "[]"),
                ("S1", "WholeGenomeDoubling", 2, "wgd", "[H1:chr1[0:100)]", "[]"),
                ("S1", "ChromDeletion", 3, "contig:7;loss", "[]", "[H1:chr1[0:100)]"),
            ]
        )

        result = find_ducking_events(events)
        coordinate_only_result = find_ducking_events(
            events,
            include_wgd_crossing=True,
        )

        self.assertTrue(result.empty)
        self.assertEqual(len(coordinate_only_result), 1)

    def test_does_not_match_a_deletion_on_a_different_contig(self):
        events = make_events(
            [
                ("S1", "InternalDuplication", 1, "contig:7;gain", "[H1:chr1[20:30)]", "[]"),
                ("S1", "InternalDeletion", 2, "contig:8;loss", "[]", "[H1:chr1[10:40)]"),
            ]
        )

        result = find_ducking_events(events)

        self.assertTrue(result.empty)

    def test_finds_two_neighboring_loss_events(self):
        events = make_events(
            [
                ("S1", "InternalDeletion", 1, "contig:7;left", "[]", "[H1:chr1[10:20)]"),
                ("S1", "InternalDeletion", 2, "contig:7;right", "[]", "[H1:chr1[20:30)]"),
            ]
        )

        result = find_ducking_events(events)

        self.assertEqual(len(result), 1)
        self.assertEqual(result.loc[0, "ducking_mode"], "neighboring_loss")
        self.assertEqual(result.loc[0, "hidden_depth"], 1)
        self.assertEqual(result.loc[0, "hiding_depth"], 2)

    def test_neighboring_mode_requires_two_losses_with_no_gap(self):
        events = make_events(
            [
                ("S1", "InternalDuplication", 1, "contig:7;gain", "[H1:chr1[10:20)]", "[]"),
                ("S1", "InternalDeletion", 2, "contig:7;loss", "[]", "[H1:chr1[20:30)]"),
                ("S2", "InternalDeletion", 1, "contig:7;left", "[]", "[H1:chr1[10:20)]"),
                ("S2", "InternalDeletion", 2, "contig:7;right", "[]", "[H1:chr1[21:30)]"),
            ]
        )

        result = find_ducking_events(events)

        self.assertTrue(result.empty)

    def test_finds_two_neighboring_gain_events(self):
        events = make_events(
            [
                (
                    "S1",
                    "InternalDuplication",
                    1,
                    "contig:7;left",
                    "[H1:chr1[10:20)]",
                    "[]",
                ),
                (
                    "S1",
                    "InternalDuplication",
                    2,
                    "contig:7;right",
                    "[H1:chr1[20:30)]",
                    "[]",
                ),
            ]
        )

        result = find_ducking_events(events)

        self.assertEqual(len(result), 1)
        self.assertEqual(result.loc[0, "ducking_mode"], "neighboring_gain")
        self.assertEqual(result.loc[0, "hidden_depth"], 1)
        self.assertEqual(result.loc[0, "hiding_depth"], 2)

    def test_finds_earliest_overlapping_duplication_before_hiding_deletion(self):
        events = make_events(
            [
                (
                    "S1",
                    "InternalDeletion",
                    1,
                    "contig:7;hidden",
                    "[]",
                    "[H1:chr1[20:30)]",
                ),
                (
                    "S1",
                    "InternalDuplication",
                    2,
                    "contig:7;first rescue",
                    "[H1:chr1[25:35)]",
                    "[]",
                ),
                (
                    "S1",
                    "InternalDuplication",
                    3,
                    "contig:7;second rescue",
                    "[H1:chr1[15:25)]",
                    "[]",
                ),
                (
                    "S1",
                    "ChromDeletion",
                    4,
                    "contig:7;hider",
                    "[]",
                    "[H1:chr1[0:100)]",
                ),
            ]
        )

        result = find_ducking_events(events)
        hidden_event = result.loc[result["hidden_depth"] == 1].iloc[0]

        self.assertTrue(hidden_event["is_rescued"])
        self.assertEqual(hidden_event["rescuing_depth"], 2)
        self.assertEqual(
            hidden_event["rescuing_description"],
            "contig:7;first rescue",
        )
        self.assertEqual(hidden_event["hiding_depth"], 4)

    def test_does_not_rescue_touching_or_post_deletion_duplications(self):
        events = make_events(
            [
                (
                    "S1",
                    "InternalDeletion",
                    1,
                    "contig:7;hidden",
                    "[]",
                    "[H1:chr1[20:30)]",
                ),
                (
                    "S1",
                    "InternalDuplication",
                    2,
                    "contig:7;touching",
                    "[H1:chr1[30:40)]",
                    "[]",
                ),
                (
                    "S1",
                    "ChromDeletion",
                    3,
                    "contig:7;hider",
                    "[]",
                    "[H1:chr1[0:100)]",
                ),
                (
                    "S1",
                    "InternalDuplication",
                    4,
                    "contig:7;too late",
                    "[H1:chr1[20:30)]",
                    "[]",
                ),
            ]
        )

        result = find_ducking_events(events)
        hidden_event = result.loc[result["hidden_depth"] == 1].iloc[0]

        self.assertFalse(hidden_event["is_rescued"])
        self.assertTrue(pd.isna(hidden_event["rescuing_depth"]))

    def test_does_not_rescue_duplication_only_in_gap_between_fragments(self):
        events = make_events(
            [
                (
                    "S1",
                    "InternalDeletion",
                    1,
                    "contig:7;fragmented hidden event",
                    "[]",
                    "[H1:chr1[10:20),H1:chr1[30:40)]",
                ),
                (
                    "S1",
                    "InternalDuplication",
                    2,
                    "contig:7;gap only",
                    "[H1:chr1[20:30)]",
                    "[]",
                ),
                (
                    "S1",
                    "ChromDeletion",
                    3,
                    "contig:7;hider",
                    "[]",
                    "[H1:chr1[0:100)]",
                ),
            ]
        )

        result = find_ducking_events(events)
        hidden_event = result.loc[result["hidden_depth"] == 1].iloc[0]

        self.assertFalse(hidden_event["is_rescued"])
        self.assertTrue(pd.isna(hidden_event["rescuing_depth"]))


class UnrescuedDuckingEventsTests(unittest.TestCase):
    def test_excludes_rescued_rows_and_supports_legacy_tables(self):
        ducking_events = pd.DataFrame(
            {"sample_id": ["S1", "S2"], "is_rescued": [True, False]}
        )
        legacy_ducking_events = ducking_events.drop(columns="is_rescued")

        result = unrescued_ducking_events(ducking_events)

        self.assertEqual(result["sample_id"].tolist(), ["S2"])
        self.assertEqual(len(unrescued_ducking_events(legacy_ducking_events)), 2)


if __name__ == "__main__":
    unittest.main()

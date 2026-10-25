"""出院致残队列补统计：筛选、患者去重、原排序、有效分母及分享边界。"""
from __future__ import annotations

import contextlib
import io
from pathlib import Path
import unittest
from unittest.mock import patch
import zipfile

from mrs_summary import MONTHS, aggregate_groups
from roi_summary import analyse_masks, load_atlas
from run_summary import csv_rows, digest, public, write_json
from supplement_bgic import load_original_states
from supplement_discharge_cohort import execute, group_curve, review_pages, summarize
import test_supplement_deep_union as source_tests


class DischargeSupplementTest(unittest.TestCase):
    write_image = source_tests.DeepUnionSupplementTest.write_image

    def setUp(self):
        source_tests.DeepUnionSupplementTest.setUp(self)
        self.records = {identifier: {"D_DIAG": diagnosis, "D_MRS": mrs}
                        for identifier, diagnosis, mrs in zip(self.ids[:6], ["1", "1", "1", "1", "2", "1"],
                                                             ["1", "2", "5", "6", "2", "3"])}
        self.source_audit = {"rows": 6, "id_storage": "numeric", "original_sha256_matched": True,
                             "schema": {key: {"actual": key, "type": "double", "label": key} for key in ("D_DIAG", "D_MRS")}}
        self.loader = patch("supplement_discharge_cohort.load_discharge_records", return_value=(self.records, self.source_audit))
        self.loader.start(); self.addCleanup(self.loader.stop)

    def test_cohort_and_union_deduplicated_and_original_order_preserved(self):
        result = summarize(self.run)
        self.assertEqual(result["cohort_n"], 4)
        self.assertEqual(result["union_roi"]["any_overlap_n"], 3)
        self.assertEqual(result["union_roi"]["n"], 3)
        self.assertIsNone(result["union_roi"]["label_id"])
        self.assertEqual(result["overlap_audit"]["membership_multiplicity_histogram"], {"0": 1, "1": 0, "2": 3, "3": 0, "4": 0})
        expected = [row["roi_id"] for row in sorted(self.report["roi"]["roi_rows"], key=lambda row: (-row["any_overlap_n"], row["roi_name"]))]
        self.assertEqual([row["roi_id"] for row in result["ranked_rois_13"]], expected)
        self.assertEqual([row["rank"] for row in result["ranked_rois_13"]], list(range(1, 14)))
        self.assertEqual(len(result["ranked_outcome_rates_13"]), 78)
        self.assertEqual(result["validation"]["component_and_overall_mrs_cells_checked"], 1680)
        for row in result["ranked_rois_13"]:
            self.assertEqual(row["denominator_n"], 4)
            self.assertAlmostEqual(row["any_overlap_proportion"], row["any_overlap_n"]/4)
            self.assertNotIn("total_overlap_ml", row)
        curve = {row["month"]: row for row in result["overall_curve"]}
        self.assertEqual((curve[3]["valid_n"], curve[3]["good_recovery_n"], curve[3]["functional_disability_n"], curve[3]["functional_loss_n"]), (4, 1, 1, 2))
        self.assertEqual(curve[12]["valid_n"], 3)
        self.assertEqual(curve[12]["death6_n"], 0)
        self.assertEqual(curve[12]["missing_n"], 1)
        self.assertEqual(len(review_pages(result)), 10)
        for row in result["overall_curve"]+result["union_curve"]:
            self.assertEqual(row["good_recovery_n"]+row["functional_disability_n"]+row["functional_loss_n"], row["valid_n"])
            self.assertAlmostEqual(sum(row[key] for key in ("good_recovery_proportion", "functional_disability_proportion", "functional_loss_proportion")), 1.)

    def test_discharge_death_conflicts_are_new_subset_states_not_source_edits(self):
        self.records[self.ids[2]]["D_MRS"] = "6"
        before = digest(self.run/"private/derived_states.csv")
        result = summarize(self.run)
        self.assertEqual(result["status"], "COMPLETED_WITH_FLAGS")
        self.assertEqual(result["discharge_death_conflict_audit"]["new_conflict_patient_n"], 1)
        self.assertEqual(result["discharge_death_conflict_audit"]["new_conflict_visit_n"], 6)
        self.assertEqual(result["overall_curve"][0]["valid_n"], 3)
        self.assertEqual(result["overall_curve"][0]["conflict_n"], 1)
        self.assertEqual(before, digest(self.run/"private/derived_states.csv"))

    def test_union_sensitivity_uses_union_voxels_not_thresholded_member_union(self):
        # 7个BG体素＋6个CSO体素，共0.104mL；单区均未达到0.1mL。
        import numpy as np
        identifier = self.ids[4]
        values = np.zeros(self.shape, dtype=np.uint8)
        values.ravel()[:7] = 1; values.ravel()[94:100] = 1
        self.write_image(Path(self.mapping[identifier]), values)
        self.records[identifier] = {"D_DIAG": "1", "D_MRS": "2"}
        atlas = load_atlas(self.report["config"]["atlas_path"], self.report["config"]["atlas_labels"], "independent_binary_channels")
        original_roi = analyse_masks(self.mapping, atlas, .1)
        states = load_original_states(self.run/"private/derived_states.csv")
        self.report.update(roi=public(original_roi),
            primary=aggregate_groups(states, original_roi["private_any_memberships"]),
            sensitivity=aggregate_groups(states, original_roi["private_memberships"]))
        write_json(self.report_path, self.report)
        csv_rows(self.run/"private/mask_inventory.csv", [
            {"clinical_id": key, "path": path, "status": "matched", "bytes": Path(path).stat().st_size}
            for key, path in self.mapping.items()])
        result = summarize(self.run)
        self.assertEqual(result["union_roi"]["n"], 4)
        rows = {row["roi_id"]: row for row in result["ranked_rois_13"]}
        self.assertEqual(rows["roi_0001"]["n"], 1)
        self.assertEqual(rows["roi_0012"]["n"], 1)
        self.assertEqual(result["union_roi"]["any_overlap_n"], 4)

    def test_empty_selection_stops_before_scanning_and_releases_no_stats(self):
        for row in self.records.values():
            row["D_MRS"] = "0"
        with patch("supplement_discharge_cohort.analyse_masks") as analyse:
            with contextlib.redirect_stdout(io.StringIO()):
                result, output = execute(self.run, self.root/"empty")
            analyse.assert_not_called()
        self.assertEqual(result["status"], "FAILED")
        self.assertNotIn("cohort_n", result)
        self.assertFalse((output/"share.zip").exists())

    def test_any_original_roi_drift_blocks(self):
        row = next(row for row in self.report["primary"]["distributions"] if row["group"] == "roi_0002" and row["month"] == 3 and row["state"] == "0")
        row["n"] += 1
        write_json(self.report_path, self.report)
        with self.assertRaisesRegex(ValueError, "six-visit mRS distributions"):
            summarize(self.run)

    def test_original_diagnosis_drift_blocks_before_scan(self):
        self.report["diagnosis_counts"] = {"infarction_1": 6, "TIA_2": 0, "other_or_missing": 0}
        write_json(self.report_path, self.report)
        with patch("supplement_discharge_cohort.analyse_masks") as analyse:
            with self.assertRaisesRegex(ValueError, "diagnosis counts"):
                summarize(self.run)
            analyse.assert_not_called()

    def test_source_reader_failure_blocks_without_publishing(self):
        with patch("supplement_discharge_cohort.load_discharge_records", side_effect=ValueError("early SAS fingerprint changed")):
            with contextlib.redirect_stdout(io.StringIO()):
                result, output = execute(self.run, self.root/"source_drift")
        self.assertEqual(result["status"], "FAILED")
        self.assertFalse((output/"share.zip").exists())
        self.assertTrue((output/"private/failure.txt").is_file())

    def test_zero_denominator_is_none_not_zero_percent(self):
        states = {"test": {month: "missing" for month in MONTHS}}
        rows = group_curve(aggregate_groups(states, {}), "overall")
        for row in rows:
            self.assertEqual(row["valid_n"], 0)
            self.assertIsNone(row["good_recovery_proportion"])
            self.assertIsNone(row["disability_death_proportion"])

    def test_share_privacy_original_files_unchanged_and_review_pages(self):
        before = {path: digest(path) for path in self.root.rglob("*") if path.is_file()}
        with contextlib.redirect_stdout(io.StringIO()):
            result, output = execute(self.run, self.root/"complete")
        self.assertEqual(result["status"], "COMPLETED")
        self.assertTrue((output/"private/selection.csv").is_file())
        self.assertEqual(len(load_original_states(output/"private/derived_states.csv")), 4)
        with zipfile.ZipFile(output/"share.zip") as archive:
            self.assertIsNone(archive.testzip())
            self.assertIn("discharge_cohort.json", archive.namelist())
            self.assertIn("page_10.txt", archive.namelist())
            for name in archive.namelist():
                value = archive.read(name).decode("utf-8-sig")
                for identifier in self.ids:
                    self.assertNotIn(identifier, value)
                self.assertNotIn("patient_id", value)
                self.assertNotIn(str(self.root), value)
        for path, checksum in before.items():
            self.assertEqual(digest(path), checksum)
        with self.assertRaises(FileExistsError):
            execute(self.run, output)


if __name__ == "__main__":
    unittest.main()

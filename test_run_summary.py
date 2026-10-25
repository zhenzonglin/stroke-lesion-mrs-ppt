"""整条汇总流程的合成测试；真实小影像、模拟SAS读取，绝不使用患者数据。"""
from __future__ import annotations

import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import nibabel as nib
import numpy as np

from run_summary import execute, review_pages, outcome_display_settings
from configure_workstation import create_config
from test_inputs import reader


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class EndToEndSummaryTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.atlas_dir = self.root / "atlas"
        self.atlas_dir.mkdir()
        self.mask_dir = self.root / "masks"
        self.mask_dir.mkdir()
        self.affine = np.diag([2., 2., 2., 1.])
        self.atlas = np.zeros((5, 5, 5), dtype=np.int16)
        self.atlas[:2] = 1
        self.atlas[3:] = 2
        self.atlas_path = self.atlas_dir / "test_atlas.nii.gz"
        self.write_image(self.atlas_path, self.atlas)
        self.labels = self.atlas_dir / "labels.tsv"
        self.labels.write_text("label_id\troi_name\n0\tBackground\n1\tTest region A\n2\tTest region B\n3\tAbsent region\n", encoding="utf-8")
        self.provenance = self.atlas_dir / "provenance.json"
        self.provenance.write_text(json.dumps({"atlas_sha256": sha(self.atlas_path),
            "labels_sha256": sha(self.labels), "origin": "synthetic software test only"}), encoding="utf-8")
        self.identifiers = (101234567, 202345678, 303456789)
        for identifier, overlaps in zip(self.identifiers, ({1: 13}, {1: 7, 2: 15}, {0: 5})):
            values = np.zeros(self.atlas.shape, dtype=np.uint8)
            for label, count in overlaps.items():
                coordinates = np.argwhere(self.atlas == label)[:count]
                values[tuple(coordinates.T)] = 1
            self.write_image(self.mask_dir / f"{identifier}_lesion.nii.gz", values)
        for name in ("early.sas7bdat", "late.sas7bdat"):
            (self.root / name).write_bytes(b"Synthetic placeholder; reader is mocked in this test")
        ids = [float(identifier) for identifier in self.identifiers]
        self.specs = {
            "early.sas7bdat": {"data": {"code_n": ids, "F3_MRS": [3., 1., 0.],
                "F12_MRS": [1., 6., 0.], "D_DIAG": [1., 1., 2.],
                "D_DEATH": [1., 1., 1.], "F3_DEATH": [1., 1., 1.], "F12_DEATH": [1., 2., 1.]}},
            "late.sas7bdat": {"data": {"code_n": ids, "m24_mrs": [0., None, 0.],
                "m36_mrs": [1., 6., 0.], "m48_mrs": [2., None, 0.], "m60_mrs": [6., 6., 0.],
                "F2Y_DEATH": [1., 2., 1.]}},
        }
        self.cfg = {"source_files": {"early": str(self.root / "early.sas7bdat"),
                                    "late": str(self.root / "late.sas7bdat")},
            "lesion_root": str(self.mask_dir), "id_mode": "numeric_code_n",
            "id_regex": r"^(?:sub-)?(?P<id>[0-9]+)(?:$|[_-])",
            "atlas_path": "atlas/test_atlas.nii.gz", "atlas_labels": "atlas/labels.tsv",
            "atlas_provenance": "atlas/provenance.json", "allow_unknown_units_if_exact_grid": False,
            "primary_membership": "any_positive_overlap", "sensitivity_min_overlap_ml": 0.1,
            "death_definition": "mrs_equal_6_only_at_all_visits", "carry_forward": False,
            "diagnosis_filter": None, "output_root": "outputs"}
        self.config = self.root / "config.json"

    def write_image(self, path, values, affine=None):
        image = nib.Nifti1Image(values, self.affine if affine is None else affine)
        image.header.set_xyzt_units("mm")
        nib.save(image, path)

    def run_execute(self):
        self.config.write_text(json.dumps(self.cfg), encoding="utf-8")
        self.originals = {p: sha(p) for p in self.root.rglob("*") if p.is_file()}
        columns = []
        with patch.dict(sys.modules, {"pyreadstat": reader(self.specs, columns)}), contextlib.redirect_stdout(io.StringIO()):
            report, output = execute(self.config, self.root / "result")
        return report, output, columns

    def assert_source_preservation(self):
        for path, digest in self.originals.items():
            self.assertEqual(sha(path), digest, str(path))

    def assert_share_has_only_aggregate(self, output):
        with zipfile.ZipFile(output / "share.zip") as archive:
            self.assertIsNone(archive.testzip())
            self.assertIn("ppt_data.json", archive.namelist())
            self.assertIn("manifest.json", archive.namelist())
            for name in archive.namelist():
                self.assertNotIn("private/", name)
                self.assertNotIn("mask_inventory", name)
                self.assertNotIn("derived_states", name)
                text = archive.read(name).decode("utf-8-sig")
                for identifier in self.identifiers:
                    self.assertNotIn(str(identifier), text, name)
                self.assertNotIn('"private_', text, name)
        self.assertTrue((output / "private/mask_inventory.csv").is_file())

    def test_completed_all_rois_mrs_only_and_input_hashes(self):
        report, output, columns = self.run_execute()
        self.assertEqual(report["status"], "COMPLETED", report["issues"])
        self.assertFalse(any("DEATH" in field.upper() for field in columns))
        self.assertEqual(report["roi"]["valid_n"], 3)
        self.assertEqual(report["roi"]["coverage_summary"]["no_atlas_overlap_n"], 1)
        roi_rows = {row["label_id"]: row for row in report["roi"]["roi_rows"]}
        self.assertEqual([roi_rows[label]["any_overlap_n"] for label in (1, 2, 3)], [2, 1, 0])
        self.assertEqual([roi_rows[label]["n"] for label in (1, 2, 3)], [1, 1, 0])
        self.assertTrue(all(row["denominator_n"] == 3 for row in roi_rows.values()))
        self.assertEqual(report["diagnosis_counts"]["TIA_2"], 1)
        self.assertEqual(report["mrs_audit"]["derived_counts"]["3"]["6"], 0)
        self.assertEqual(report["mrs_audit"]["derived_counts"]["12"]["6"], 1)
        self.assertEqual(report["mrs_audit"]["derived_counts"]["24"]["missing"], 1)
        self.assertEqual(report["mrs_audit"]["derived_counts"]["60"]["6"], 2)
        for version in ("primary", "sensitivity"):
            self.assertTrue(report[version]["validation"]["all_passed"])
            for group in report[version]["groups"]:
                for month in (3, 12, 24, 36, 48, 60):
                    rows = [r for r in report[version]["distributions"] if r["group"] == group["group"] and r["month"] == month]
                    self.assertEqual(sum(row["n"] for row in rows), group["n"])
                for month in (3, 12, 24, 36, 48):
                    rows = [r for r in report[version]["transitions"] if r["group"] == group["group"] and r["from_month"] == month]
                    self.assertEqual(sum(row["n"] for row in rows), group["n"])
        self.assert_share_has_only_aggregate(output)
        self.assert_source_preservation()

    def test_missing_mask_path_returns_input_report_without_primary(self):
        self.cfg["lesion_root"] = str(self.root / "missing_mask_folder")
        report, output, _ = self.run_execute()
        self.assertEqual(report["status"], "INPUTS_REQUIRED")
        self.assertNotIn("primary", report)
        self.assertNotIn("roi", report)
        self.assertTrue((output / "share/ppt_data.json").is_file())
        self.assertFalse((output / "share/primary_distributions.csv").exists())
        self.assert_share_has_only_aggregate(output)
        self.assert_source_preservation()

    def test_all_masks_wrong_space_never_claims_completed(self):
        shifted = self.affine.copy()
        shifted[0, 3] = 0.5
        for path in self.mask_dir.iterdir():
            self.write_image(path, np.ones(self.atlas.shape, dtype=np.uint8), shifted)
        report, output, _ = self.run_execute()
        self.assertEqual(report["status"], "NO_VALID_MASKS", report["issues"])
        self.assertEqual(report["roi"]["exclusion_counts"], {"grid_mismatch": 3})
        self.assertEqual(report["roi"]["valid_n"], 0)
        self.assertNotIn("primary", report)
        self.assertNotIn("mrs_audit", report)
        self.assertFalse((output / "share/primary_distributions.csv").exists())
        self.assertTrue(all(row["proportion"] is None for row in report["roi"]["roi_rows"]))
        self.assert_share_has_only_aggregate(output)
        self.assert_source_preservation()

    def test_v2_overlapping_channels_and_all_visit_outcome_csv(self):
        masks = np.stack([self.atlas == 1, (self.atlas == 1) | (self.atlas == 2)], axis=3).astype(np.uint8)
        self.write_image(self.atlas_path, masks)
        self.labels.write_text("label_id\troi_name\n1\tBasal Ganglia\n2\tInternal Capsule\n", encoding="utf-8")
        self.provenance.write_text(json.dumps({"atlas_sha256": sha(self.atlas_path),
            "labels_sha256": sha(self.labels), "version": "synthetic_v2"}), encoding="utf-8")
        self.cfg.update(atlas_format="independent_binary_channels",
                        display_outcome={"definition": "mrs_2_6", "months": [3, 12, 24, 36, 48, 60]})
        report, output, _ = self.run_execute()
        self.assertEqual(report["status"], "COMPLETED", report["issues"])
        self.assertEqual(report["atlas"]["atlas_format"], "independent_binary_channels")
        self.assertEqual(report["roi"]["coverage_summary"]["covered_voxels_total"], 35)
        self.assertEqual([r["any_overlap_n"] for r in report["roi"]["roi_rows"]], [2, 2])
        self.assertEqual(len(report["roi_outcome_rates"]), 2*6*3)
        selected = [r for r in report["roi_outcome_rates"] if r["group"] == "roi_0001"
                    and r["month"] == 24 and r["definition"] == "mrs_2_6"][0]
        self.assertEqual((selected["numerator_n"], selected["denominator_n"], selected["total_n"]), (0, 1, 2))
        self.assertEqual(selected["lesion_n"], 2)
        self.assertEqual(selected["lesion_denominator_n"], 3)
        self.assertEqual(selected["death6_n"], 0)  # 既往6分后缺失不得前填死亡。
        self.assertTrue((output/"share/roi_outcome_rates.csv").is_file())
        self.assertTrue((output/"share/primary_recovery_distributions.csv").is_file())
        pages = review_pages(report)
        self.assertEqual(sum("个月脑图标注" in title for title, _ in pages), 6)
        self.assert_share_has_only_aggregate(output)
        self.assert_source_preservation()

    def test_configuration_v1_preserved_v2_explicit_and_exclusive(self):
        v1_path = self.root/"config.v1.local.json"
        v1 = create_config(self.root, v1_path)
        self.assertEqual(v1["atlas_format"], "discrete_labels")
        self.assertIsNone(v1["display_outcome"])
        self.assertIn("display_MNI9_HObrainstem_2mm", v1["atlas_path"])
        v2 = create_config(self.root, self.root/"config.v2.local.json", atlas_version="v2")
        self.assertEqual(v2["atlas_format"], "independent_binary_channels")
        self.assertTrue(Path(v2["atlas_path"]).is_absolute())
        self.assertEqual(outcome_display_settings(v2), ([3, 12, 24, 36, 48, 60], "mrs_2_6"))
        before = v1_path.read_bytes()
        with self.assertRaises(FileExistsError):
            create_config(self.root, v1_path, atlas_version="v2")
        self.assertEqual(before, v1_path.read_bytes())
        with self.assertRaises(ValueError):
            create_config(self.root, self.root/"bad.json", display_month=60)
        with self.assertRaises(ValueError):
            outcome_display_settings({"display_outcome": {"months": [3, 3], "definition": "mrs_2_6"}})


if __name__ == "__main__":
    unittest.main()

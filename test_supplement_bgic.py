"""并集补统计的合成测试；不读取真实患者数据。"""
from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

import nibabel as nib
import numpy as np

from mrs_summary import MONTHS, aggregate_groups
from roi_summary import analyse_masks, load_atlas
from run_summary import csv_rows, digest, public, write_json
from supplement_bgic import MERGED_ID, execute, summarize


class BasalGangliaInternalCapsuleSupplementTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root/"original"
        (self.source/"share").mkdir(parents=True)
        (self.source/"private").mkdir()
        self.affine = np.diag([2., 2., 2., 1.])
        self.shape = (5, 5, 5)
        bg = np.zeros(self.shape, dtype=bool); bg[:2] = True
        ic = np.zeros(self.shape, dtype=bool); ic[1:3] = True
        other = np.zeros(self.shape, dtype=bool); other[3] = True
        self.atlas_path = self.root/"atlas.nii.gz"
        self.write_image(self.atlas_path, np.stack((bg, ic, other), axis=3).astype(np.uint8))
        self.labels = self.root/"labels.tsv"
        self.labels.write_text("label_id\troi_name\n1\tBasal Ganglia\n2\tInternal Capsule\n3\tOther region\n", encoding="utf8")
        self.atlas = load_atlas(self.atlas_path, self.labels, "independent_binary_channels")
        self.ids = ["700111001", "700222002", "700333003", "700444004", "700555005"]
        self.mapping = {}
        arrangements = [{0:13}, {2:13}, {0:10, 2:10}, {4:2}, {}]
        for identifier, arrangement in zip(self.ids, arrangements):
            values = np.zeros(self.shape, dtype=np.uint8)
            for axis, count in arrangement.items():
                values[axis].ravel()[:count] = 1
            path = self.root/(identifier+"_lesion.nii.gz")
            self.write_image(path, values)
            self.mapping[identifier] = str(path)
        self.states = {identifier: dict(zip(MONTHS, values)) for identifier, values in zip(self.ids[:4], [
            ["2", "6", "missing", "6", "6", "6"],
            ["0", "1", "2", "2", "3", "6"],
            ["4", "3", "2", "missing", "5", "5"],
            ["1", "1", "1", "1", "1", "1"],
        ])}
        roi = analyse_masks(self.mapping, self.atlas, .1)
        self.report = {"status": "COMPLETED_WITH_FLAGS", "run_id": "synthetic_original",
            "config": {"atlas_path": str(self.atlas_path), "atlas_labels": str(self.labels),
                "death_definition": "mrs_equal_6_only_at_all_visits", "carry_forward": False,
                "primary_membership": "any_positive_overlap", "sensitivity_min_overlap_ml": .1,
                "allow_unknown_units_if_exact_grid": False},
            "atlas": public(self.atlas), "roi": public(roi),
            "primary": aggregate_groups(self.states, roi["private_any_memberships"]),
            "sensitivity": aggregate_groups(self.states, roi["private_memberships"])}
        self.report_path = self.source/"share/ppt_data.json"
        write_json(self.report_path, self.report)
        csv_rows(self.source/"private/derived_states.csv", [{"patient_id": i, **{f"m{m}":v for m,v in visits.items()}}
                                                        for i,visits in self.states.items()])
        csv_rows(self.source/"private/mask_inventory.csv", [{"id":i, "clinical_id":i, "path":p,
            "status":"matched", "bytes":Path(p).stat().st_size} for i,p in self.mapping.items()])
        csv_rows(self.source/"private/mask_exclusions.csv", [{"patient_id":i,"path":self.mapping[i],"reason":reason}
                    for i,reason in roi["private_exclusions"].items()])

    def write_image(self, path: Path, values: np.ndarray) -> None:
        image = nib.Nifti1Image(values, self.affine)
        image.header.set_xyzt_units("mm")
        nib.save(image, path)

    def test_patient_union_not_sum_and_union_volume_threshold(self):
        result = summarize(self.source)
        self.assertEqual(result["status"], "COMPLETED")
        self.assertEqual(result["cohort_n"], 4)
        self.assertEqual(result["inclusion_exclusion"], {
            "basal_ganglia_n":2, "internal_capsule_n":2, "intersection_n":1, "union_n":3,
            "formula":"union_n = basal_ganglia_n + internal_capsule_n - intersection_n"})
        self.assertEqual(result["merged_roi"]["any_overlap_n"], 3)
        self.assertIsNone(result["merged_roi"]["label_id"])
        self.assertEqual(result["merged_roi"]["component_label_ids"], [1,2])
        self.assertEqual(result["merged_roi"]["n"], 3)  # 双侧各0.08mL，并集0.16mL，不能用阈值患者集合OR。
        self.assertEqual(result["merged_roi"]["any_overlap_proportion"], .75)
        self.assertEqual([r["roi_id"] for r in result["ranked_rois"]], [MERGED_ID,"roi_0003"])
        self.assertEqual([r["rank"] for r in result["ranked_rois"]], [1,2])
        self.assertEqual(len(result["ranked_outcome_rates"]), 12)
        selected = {r["month"]:r for r in result["merged_outcome_rates"] if r["definition"]=="mrs_2_6"}
        self.assertEqual((selected[3]["numerator_n"], selected[3]["denominator_n"]), (2,3))
        self.assertEqual((selected[24]["numerator_n"], selected[24]["denominator_n"], selected[24]["death6_n"]), (2,2,0))
        self.assertEqual((selected[60]["numerator_n"], selected[60]["denominator_n"], selected[60]["death6_n"]), (3,3,2))
        self.assertTrue(result["validation"]["valid_patient_set_reproduced"])
        self.assertFalse(result["validation"]["original_mask_sha256_available"])

    def test_empty_mask_remains_excluded_and_changed_valid_set_rejected(self):
        self.assertEqual(self.report["roi"]["exclusion_counts"], {"empty_mask":1})
        self.assertEqual(summarize(self.source)["cohort_n"],4)
        # 同一文件保持可读但空掩膜变非空，可能首先由size审计或valid IDs检查阻断。
        self.write_image(Path(self.mapping[self.ids[-1]]), np.ones(self.shape,dtype=np.uint8))
        with self.assertRaisesRegex(ValueError, "mask size changed|valid patient set"):
            summarize(self.source)

    def test_component_report_drift_is_rejected(self):
        self.report["roi"]["roi_rows"][0]["any_overlap_n"] += 1
        write_json(self.report_path,self.report)
        with self.assertRaisesRegex(ValueError,"Component any-overlap count"):
            summarize(self.source)

    def test_component_mrs_drift_is_rejected(self):
        row=next(r for r in self.report["primary"]["distributions"] if r["group"]=="roi_0001" and r["month"]==3 and r["state"]=="0")
        row["n"] += 1
        write_json(self.report_path,self.report)
        with self.assertRaisesRegex(ValueError,"six-visit mRS distributions"):
            summarize(self.source)

    def test_share_aggregates_only_source_preserved_and_no_overwrite(self):
        original={p:digest(p) for p in self.root.rglob('*') if p.is_file()}
        with contextlib.redirect_stdout(io.StringIO()):
            result,out=execute(self.source,self.root/'supplement')
        self.assertEqual(result["status"],"COMPLETED")
        with zipfile.ZipFile(out/'share.zip') as archive:
            self.assertIsNone(archive.testzip())
            self.assertIn('bgic_supplement.json',archive.namelist())
            self.assertIn('ranked_outcome_rates.csv',archive.namelist())
            for name in archive.namelist():
                value=archive.read(name).decode('utf-8-sig')
                for identifier in self.ids:self.assertNotIn(identifier,value,name)
                self.assertNotIn('patient_id',value,name)
                self.assertNotIn(str(self.root),value,name)
        for path,old in original.items():self.assertEqual(digest(path),old,str(path))
        with self.assertRaises(FileExistsError):execute(self.source,out)

    def test_failure_does_not_release_partial_merged_counts(self):
        self.report['roi']['roi_rows'][0]['any_overlap_n'] += 1
        write_json(self.report_path,self.report)
        with contextlib.redirect_stdout(io.StringIO()):
            result,out=execute(self.source,self.root/'failed_supplement')
        self.assertEqual(result['status'],'FAILED')
        self.assertNotIn('merged_roi',result)
        self.assertFalse((out/'share.zip').exists())
        self.assertTrue((out/'private/failure.txt').exists())


if __name__ == '__main__':
    unittest.main()

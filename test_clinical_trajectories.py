"""仅SAS队列测试，所有影像路径故意不存在。"""
from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import zipfile

from run_clinical_trajectories import ClinicalInputsError, execute, review_pages, summarize, digest
import test_inputs as fixtures


class ClinicalTrajectoriesTest(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        self.root=Path(temporary.name)
        self.ids=[810000001+i for i in range(8)]
        self.paths={role:self.root/(role+".sas7bdat") for role in ("early","late")}
        for path in self.paths.values():path.write_bytes(b"synthetic_source")
        self.specs={
            "early.sas7bdat":{"data":{
                "code_n":self.ids[:7],"D_DIAG":[1,1,1,1,2,1,1],"D_MRS":[0,2,99,None,4,6,5],
                "F3_MRS":[0,2,4,1,5,1,6],"F12_MRS":[1,None,3,2,5,None,6],"D_DEATH":[2]*7}},
            "late.sas7bdat":{"data":{
                "code_n":self.ids[:6]+[self.ids[7]],"m24_mrs":[2,3,1,None,0,6,4],
                "m36_mrs":[1,2,1,3,0,None,3],"m48_mrs":[0,3,2,3,0,0,2],
                "m60_mrs":[0,4,3,4,0,6,1],"F5Y_DEATH":[2]*7}},
        }
        self.config={"source_files":{key:str(value) for key,value in self.paths.items()},
            "lesion_root":str(self.root/"forbidden_masks"),"atlas_path":str(self.root/"forbidden_atlas.nii.gz"),
            "atlas_labels":str(self.root/"forbidden_labels.tsv"),"diagnosis_filter":"ignored_old_config",
            "output_root":"old_imaging_outputs"}
        self.config_path=self.root/"config.v2.local.json"
        self.config_path.write_text(json.dumps(self.config),encoding="utf8")
        self.observed=[]

    def mocked(self):
        return patch.dict(sys.modules,{"pyreadstat":fixtures.reader(self.specs,self.observed)})

    def run_summary(self):
        with self.mocked():return summarize(self.config,self.root)

    def test_no_imaging_io_and_cohorts_include_early_only(self):
        original_stat=Path.stat
        def guarded_stat(path,*args,**kwargs):
            if "forbidden_" in str(path):raise AssertionError("Imaging path was accessed")
            return original_stat(path,*args,**kwargs)
        with patch.object(Path,"stat",guarded_stat),patch("inputs.scan_masks",side_effect=AssertionError("scan called")),patch("roi_summary.analyse_masks",side_effect=AssertionError("image analysis called")):
            result=self.run_summary()
        self.assertEqual((result["overall_n"],result["disabled_n"]),(6,3))
        flow=result["selection"]
        self.assertEqual(flow["clinical_union_n"],8)
        self.assertEqual(flow["diagnosis_counts"]["no_early_record"],1)
        self.assertEqual(flow["diagnosis_counts"]["TIA_2"],1)
        self.assertEqual(flow["disabled_exclusions"],{"discharge_mrs_0_1":1,"discharge_mrs_missing":1,"discharge_mrs_invalid":1})
        self.assertEqual(result["clinical_audit"]["early_only_n"],1)
        self.assertTrue(result["validation"]["no_imaging_filter"])
        self.assertIn("D_MRS",self.observed)
        self.assertFalse(any("DEATH" in name for name in self.observed))

    def test_group_counts_denominators_death_conflicts_no_imputation(self):
        result=self.run_summary()
        overall={r["month"]:r for r in result["overall_curve"]}
        disabled={r["month"]:r for r in result["disabled_curve"]}
        self.assertEqual((overall[3]["valid_n"],overall[3]["good_n"],overall[3]["disability_n"],overall[3]["loss_n"]),(5,2,1,2))
        self.assertEqual((disabled[3]["valid_n"],disabled[3]["good_n"],disabled[3]["disability_n"],disabled[3]["loss_n"]),(2,0,1,1))
        self.assertEqual(overall[3]["conflict_n"],1)
        self.assertEqual(disabled[12]["missing_n"],2)
        self.assertEqual(disabled[12]["death6_n"],1)
        self.assertEqual(disabled[24]["missing_n"],1)
        self.assertEqual(disabled[48]["conflict_n"],1)
        # 48月已由24月的6分标冲突，出院6分只新标记3月那次。
        self.assertEqual(result["discharge_death_conflict_audit"]["new_conflict_visit_n"],1)
        for key in ("overall_curve","disabled_curve"):
            for row in result[key]:
                self.assertEqual(row["good_n"]+row["disability_n"]+row["loss_n"],row["valid_n"])
                self.assertAlmostEqual(row["good_percent"]+row["disability_percent"]+row["loss_percent"],100)
                self.assertEqual(row["valid_n"]+row["missing_n"]+row["invalid_n"]+row["conflict_n"],row["total_n"])

    def test_missing_required_discharge_field_blocks(self):
        del self.specs["early.sas7bdat"]["data"]["D_MRS"]
        with self.mocked(),contextlib.redirect_stdout(io.StringIO()):
            result,out=execute(self.config_path,self.root/"failed")
        self.assertEqual(result["status"],"INPUTS_REQUIRED")
        self.assertNotIn("overall_curve",result)
        self.assertFalse((out/"share.zip").exists())
        self.assertIn("D_MRS",str(result["clinical_audit"]))

    def test_duplicate_ids_block_without_partial_data(self):
        self.specs["early.sas7bdat"]["data"]["code_n"][1]=self.ids[0]
        with self.assertRaises(ClinicalInputsError):self.run_summary()

    def test_mixed_id_storage_blocks(self):
        self.specs["late.sas7bdat"]["numeric"]=False
        self.specs["late.sas7bdat"]["data"]["code_n"]=[str(i) for i in self.specs["late.sas7bdat"]["data"]["code_n"]]
        with self.assertRaisesRegex(ValueError,"storage types differ"):self.run_summary()

    def test_source_mutation_during_read_blocks(self):
        source=fixtures.reader(self.specs,self.observed)
        def read(path,**kwargs):
            result=source.read_sas7bdat(path,**kwargs)
            if Path(path).name=="late.sas7bdat" and kwargs.get("usecols"):
                self.paths["late"].write_bytes(b"modified_after_read")
            return result
        with patch.dict(sys.modules,{"pyreadstat":SimpleNamespace(read_sas7bdat=read)}):
            with self.assertRaisesRegex(ValueError,"SAS changed"):summarize(self.config,self.root)

    def test_zero_subgroup_or_visit_is_na_not_zero_percent(self):
        self.specs["early.sas7bdat"]["data"]["D_MRS"]=[0]*7
        result=self.run_summary()
        self.assertEqual(result["disabled_n"],0)
        self.assertEqual(result["status"],"COMPLETED_WITH_FLAGS")
        for row in result["disabled_curve"]:
            self.assertIsNone(row["good_percent"])
            self.assertEqual(row["valid_n"],0)
        self.specs["early.sas7bdat"]["data"]["F12_MRS"]=[None]*7
        result=self.run_summary()
        row=next(r for r in result["overall_curve"] if r["month"]==12)
        self.assertIsNone(row["good_percent"])
        self.assertEqual(row["valid_n"],0)

    def test_no_ischemic_patients_blocks(self):
        self.specs["early.sas7bdat"]["data"]["D_DIAG"]=[2]*7
        with self.assertRaisesRegex(ValueError,"No patients meet"):self.run_summary()

    def test_share_privacy_inputs_preserved_and_new_output_only(self):
        before={p:digest(p) for p in self.root.rglob("*") if p.is_file()}
        with self.mocked(),contextlib.redirect_stdout(io.StringIO()):
            result,out=execute(self.config_path,self.root/"new_clinical")
        self.assertEqual(result["overall_n"],6)
        self.assertTrue((out/"private/selection.csv").is_file())
        with zipfile.ZipFile(out/"share.zip") as archive:
            self.assertIsNone(archive.testzip())
            self.assertIn("clinical_trajectories.json",archive.namelist())
            self.assertIn("overall_curve.csv",archive.namelist())
            self.assertIn("disabled_curve.csv",archive.namelist())
            for name in archive.namelist():
                value=archive.read(name).decode("utf-8-sig")
                self.assertNotIn(str(self.root),value)
                self.assertNotIn("patient_id",value)
                self.assertNotIn("private_",value)
                for identifier in self.ids:self.assertNotIn(str(identifier),value)
        for path,value in before.items():self.assertEqual(digest(path),value)
        self.assertTrue(all(len(lines)<=24 for _,lines in review_pages(result)))
        with self.assertRaises(FileExistsError):execute(self.config_path,out)


if __name__=="__main__":unittest.main()

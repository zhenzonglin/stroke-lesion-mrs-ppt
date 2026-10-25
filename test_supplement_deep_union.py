"""四组成区联合人群补统计：保留13区，验证去重、分母、复现与分享范围。"""
from __future__ import annotations

import contextlib
import io
from pathlib import Path
import tempfile
import unittest
import zipfile

import nibabel as nib
import numpy as np

from mrs_summary import MONTHS, aggregate_groups
from roi_summary import analyse_masks, load_atlas
from run_summary import csv_rows, digest, public, write_json
from supplement_deep_union import UNION_ID, execute, summarize


class DeepUnionSupplementTest(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup)
        self.root=Path(temp.name);self.run=self.root/"original"
        (self.run/"share").mkdir(parents=True);(self.run/"private").mkdir()
        self.affine=np.diag([2.,2.,2.,1.]);self.shape=(6,6,6)
        masks=np.zeros((*self.shape,13),dtype=np.uint8)
        for channel,start in ((0,0),(6,20),(10,40),(11,60)):
            values=np.zeros(self.shape,dtype=np.uint8)
            values.ravel()[start:start+40]=1;values.ravel()[190]=1
            masks[...,channel]=values
        for channel in set(range(13))-{0,6,10,11}:
            values=np.zeros(self.shape,dtype=np.uint8);values.ravel()[120:125]=1;masks[...,channel]=values
        atlas_path=self.root/"atlas.nii.gz";self.write_image(atlas_path,masks)
        names={1:"Basal Ganglia",7:"Internal Capsule",11:"Corona Radiata",12:"Centrum Semiovale"}
        labels=self.root/"labels.tsv"
        labels.write_text("label_id\troi_name\n"+"".join(f"{i}\t{names.get(i,'Other region '+str(i))}\n" for i in range(1,14)),encoding="utf8")
        atlas=load_atlas(atlas_path,labels,"independent_binary_channels")
        self.ids=[f"81000000{i}" for i in range(1,8)]
        self.mapping={}
        for identifier,indices in zip(self.ids,[range(0,13),range(25,38),range(45,58),range(65,78),[190],range(120,125),[]]):
            values=np.zeros(self.shape,dtype=np.uint8);values.ravel()[list(indices)]=1
            path=self.root/(identifier+"_lesion.nii.gz");self.write_image(path,values);self.mapping[identifier]=str(path)
        states={identifier:dict(zip(MONTHS,values)) for identifier,values in zip(self.ids[:6],[
            ["0","1","2","6","missing","6"],
            ["2","3","2","2","3","4"],
            ["4","4","3","2","1","0"],
            ["6","missing","6","6","6","6"],
            ["1","0","missing","3","4","5"],
            ["0","0","0","0","0","0"],
        ])}
        roi=analyse_masks(self.mapping,atlas,.1)
        self.report={"status":"COMPLETED_WITH_FLAGS","run_id":"synthetic_four_union",
            "config":{"atlas_path":str(atlas_path),"atlas_labels":str(labels),
                "death_definition":"mrs_equal_6_only_at_all_visits","carry_forward":False,
                "primary_membership":"any_positive_overlap","sensitivity_min_overlap_ml":.1,
                "allow_unknown_units_if_exact_grid":False},
            "atlas":public(atlas),"roi":public(roi),
            "primary":aggregate_groups(states,roi["private_any_memberships"]),
            "sensitivity":aggregate_groups(states,roi["private_memberships"])}
        self.report_path=self.run/"share/ppt_data.json";write_json(self.report_path,self.report)
        csv_rows(self.run/"private/derived_states.csv",[{"patient_id":i,**{f"m{m}":v for m,v in visits.items()}} for i,visits in states.items()])
        csv_rows(self.run/"private/mask_inventory.csv",[{"clinical_id":i,"path":p,"status":"matched","bytes":Path(p).stat().st_size} for i,p in self.mapping.items()])
        csv_rows(self.run/"private/mask_exclusions.csv",[{"patient_id":i,"reason":r} for i,r in roi["private_exclusions"].items()])

    def write_image(self,path:Path,values:np.ndarray):
        im=nib.Nifti1Image(values,self.affine);im.header.set_xyzt_units("mm");nib.save(im,path)

    def test_four_component_union_and_original_13_rois_preserved(self):
        result=summarize(self.run)
        self.assertEqual(result["cohort_n"],6)
        self.assertEqual(result["union_roi"]["any_overlap_n"],5)
        self.assertEqual(result["union_roi"]["n"],4)
        self.assertIsNone(result["union_roi"]["label_id"])
        self.assertEqual(set(result["union_roi"]["component_label_ids"]),{1,7,11,12})
        self.assertEqual(result["overlap_audit"]["membership_multiplicity_histogram"],{"0":1,"1":1,"2":3,"3":0,"4":1})
        self.assertEqual(result["overlap_audit"]["sum_component_memberships_n"],11)
        self.assertEqual(result["overlap_audit"]["overcount_if_added_n"],6)
        self.assertEqual(result["overlap_audit"]["four_way_inclusion_exclusion_n"],5)
        original={r["roi_id"] for r in self.report["roi"]["roi_rows"]}
        self.assertEqual({r["roi_id"] for r in result["ranked_rois_13"]},original)
        self.assertNotIn(UNION_ID,original)
        self.assertEqual(len(result["ranked_outcome_rates_13"]),78)
        self.assertEqual([r["rank"] for r in result["ranked_rois_13"]],list(range(1,14)))
        self.assertEqual(result["validation"]["component_and_overall_mrs_cells_checked"],600)

    def test_curve_partitions_valid_denominator_and_no_death_fill(self):
        result=summarize(self.run);curve={r["month"]:r for r in result["union_curve"]}
        self.assertEqual((curve[3]["valid_n"],curve[3]["good_recovery_n"],curve[3]["functional_disability_n"],curve[3]["functional_loss_n"]),(5,2,1,2))
        self.assertEqual(curve[12]["valid_n"],4)
        self.assertEqual(curve[12]["death6_n"],0)
        self.assertEqual(curve[48]["valid_n"],4)
        self.assertEqual(curve[48]["death6_n"],1)
        for row in curve.values():
            self.assertEqual(sum(row[k] for k in ("good_recovery_n","functional_disability_n","functional_loss_n")),row["valid_n"])
            self.assertAlmostEqual(sum(row[k] for k in ("good_recovery_proportion","functional_disability_proportion","functional_loss_proportion")),1.)
        self.assertEqual(len(result["union_distributions"]),60)
        self.assertEqual(len(result["union_recovery_distributions"]),18)

    def test_component_drift_fails_without_partial_union_release(self):
        self.report["roi"]["roi_rows"][10]["any_overlap_n"]+=1
        write_json(self.report_path,self.report)
        with contextlib.redirect_stdout(io.StringIO()):result,out=execute(self.run,self.root/"failed")
        self.assertEqual(result["status"],"FAILED")
        self.assertNotIn("union_roi",result)
        self.assertFalse((out/"share.zip").exists())

    def test_fourth_component_sensitivity_mrs_drift_fails(self):
        row=next(r for r in self.report["sensitivity"]["distributions"] if r["group"]=="roi_0012" and r["month"]==60 and r["state"]=="0")
        row["n"]+=1;write_json(self.report_path,self.report)
        with self.assertRaisesRegex(ValueError,"six-visit mRS distributions"):summarize(self.run)

    def test_share_has_no_patient_ids_and_sources_are_unchanged(self):
        before={p:digest(p) for p in self.root.rglob('*') if p.is_file()}
        with contextlib.redirect_stdout(io.StringIO()):result,out=execute(self.run,self.root/"complete")
        self.assertEqual(result["status"],"COMPLETED")
        with zipfile.ZipFile(out/"share.zip") as archive:
            self.assertIsNone(archive.testzip())
            self.assertIn("deep_union_supplement.json",archive.namelist())
            self.assertIn("union_curve.csv",archive.namelist())
            for name in archive.namelist():
                value=archive.read(name).decode('utf-8-sig')
                for identifier in self.ids:self.assertNotIn(identifier,value)
                self.assertNotIn('patient_id',value)
                self.assertNotIn(str(self.root),value)
        for path,checksum in before.items():self.assertEqual(digest(path),checksum)
        with self.assertRaises(FileExistsError):execute(self.run,out)


if __name__ == '__main__':unittest.main()

"""出院致残性缺血性卒中筛选：使用合成字段，不读取患者数据。"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from discharge_selection import (apply_discharge_death_conflicts, load_discharge_records,
                                 select_discharge_cohort)
from mrs_summary import MONTHS
from run_summary import digest
from test_inputs import Column, Frame


class DischargeSelectionTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.source = self.root / "early.sas7bdat"
        self.source.write_bytes(b"synthetic early SAS")
        self.data = {"code_n": [810001., 810002., 810003.], "D_DIAG": [1., 2., 1.],
                     "D_MRS": [2., 6., None], "D_DEATH": [2, 1, 2]}
        self.numeric = True
        self.observed = []
        self.report = {
            "config": {"source_files": {"early": str(self.source)}, "id_width": None},
            "source_fingerprints": [{"path": str(self.source), "bytes": self.source.stat().st_size,
                                     "sha256": digest(self.source)}],
            "clinical": {"sources": {"early": {"rows": 3, "id_storage": "numeric"}}},
        }

    def reader(self, path, metadataonly=False, usecols=None, **kwargs):
        self.assertEqual(Path(path), self.source)
        types = {field: "string" if field.casefold() == "code_n" and not self.numeric else "double"
                 for field in self.data}
        meta = SimpleNamespace(column_names=list(self.data), readstat_variable_types=types,
                               column_names_to_labels={"D_MRS": "Discharge modified Rankin Scale"})
        if metadataonly:
            return None, meta
        self.observed.extend(usecols)
        self.assertTrue(kwargs["user_missing"])
        self.assertTrue(kwargs["disable_datetime_conversion"])
        return Frame({field: Column(self.data[field], "float64" if self.numeric else "object")
                      for field in usecols}), meta

    def load(self, reader=None):
        with patch.dict(sys.modules, {"pyreadstat": SimpleNamespace(read_sas7bdat=reader or self.reader)}):
            return load_discharge_records(self.report)

    def test_frozen_source_whitelist_and_anonymous_audit(self):
        records, audit = self.load()
        self.assertEqual(records["810001"], {"D_DIAG": "1", "D_MRS": "2"})
        self.assertEqual(records["810003"]["D_MRS"], "missing")
        self.assertEqual(set(self.observed), {"code_n", "D_DIAG", "D_MRS"})
        self.assertEqual(audit["schema"]["D_MRS"]["label"], "Discharge modified Rankin Scale")
        self.assertTrue(audit["original_sha256_matched"])
        self.assertTrue(audit["unchanged_during_read"])
        self.assertEqual(audit["id_unique_n"], 3)
        serialized = json.dumps(audit)
        for token in (str(self.root), "810001", "810002", "810003", "D_DEATH"):
            self.assertNotIn(token, serialized)

    def test_source_changed_before_read_blocks(self):
        self.source.write_bytes(b"different source")
        with self.assertRaisesRegex(ValueError, "original run fingerprint"):
            self.load()
        self.assertEqual(self.observed, [])

    def test_source_changed_during_read_blocks(self):
        def changed(*args, **kwargs):
            result = self.reader(*args, **kwargs)
            if not kwargs.get("metadataonly"):
                self.source.write_bytes(b"changed during data read")
            return result
        with self.assertRaisesRegex(ValueError, "changed while"):
            self.load(changed)

    def test_missing_or_ambiguous_fingerprint_blocks(self):
        self.report["source_fingerprints"] *= 2
        with self.assertRaisesRegex(ValueError, "missing or ambiguous"):
            self.load()
        self.report["source_fingerprints"] = []
        with self.assertRaisesRegex(ValueError, "missing or ambiguous"):
            self.load()

    def test_case_insensitive_fields_and_numeric_width(self):
        self.data["CODE_N"] = self.data.pop("code_n")
        self.data["d_mrs"] = self.data.pop("D_MRS")
        self.report["config"]["id_width"] = 8
        records, audit = self.load()
        self.assertIn("00810001", records)
        self.assertEqual(audit["schema"]["D_MRS"]["actual"], "d_mrs")

    def test_missing_or_case_duplicate_discharge_field_blocks(self):
        missing = self.data.pop("D_MRS")
        with self.assertRaisesRegex(ValueError, "exist exactly once"):
            self.load()
        self.data.update(D_MRS=missing, d_mrs=missing)
        with self.assertRaisesRegex(ValueError, "exist exactly once"):
            self.load()

    def test_bad_ids_and_storage_drift_block(self):
        for values in ([810001, 810001, 810003], [810001, None, 810003],
                       [810001, 2.5, 810003]):
            with self.subTest(values=values):
                self.data["code_n"] = values
                with self.assertRaisesRegex(ValueError, "patient IDs"):
                    self.load()
        self.numeric = False
        self.data["code_n"] = ["00810001", "00810002", "00810003"]
        with self.assertRaisesRegex(ValueError, "ID storage differs"):
            self.load()
        self.report["clinical"]["sources"]["early"]["id_storage"] = "character"
        records, _ = self.load()
        self.assertIn("00810001", records)

    def test_invalid_discharge_values_are_not_threshold_selected(self):
        self.data["D_MRS"] = [99, 2.5, True]
        records, _ = self.load()
        self.assertEqual({row["D_MRS"] for row in records.values()}, {"invalid"})
        selected, audit = select_discharge_cohort({key: {} for key in records}, records)
        self.assertFalse(selected)
        self.assertEqual(audit["mutually_exclusive_exclusions"]["ischemic_discharge_mrs_invalid"], 2)

    def test_selection_flow_keeps_discharge6_and_ignores_followup_missingness(self):
        records = {f"include{n}": {"D_DIAG": "1", "D_MRS": str(n)} for n in range(2, 7)}
        records.update({
            "tia": {"D_DIAG": "2", "D_MRS": "6"},
            "diag_missing": {"D_DIAG": "missing", "D_MRS": "5"},
            "diag_invalid": {"D_DIAG": "invalid", "D_MRS": "3"},
            "diag_other": {"D_DIAG": "0", "D_MRS": "4"},
            "discharge_missing": {"D_DIAG": "1", "D_MRS": "missing"},
            "discharge_invalid": {"D_DIAG": "1", "D_MRS": "invalid"},
            "discharge0": {"D_DIAG": "1", "D_MRS": "0"},
            "discharge1": {"D_DIAG": "1", "D_MRS": "1"},
            "outside_imaging": {"D_DIAG": "1", "D_MRS": "4"},
        })
        states = {key: {m: "missing" for m in MONTHS} for key in records if key != "outside_imaging"}
        states["no_early_record"] = {m: "0" for m in MONTHS}
        selected, audit = select_discharge_cohort(states, records)
        self.assertEqual(selected, {f"include{n}" for n in range(2, 7)})
        self.assertEqual(audit["original_valid_n"], 14)
        self.assertEqual(audit["included_n"], 5)
        self.assertEqual(audit["ischemic_n"], 9)
        self.assertEqual(audit["included_discharge_mrs_counts"], {str(n): 1 for n in range(2, 7)})
        self.assertEqual(audit["discharge_mrs6_included_n"], 1)
        self.assertEqual(sum(audit["mutually_exclusive_exclusions"].values()), 9)
        self.assertTrue(audit["flow_conservation_passed"])
        self.assertFalse(audit["followup_completeness_required"])

    def test_discharge6_conflicts_are_visit_local_and_copy_only(self):
        states = {"dead": dict(zip(MONTHS, ["0", "missing", "invalid", "6", "conflict", "5"])),
                  "disabled": dict(zip(MONTHS, ["0", "1", "2", "3", "4", "5"]))}
        records = {"dead": {"D_DIAG": "1", "D_MRS": "6"},
                   "disabled": {"D_DIAG": "1", "D_MRS": "2"}}
        before = copy.deepcopy(states)
        updated, audit = apply_discharge_death_conflicts(states, records)
        self.assertEqual(states, before)
        self.assertEqual(list(updated["dead"].values()), ["conflict", "missing", "invalid", "6", "conflict", "conflict"])
        self.assertEqual(updated["disabled"], before["disabled"])
        self.assertEqual(audit["new_conflict_patient_n"], 1)
        self.assertEqual(audit["new_conflict_visit_n"], 2)
        self.assertEqual(audit["new_conflict_visit_counts"], {"3": 1, "12": 0, "24": 0, "36": 0, "48": 0, "60": 1})
        self.assertFalse(audit["carry_forward"])

    def test_conflicts_reject_nonselected_or_incomplete_states(self):
        states = {"a": {m: "0" for m in MONTHS}}
        with self.assertRaisesRegex(ValueError, "selected ischemic cohort"):
            apply_discharge_death_conflicts(states, {"a": {"D_DIAG": "2", "D_MRS": "6"}})
        with self.assertRaisesRegex(ValueError, "six validated"):
            apply_discharge_death_conflicts({"a": {3: "0"}}, {"a": {"D_DIAG": "1", "D_MRS": "6"}})


if __name__ == "__main__":
    unittest.main()

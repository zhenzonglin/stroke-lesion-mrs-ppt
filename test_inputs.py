"""合成软件测试，不含患者记录。"""
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch

from inputs import SOURCE_FIELDS, _anonymous_pattern, load_clinical, scan_masks


class Column(list):
    def __init__(self, values, dtype="object"):
        super().__init__(values)
        self.dtype = dtype

    def tolist(self):
        return list(self)


class Frame(dict):
    def __len__(self):
        return len(next(iter(self.values()), []))


def reader(specs, observed):
    def read(path, metadataonly=False, usecols=None, **kwargs):
        source = specs[Path(path).name]
        data = source["data"]
        numeric = source.get("numeric", True)
        types = {field: ("double" if field.casefold() != "code_n" or numeric else "string") for field in data}
        meta = SimpleNamespace(column_names=list(data), readstat_variable_types=types,
                               column_names_to_labels={})
        if metadataonly:
            return None, meta
        observed.extend(usecols)
        return Frame({field: Column(data[field], "float64" if numeric else "object") for field in usecols}), meta
    return SimpleNamespace(read_sas7bdat=read)


class InputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for name in ("early.sas7bdat", "late.sas7bdat"):
            (self.root / name).write_bytes(b"mock_sas")
        self.specs = {
            "early.sas7bdat": {"data": {"code_n": [1., 2.], "F3_MRS": [0., 5.], "F12_MRS": [1., 6.], "D_DIAG": [1., 1.], "F3_DEATH": [2, 1]}},
            "late.sas7bdat": {"data": {"code_n": [1., 3.], "m24_mrs": [2., 6.], "m36_mrs": [1., 6.], "m48_mrs": [0., 6.], "m60_mrs": [0., 6.], "F2Y_DEATH": [2, 1]}},
        }
        self.config = {"source_files": {"early": str(self.root / "early.sas7bdat"), "late": str(self.root / "late.sas7bdat")}}

    def load(self, **kwargs):
        observed = []
        with patch.dict(sys.modules, {"pyreadstat": reader(self.specs, observed)}):
            records, audit = load_clinical(self.config, **kwargs)
        return records, audit, observed

    def mask(self, relative):
        path = self.root / "masks" / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"nonempty_fake_nii")
        return path

    def test_fixed_sources_outer_union_and_death_whitelist(self):
        records, audit, columns = self.load()
        self.assertEqual(set(records), {"1", "2", "3"})
        self.assertEqual(audit["join"], "outer_union")
        self.assertEqual(audit["joined_n"], 3)
        self.assertEqual(audit["early_only_n"], 1)
        self.assertEqual(audit["late_only_n"], 1)
        self.assertEqual(audit["early_missing_late_n"], 1)
        self.assertNotIn("m24_mrs", records["2"])
        self.assertNotIn("F3_MRS", records["3"])
        self.assertNotIn("D_DIAG", records["3"])
        self.assertEqual(records["3"]["m24_mrs"], 6.)
        self.assertTrue(audit["numeric_id_compatible"])
        self.assertTrue(audit["clinical_allowed"])
        self.assertFalse(any("DEATH" in name for name in columns))

    def test_discharge_mrs_is_read_only_when_requested(self):
        original = dict(SOURCE_FIELDS)
        self.specs["early.sas7bdat"]["data"]["d_mrs"] = [2., 6.]
        self.specs["early.sas7bdat"]["data"]["D_DEATH"] = [2, 1]
        records, audit, columns = self.load()
        self.assertTrue(audit["clinical_allowed"])
        self.assertNotIn("D_MRS", records["1"])
        self.assertNotIn("d_mrs", columns)
        records, audit, columns = self.load(include_discharge_mrs=True)
        self.assertTrue(audit["clinical_allowed"])
        self.assertEqual(records["1"]["D_MRS"], 2.)
        self.assertEqual(records["2"]["D_MRS"], 6.)
        self.assertNotIn("D_MRS", records["3"])
        self.assertIn("D_MRS", audit["sources"]["early"]["required_fields"])
        self.assertEqual(audit["sources"]["early"]["schema"]["D_MRS"]["actual"], "d_mrs")
        self.assertNotIn("D_MRS", audit["sources"]["late"]["required_fields"])
        self.assertIn("d_mrs", columns)
        self.assertFalse(any("DEATH" in name.upper() for name in columns))
        self.assertEqual(SOURCE_FIELDS, original)
        _, default_audit, default_columns = self.load()
        self.assertNotIn("D_MRS", default_audit["sources"]["early"]["required_fields"])
        self.assertNotIn("d_mrs", default_columns)

    def test_missing_discharge_mrs_blocks_only_requested_read(self):
        records, audit, _ = self.load(include_discharge_mrs=True)
        self.assertEqual(records, {})
        self.assertFalse(audit["clinical_allowed"])
        self.assertEqual(audit["sources"]["early"]["missing_fields"], ["D_MRS"])
        records, audit, _ = self.load()
        self.assertTrue(audit["clinical_allowed"])
        self.assertEqual(set(records), {"1", "2", "3"})

    def test_discharge_mrs_requires_unique_case_insensitive_column(self):
        self.specs["early.sas7bdat"]["data"].update(D_MRS=[2., 6.], d_mrs=[3., 5.])
        records, audit, columns = self.load(include_discharge_mrs=True)
        self.assertEqual(records, {})
        self.assertFalse(audit["clinical_allowed"])
        self.assertEqual(audit["sources"]["early"]["duplicate_columns"], ["D_MRS"])
        self.assertFalse(any(name.upper() == "D_MRS" for name in columns))
        _, default_audit, _ = self.load()
        self.assertTrue(default_audit["clinical_allowed"])

    def test_duplicate_or_missing_sas_id_blocks(self):
        self.specs["early.sas7bdat"]["data"]["code_n"] = [1., 1.]
        records, audit, _ = self.load()
        self.assertEqual(records, {})
        self.assertFalse(audit["clinical_allowed"])
        self.specs["early.sas7bdat"]["data"]["code_n"] = [1., None]
        records, audit, _ = self.load()
        self.assertFalse(audit["clinical_allowed"])

    def test_character_ids_preserve_zeroes_and_case(self):
        for source in self.specs.values():
            source["numeric"] = False
            source["data"]["code_n"] = ["001", "Ab2"]
        records, audit, _ = self.load()
        self.assertEqual(set(records), {"001", "Ab2"})
        self.assertFalse(audit["numeric_id_compatible"])

    def test_case_insensitive_column_names(self):
        data = self.specs["early.sas7bdat"]["data"]
        data["CODE_N"] = data.pop("code_n")
        data["f3_mrs"] = data.pop("F3_MRS")
        records, audit, _ = self.load()
        self.assertTrue(audit["clinical_allowed"])
        self.assertIn("F3_MRS", records["1"])

    def test_empty_sas_source_is_blocked(self):
        self.specs["early.sas7bdat"]["data"] = {field: [] for field in self.specs["early.sas7bdat"]["data"]}
        records, audit, _ = self.load()
        self.assertFalse(audit["clinical_allowed"])
        self.assertEqual(records, {})
        self.assertEqual(audit["sources"]["early"]["status"], "empty_source")

    def test_numeric_zeroes_and_explicit_width_match(self):
        self.mask("nested/sub-00001_lesion.nii.gz")
        self.mask("00002-mask.nii")
        mapping, audit, inventory = scan_masks(str(self.root / "masks"), {"00001", "00002"})
        self.assertTrue(audit["mapping_allowed"])
        self.assertEqual(set(mapping), {"00001", "00002"})
        self.assertEqual({row["id"] for row in inventory}, {"1", "2"})

    def test_duplicate_unmatched_mask_ids_also_block(self):
        self.mask("1_mask.nii.gz")
        self.mask("nested/099_mask.nii.gz")
        self.mask("99_other.nii.gz")
        mapping, audit, inventory = scan_masks(str(self.root / "masks"), {"1"})
        self.assertFalse(audit["mapping_allowed"])
        self.assertEqual(audit["duplicate_id_n"], 1)

    def test_unknown_pattern_is_not_dropped_or_guessed(self):
        self.mask("1_mask.nii.gz")
        self.mask("patient_Maria_2_lesion.nii.gz")
        _, audit, inventory = scan_masks(str(self.root / "masks"), {"1", "2"})
        self.assertFalse(audit["mapping_allowed"])
        self.assertEqual(len(inventory), 2)
        self.assertNotIn("Maria", str(audit))
        self.assertEqual(audit["unrecognized_n"], 1)

    def test_anonymous_patterns_remove_unicode_names(self):
        self.mask("1_mask.nii.gz")
        self.mask("张三_１２_lesion.nii.gz")
        self.mask("José_éLatin_٢_mask.nii.gz")
        _, audit, inventory = scan_masks(str(self.root / "masks"), {"1"})
        shared = str(audit["filename_patterns"])
        self.assertEqual(audit["unrecognized_n"], 2)
        self.assertFalse(audit["mapping_allowed"])
        self.assertTrue(shared.isascii())
        for secret in ("张三", "１２", "José", "éLatin", "٢"):
            self.assertNotIn(secret, shared)
        self.assertTrue(any("张三" in row["path"] for row in inventory))
        self.assertEqual(_anonymous_pattern("José_éLatin_٢_mask.nii.gz"),
                         "<text>_<text>_<text>_mask.nii.gz")

    def test_explicit_regex_character_ids(self):
        self.mask("patient-Ab02_lesion.nii.gz")
        mapping, audit, _ = scan_masks(str(self.root / "masks"), {"Ab02"},
                                      r"^patient-(?P<id>[A-Za-z0-9]+)_", "character")
        self.assertEqual(set(mapping), {"Ab02"})
        self.assertTrue(audit["mapping_allowed"])

    def test_empty_file_and_unreadable_directory_block(self):
        path = self.mask("1_mask.nii.gz")
        path.write_bytes(b"")
        _, audit, _ = scan_masks(str(self.root / "masks"), {"1"})
        self.assertFalse(audit["mapping_allowed"])
        self.assertEqual(audit["empty_file_n"], 1)

        def fail(root, onerror, followlinks):
            onerror(PermissionError("private_patient_path"))

        with patch("inputs.os.walk", side_effect=fail):
            _, audit, _ = scan_masks(str(self.root / "masks"), {"1"})
        self.assertFalse(audit["mapping_allowed"])
        self.assertIn("directory_scan_PermissionError", audit["errors"])
        self.assertNotIn("private_patient", str(audit))

    def test_clinical_numeric_alias_collision_blocks(self):
        self.mask("1_mask.nii.gz")
        _, audit, _ = scan_masks(str(self.root / "masks"), {"1", "001"})
        self.assertFalse(audit["mapping_allowed"])
        self.assertIn("clinical_id_collision_after_numeric_normalization", audit["errors"])


if __name__ == "__main__":
    unittest.main()

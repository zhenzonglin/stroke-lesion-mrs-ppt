"""只验证算法；所有影像均为临时合成数据，不代表患者结果。"""
import hashlib
from pathlib import Path
import tempfile
import unittest

import nibabel as nib
import numpy as np

from roi_summary import analyse_masks, load_atlas


class ROISummaryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.affine = np.diag([2., 2., 2., 1.])
        self.data = np.zeros((5, 5, 5), dtype=np.int16)
        self.data[:2] = 1
        self.data[3:] = 2
        self.atlas_path = self.save("atlas", self.data)
        self.labels_path = self.root / "labels.tsv"
        self.labels_path.write_text("label_id\troi_name\n0\tBackground\n1\tLeft test\n2\tRight test\n3\tAbsent region\n", encoding="utf-8")
        self.atlas = load_atlas(self.atlas_path, self.labels_path)

    def save(self, name, data, affine=None, unit="mm"):
        image = nib.Nifti1Image(data, self.affine if affine is None else affine)
        image.header.set_xyzt_units(unit)
        path = self.root / f"{name}.nii.gz"
        nib.save(image, path)
        return path

    def mask(self, name, counts):
        values = np.zeros(self.data.shape, dtype=np.uint8)
        for label, count in counts.items():
            coordinates = np.argwhere(self.data == label)[:count]
            values[tuple(coordinates.T)] = 1
        return self.save(name, values)

    def test_counts_threshold_zero_regions_and_denominator(self):
        mapping = {"a": str(self.mask("a", {1: 13})),
                   "b": str(self.mask("b", {1: 7, 2: 15})),
                   "c": str(self.mask("c", {0: 9}))}
        result = analyse_masks(mapping, self.atlas)
        self.assertEqual(result["valid_n"], 3)
        left, right, absent = result["roi_rows"]
        self.assertEqual((left["n"], left["any_overlap_n"], left["denominator_n"]), (1, 2, 3))
        self.assertEqual((right["n"], right["any_overlap_n"]), (1, 1))
        self.assertEqual((absent["n"], absent["roi_voxels"]), (0, 0))
        self.assertAlmostEqual(left["total_overlap_ml"], 20 * .008)
        self.assertEqual(result["coverage_summary"]["no_atlas_overlap_n"], 1)
        self.assertEqual(result["coverage_summary"]["covered_voxels_total"], 35)
        self.assertEqual(result["coverage_summary"]["lesion_voxels_total"], 44)
        self.assertTrue(result["validation"]["input_conservation"])
        self.assertTrue(result["validation"]["all_roi_counts_within_denominator"])
        self.assertEqual(result["validation"]["resampled_or_registered_masks_n"], 0)

    def test_zero_threshold_still_requires_positive_overlap(self):
        result = analyse_masks({"a": str(self.mask("a", {1: 1}))}, self.atlas, 0)
        self.assertEqual([row["n"] for row in result["roi_rows"]], [1, 0, 0])

    def test_one_person_can_belong_to_two_rois(self):
        result = analyse_masks({"a": str(self.mask("a", {1: 13, 2: 13}))}, self.atlas)
        self.assertEqual(sum(row["n"] for row in result["roi_rows"]), 2)
        self.assertEqual(result["valid_n"], 1)

    def test_flip_axes_preserves_anatomical_side(self):
        raw = np.zeros(self.data.shape, dtype=np.uint8)
        raw[0, 1, 1] = 1
        orientation = np.eye(4)
        orientation[0, 0] = -1
        orientation[0, 3] = self.data.shape[0] - 1
        flipped = self.save("flipped", raw[::-1], self.affine @ orientation)
        result = analyse_masks({"a": str(flipped)}, self.atlas, 0)
        self.assertEqual([row["n"] for row in result["roi_rows"]], [1, 0, 0])
        self.assertEqual(result["orientation_reordered_n"], 1)

    def test_mismatched_grid_is_excluded_without_interpolation(self):
        shifted = self.affine.copy()
        shifted[0, 3] += 1
        path = self.save("shifted", np.ones(self.data.shape, dtype=np.uint8), shifted)
        result = analyse_masks({"a": str(path)}, self.atlas)
        self.assertEqual(result["grid_mismatch_n"], 1)
        self.assertEqual(result["valid_n"], 0)
        self.assertTrue(all(row["proportion"] is None for row in result["roi_rows"]))

    def test_empty_nonbinary_nan_unknown_units_and_qform_conflict(self):
        empty = self.save("empty", np.zeros(self.data.shape, dtype=np.uint8))
        fractional = self.save("fractional", np.full(self.data.shape, .5, dtype=np.float32))
        nonfinite = self.save("nonfinite", np.full(self.data.shape, np.nan, dtype=np.float32))
        unknown = self.save("unknown", np.ones(self.data.shape, dtype=np.uint8), unit="unknown")
        image = nib.Nifti1Image(np.ones(self.data.shape, dtype=np.uint8), self.affine)
        shifted = self.affine.copy()
        shifted[0, 3] = 8
        image.set_qform(shifted, code=1)
        image.set_sform(self.affine, code=2)
        image.header.set_xyzt_units("mm")
        conflict = self.root / "conflict.nii.gz"
        nib.save(image, conflict)
        mapping = {key: str(value) for key, value in locals().copy().items()
                   if key in {"empty", "fractional", "nonfinite", "unknown", "conflict"}}
        result = analyse_masks(mapping, self.atlas)
        self.assertEqual(result["valid_n"], 0)
        self.assertEqual(result["excluded_n"], 5)
        self.assertEqual(result["empty_n"], 1)
        self.assertEqual(result["invalid_n"], 4)
        self.assertEqual(result["exclusion_counts"]["qform_sform_conflict"], 1)

    def test_unknown_unit_requires_explicit_opt_in_and_exact_grid(self):
        values = np.zeros(self.data.shape, dtype=np.uint8)
        values[0, 0, 0] = 1
        path = self.save("unknown_exact", values, unit="unknown")
        rejected = analyse_masks({"a": str(path)}, self.atlas, 0)
        self.assertEqual(rejected["valid_n"], 0)
        accepted = analyse_masks({"a": str(path)}, self.atlas, 0,
                                 allow_unknown_units_if_exact_grid=True)
        self.assertEqual(accepted["valid_n"], 1)
        self.assertEqual(accepted["units_assumed_mm_n"], 1)
        self.assertEqual(accepted["roi_rows"][0]["n"], 1)

    def test_unknown_mismatched_grid_and_declared_other_units_never_assumed_mm(self):
        values = np.ones(self.data.shape, dtype=np.uint8)
        shifted = self.affine.copy()
        shifted[0, 3] += 1
        unknown_shifted = self.save("unknown_shifted", values, shifted, unit="unknown")
        meter = self.save("meter", values, unit="meter")
        micron = self.save("micron", values, unit="micron")
        result = analyse_masks({"a": str(unknown_shifted), "b": str(meter), "c": str(micron)},
                               self.atlas, allow_unknown_units_if_exact_grid=True)
        self.assertEqual(result["valid_n"], 0)
        self.assertEqual(result["units_assumed_mm_n"], 0)
        self.assertEqual(result["exclusion_counts"]["grid_mismatch"], 1)
        self.assertEqual(result["exclusion_counts"]["spatial_unit_not_mm"], 2)

    def test_atlas_missing_label_name_or_fractional_values_are_rejected(self):
        self.labels_path.write_text("label_id\troi_name\n1\tOnly one\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "lack TSV names"):
            load_atlas(self.atlas_path, self.labels_path)
        fractional = self.save("fractional_atlas", np.full(self.data.shape, 1.2, dtype=np.float32))
        with self.assertRaisesRegex(ValueError, "integer labels"):
            load_atlas(fractional, self.labels_path)

    def test_all_input_files_remain_unchanged_and_progress_final(self):
        path = self.mask("a", {1: 13})
        paths = [self.atlas_path, self.labels_path, path]
        before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
        calls = []
        analyse_masks({"SECRET_PATIENT": str(path)}, self.atlas, progress=lambda n, total: calls.append((n, total)))
        self.assertEqual(calls, [(1, 1)])
        self.assertEqual(before, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

    def test_independent_channels_keep_overlap_but_coverage_counts_union_once(self):
        channels = np.zeros((*self.data.shape, 3), dtype=np.uint8)
        # 第一个ROI代表已合并的基底节；同一受试者/体素都只计一次。
        channels[:2, :, :, 0] = 1
        channels[1:3, :, :, 1] = 1
        path = self.save("independent", channels)
        self.labels_path.write_text("label_id\troi_name\n1\tBasal Ganglia\n2\tInternal Capsule\n3\tEmpty region\n", encoding="utf-8")
        atlas = load_atlas(path, self.labels_path, "independent_binary_channels")
        values = np.zeros(self.data.shape, dtype=np.uint8)
        values[1, :, :] = 1  # 25个共享体素，不得算成50个被覆盖体素。
        lesion = self.save("overlap", values)
        result = analyse_masks({"a": str(lesion)}, atlas)
        self.assertEqual([row["any_overlap_n"] for row in result["roi_rows"]], [1, 1, 0])
        self.assertEqual([row["n"] for row in result["roi_rows"]], [1, 1, 0])
        self.assertEqual(result["coverage_summary"]["covered_voxels_total"], 25)
        self.assertEqual(result["coverage_summary"]["pooled_covered_fraction"], 1)
        self.assertEqual(atlas["overlapping_atlas_voxels_n"], 25)
        self.assertEqual(atlas["shape"], [5, 5, 5])
        self.assertEqual(atlas["source_shape"], [5, 5, 5, 3])
        self.assertAlmostEqual(result["roi_rows"][0]["total_overlap_ml"], 25*.008)
        self.assertEqual(result["private_any_memberships"]["roi_0001"], {"a"})

    def test_independent_channel_atlas_requires_explicit_mode_binary_and_exact_labels(self):
        channels = np.zeros((*self.data.shape, 2), dtype=np.uint8)
        channels[0, 0, 0, :] = 1
        path = self.save("independent_invalid", channels)
        with self.assertRaisesRegex(ValueError, "not_3d"):
            load_atlas(path, self.labels_path)
        with self.assertRaisesRegex(ValueError, "match channels"):
            load_atlas(path, self.labels_path, "independent_binary_channels")
        self.labels_path.write_text("label_id\troi_name\n1\tA\n2\tB\n", encoding="utf-8")
        channels[0, 0, 0, 1] = 2
        path = self.save("nonbinary_channels", channels)
        with self.assertRaisesRegex(ValueError, "binary"):
            load_atlas(path, self.labels_path, "independent_binary_channels")

    def test_independent_channels_axis_flip_preserves_channel_identity(self):
        channels = np.zeros((*self.data.shape, 2), dtype=np.uint8)
        channels[0, :, :, 0] = 1
        channels[-1, :, :, 1] = 1
        flip = np.eye(4)
        flip[0, 0] = -1
        flip[0, 3] = self.data.shape[0]-1
        path = self.save("independent_flipped", channels[::-1], self.affine @ flip)
        self.labels_path.write_text("label_id\troi_name\n1\tA\n2\tB\n", encoding="utf-8")
        atlas = load_atlas(path, self.labels_path, "independent_binary_channels")
        values = np.zeros(self.data.shape, dtype=np.uint8)
        values[0, 0, 0] = 1
        result = analyse_masks({"a": str(self.save("side", values))}, atlas, 0)
        self.assertEqual([row["n"] for row in result["roi_rows"]], [1, 0])


if __name__ == "__main__":
    unittest.main()

"""MNI病灶与离散ROI的只读重叠统计；不做配准或插值重采样。

所有 private_* 值只供进程内连接，保存共享报告前必须移除。
"""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import math
from pathlib import Path
from typing import Callable

import nibabel as nib
import numpy as np

AFFINE_ATOL = 1e-4


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _header_problem(image: nib.Nifti1Image, allow_unknown_unit: bool = False,
                    expected_dimensions: int = 3) -> str | None:
    if len(image.shape) != expected_dimensions:
        return f"not_{expected_dimensions}d"
    unit = image.header.get_xyzt_units()[0]
    if unit != "mm" and not (allow_unknown_unit and unit == "unknown"):
        return "spatial_unit_not_mm"
    affine = np.asarray(image.affine, dtype=float)
    if affine.shape != (4, 4) or not np.isfinite(affine).all():
        return "nonfinite_affine"
    if not np.allclose(affine[3], [0, 0, 0, 1], atol=AFFINE_ATOL, rtol=0):
        return "invalid_affine_last_row"
    if abs(float(np.linalg.det(affine[:3, :3]))) <= 1e-12:
        return "singular_affine"
    qform, qcode = image.get_qform(coded=True)
    sform, scode = image.get_sform(coded=True)
    if not qcode and not scode:
        return "qform_sform_both_undefined"
    if qcode and scode and not np.allclose(qform, sform, atol=AFFINE_ATOL, rtol=0):
        return "qform_sform_conflict"
    return None


def load_atlas(atlas_path: str | Path, labels_path: str | Path,
               atlas_format: str = "discrete_labels") -> dict:
    """读取显式指定的3D标签或4D独立二值ROI，不按表行号猜标签。

    独立ROI的TSV label_id必须恰为1..通道数。跨ROI体素允许重叠，
    每个ROI内部为二值并集，整体覆盖率另按所有通道的并集计算。
    """
    if atlas_format not in {"discrete_labels", "independent_binary_channels"}:
        raise ValueError("Unsupported atlas_format")
    independent = atlas_format == "independent_binary_channels"
    atlas_path, labels_path = Path(atlas_path), Path(labels_path)
    image = nib.load(str(atlas_path))
    problem = _header_problem(image, expected_dimensions=4 if independent else 3)
    if problem:
        raise ValueError("Invalid atlas header: " + problem)
    source_shape = list(image.shape)
    source_affine = np.asarray(image.affine).tolist()
    # 仅交换/翻转数组轴，保留原来物理坐标；不是重新配准或重采样。
    image = nib.as_closest_canonical(image)
    values = np.asanyarray(image.dataobj)
    if not np.isfinite(values).all() or np.any(values < 0) or not np.all(values == np.rint(values)):
        raise ValueError("Atlas must contain finite nonnegative integer labels")
    if float(values.max(initial=0)) > np.iinfo(np.int32).max:
        raise ValueError("Atlas label values exceed int32 range")
    if independent and not np.all((values == 0) | (values == 1)):
        raise ValueError("Independent atlas channels must be binary 0/1")
    channel_masks = values.astype(bool, copy=False) if independent else None
    data = np.any(channel_masks, axis=3) if independent else values.astype(np.int32, copy=False)
    names = {}
    with labels_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"label_id", "roi_name"}.issubset(reader.fieldnames or []):
            raise ValueError("ROI TSV requires label_id and roi_name")
        for row in reader:
            raw_id = row["label_id"].strip()
            try:
                label_id = int(raw_id)
            except (TypeError, ValueError) as exc:
                raise ValueError("ROI TSV label_id must be an integer") from exc
            if label_id < 0 or str(label_id) != raw_id:
                raise ValueError("ROI TSV label_id must be a canonical nonnegative integer")
            if label_id in names:
                raise ValueError("Duplicate ROI label_id in TSV")
            name = row["roi_name"].strip()
            if not name:
                raise ValueError("Empty ROI name in TSV")
            names[label_id] = name
    if independent:
        if set(names).difference({0}) != set(range(1, image.shape[3]+1)):
            raise ValueError("Independent ROI TSV labels must match channels 1..N exactly")
        voxel_counts = {index+1: int(np.count_nonzero(channel_masks[..., index]))
                        for index in range(image.shape[3])}
    else:
        present_labels, counts = np.unique(data, return_counts=True)
        voxel_counts = {int(label): int(count) for label, count in zip(present_labels, counts)}
    missing_names = set(voxel_counts).difference({0}).difference(names)
    if missing_names:
        raise ValueError("Positive atlas labels lack TSV names: " + ",".join(map(str, sorted(missing_names))))
    if not np.any(data):
        raise ValueError("Atlas has no positive labels")
    voxel_volume_ml = abs(float(np.linalg.det(image.affine[:3, :3]))) / 1000.0
    labels = [{"label_id": label_id, "roi_id": f"roi_{label_id:04d}",
               "roi_name": names[label_id], "roi_voxels": voxel_counts.get(label_id, 0),
               "roi_volume_ml": voxel_counts.get(label_id, 0) * voxel_volume_ml}
              for label_id in sorted(names) if label_id > 0]
    return {
        "atlas_path": str(atlas_path.resolve()), "labels_path": str(labels_path.resolve()),
        "sha256": _sha256(atlas_path), "labels_sha256": _sha256(labels_path),
        "shape": list(image.shape[:3]), "affine": image.affine.tolist(),
        "source_shape": source_shape, "source_affine": source_affine,
        "axis_codes": list(nib.aff2axcodes(image.affine)), "spatial_unit": "mm",
        "voxel_volume_ml": voxel_volume_ml, "labels": labels,
        "atlas_format": atlas_format,
        "atlas_interpretation": ("Independent binary ROIs may overlap across channels; coverage uses their union."
                                 if independent else "Discrete label atlas.")
                                + " Geometric compatibility does not verify anatomical registration.",
        "overlapping_atlas_voxels_n": (int(np.count_nonzero(np.sum(channel_masks, axis=3) > 1))
                                       if independent else 0),
        "private_data": data, "private_channel_masks": channel_masks, "private_image": image,
    }


def analyse_masks(mapping: dict[str, str], atlas: dict, min_overlap_ml: float = 0.1,
                  progress: Callable[[int, int], None] | None = None,
                  allow_unknown_units_if_exact_grid: bool = False) -> dict:
    """每个ROI分母均为有效非空病灶人数；ROI人群可相互重叠。

    min_overlap_ml 是可配置的展示噪声阈值，不是临床或显著性阈值。
    any_overlap_n 同时保留任意正重叠定义。空掩膜不当作无病灶病例。
    单位缺失仅在显式允许且与毫米图谱网格完全相容时解释为mm；这不是配准验证。
    """
    if (isinstance(min_overlap_ml, bool) or not math.isfinite(min_overlap_ml)
            or min_overlap_ml < 0):
        raise ValueError("min_overlap_ml must be a finite nonnegative number")
    if any(not isinstance(key, str) or not key for key in mapping):
        raise ValueError("mapping must use nonempty unique string IDs")
    if not isinstance(allow_unknown_units_if_exact_grid, bool):
        raise ValueError("allow_unknown_units_if_exact_grid must be a bool")
    atlas_data = atlas["private_data"]
    channel_masks = atlas.get("private_channel_masks")
    atlas_affine = np.asarray(atlas["affine"], dtype=float)
    atlas_shape = tuple(atlas["shape"])
    voxel_ml = float(atlas["voxel_volume_ml"])
    labels = atlas["labels"]
    memberships = {row["roi_id"]: set() for row in labels}
    any_memberships = {row["roi_id"]: set() for row in labels}
    totals = Counter()
    valid_ids, exclusions = set(), {}
    coverage, lesion_volumes = [], []
    covered_voxels_total, lesion_voxels_total = 0, 0
    reordered_n = 0
    units_assumed_mm_n = 0
    no_coverage_n = 0
    partial_coverage_n = 0
    for processed, (identifier, path) in enumerate(sorted(mapping.items()), 1):
        try:
            image = nib.load(str(path))
            unknown_unit = image.header.get_xyzt_units()[0] == "unknown"
            problem = _header_problem(image, allow_unknown_unit=allow_unknown_units_if_exact_grid)
            if problem:
                raise ValueError(problem)
            source_affine = np.asarray(image.affine)
            source_shape = tuple(image.shape)
            image = nib.as_closest_canonical(image)
            if tuple(image.shape) != atlas_shape or not np.allclose(image.affine, atlas_affine, atol=AFFINE_ATOL, rtol=0):
                raise ValueError("grid_mismatch")
            values = np.asanyarray(image.dataobj)
            if not np.isfinite(values).all():
                raise ValueError("nonfinite_mask")
            if not np.all((values == 0) | (values == 1)):
                raise ValueError("not_binary_0_1")
            lesion = values == 1
            lesion_count = int(np.count_nonzero(lesion))
            if lesion_count == 0:
                raise ValueError("empty_mask")
            # 只有通过网格与数据检查、真正进入分母的病例才计入单位解释人数。
            units_assumed_mm_n += int(unknown_unit)
            if source_shape != tuple(image.shape) or not np.allclose(source_affine, image.affine, atol=AFFINE_ATOL, rtol=0):
                reordered_n += 1
            if channel_masks is not None:
                overlap_counts = np.count_nonzero(channel_masks[lesion], axis=0)
                overlaps = {index+1: int(count) for index, count in enumerate(overlap_counts)}
                covered = int(np.count_nonzero(atlas_data[lesion]))
            else:
                overlap_labels, overlap_counts = np.unique(atlas_data[lesion], return_counts=True)
                overlaps = {int(k): int(v) for k, v in zip(overlap_labels, overlap_counts)}
                covered = lesion_count - overlaps.get(0, 0)
            no_coverage_n += int(covered == 0)
            partial_coverage_n += int(0 < covered < lesion_count)
            coverage.append(covered / lesion_count)
            lesion_volumes.append(lesion_count * voxel_ml)
            covered_voxels_total += covered
            lesion_voxels_total += lesion_count
            valid_ids.add(identifier)
            for row in labels:
                count = overlaps.get(row["label_id"], 0)
                volume = count * voxel_ml
                totals[row["label_id"]] += count
                if count:
                    any_memberships[row["roi_id"]].add(identifier)
                    if volume + 1e-12 >= min_overlap_ml:
                        memberships[row["roi_id"]].add(identifier)
        except ValueError as exc:
            # 仅内部已知错误码可以进入共享报告，不能输出可能带患者路径的异常。
            allowed = {"not_3d", "spatial_unit_not_mm", "nonfinite_affine", "singular_affine",
                       "invalid_affine_last_row", "qform_sform_both_undefined", "qform_sform_conflict",
                       "grid_mismatch", "nonfinite_mask", "not_binary_0_1", "empty_mask"}
            exclusions[identifier] = str(exc) if str(exc) in allowed else "read_error_ValueError"
        except Exception as exc:
            exclusions[identifier] = "read_error_" + type(exc).__name__
        if progress is not None and (processed % 100 == 0 or processed == len(mapping)):
            progress(processed, len(mapping))
    n = len(valid_ids)
    roi_rows = [{**row, "n": len(memberships[row["roi_id"]]), "denominator_n": n,
                 "proportion": len(memberships[row["roi_id"]]) / n if n else None,
                 "any_overlap_n": len(any_memberships[row["roi_id"]]),
                 "any_overlap_proportion": len(any_memberships[row["roi_id"]]) / n if n else None,
                 "total_overlap_ml": totals[row["label_id"]] * voxel_ml}
                for row in labels]
    exclusion_counts = dict(sorted(Counter(exclusions.values()).items()))
    coverage_summary = {
        "denominator_n": n, "no_atlas_overlap_n": no_coverage_n,
        "partial_atlas_overlap_n": partial_coverage_n,
        "complete_atlas_overlap_n": n - no_coverage_n - partial_coverage_n,
        "covered_voxels_total": covered_voxels_total, "lesion_voxels_total": lesion_voxels_total,
        "pooled_covered_fraction": covered_voxels_total / lesion_voxels_total if lesion_voxels_total else None,
        "coverage_rule": "union_of_all_positive_ROI_voxels; overlapping_ROI_voxels_count_once",
    }
    for name, values in (("patient_coverage_fraction", coverage), ("lesion_volume_ml", lesion_volumes)):
        coverage_summary[name] = ({"min": float(np.min(values)), "q25": float(np.quantile(values, .25)),
                                   "median": float(np.median(values)), "q75": float(np.quantile(values, .75)),
                                   "max": float(np.max(values)), "mean": float(np.mean(values))} if values else None)
    return {
        "input_n": len(mapping), "valid_n": n, "excluded_n": len(exclusions),
        "empty_n": exclusion_counts.get("empty_mask", 0),
        "grid_mismatch_n": exclusion_counts.get("grid_mismatch", 0),
        "invalid_n": len(exclusions) - exclusion_counts.get("empty_mask", 0) - exclusion_counts.get("grid_mismatch", 0),
        "exclusion_counts": exclusion_counts, "orientation_reordered_n": reordered_n,
        "units_assumed_mm_n": units_assumed_mm_n,
        "allow_unknown_units_if_exact_grid": allow_unknown_units_if_exact_grid,
        "units_assumption_rule": ("Unknown unit interpreted as mm only after exact atlas grid match; not anatomical registration verification."
                                  if allow_unknown_units_if_exact_grid else "Unknown spatial units are rejected."),
        "min_overlap_ml": min_overlap_ml, "membership_rule": "positive_overlap_and_volume_at_least_min_overlap_ml",
        "denominator_rule": "all_valid_nonempty_masks_supplied_by_caller",
        "threshold_interpretation": "Configurable display threshold; not a clinical threshold; any-overlap sensitivity retained.",
        "roi_rows": roi_rows, "coverage_summary": coverage_summary,
        "validation": {"input_conservation": len(mapping) == n + len(exclusions),
                       "all_roi_counts_within_denominator": all(0 <= row["n"] <= row["any_overlap_n"] <= n for row in roi_rows),
                       "roi_people_overlap_allowed": True, "resampled_or_registered_masks_n": 0},
        "private_valid_ids": valid_ids, "private_memberships": memberships,
        "private_any_memberships": any_memberships, "private_exclusions": exclusions,
    }

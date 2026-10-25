#!/usr/bin/env python3
"""从既有运行补算基底节-内囊患者并集；不重读SAS，不覆盖旧结果。"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

from mrs_summary import MONTHS, STATES, aggregate_groups
from roi_summary import analyse_masks, load_atlas
from run_summary import csv_rows, digest, write_json

MERGED_ID = "basal_ganglia_internal_capsule"
MERGED_NAME = "Basal Ganglia-Internal Capsule"
MERGED_NAME_ZH = "基底节-内囊"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def load_original_states(path: Path) -> dict[str, dict[int, str]]:
    states: dict[str, dict[int, str]] = {}
    for row in read_csv(path):
        identifier = row.get("patient_id", "")
        if not identifier or identifier in states:
            raise ValueError("Original derived states have missing or duplicate IDs")
        visits = {month: row.get(f"m{month}", "") for month in MONTHS}
        if any(value not in STATES for value in visits.values()):
            raise ValueError("Original derived states contain unsupported visit states")
        states[identifier] = visits
    if not states:
        raise ValueError("Original derived states are empty")
    # 原文件保存的是已处理死亡冲突的状态，不重新解释或插补。
    aggregate_groups(states, {})
    return states


def original_mapping(path: Path, expected_n: int) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for row in read_csv(path):
        if row.get("status") != "matched":
            continue
        identifier = row.get("clinical_id", "")
        source = Path(row.get("path", ""))
        if not identifier or identifier in mapping:
            raise ValueError("Original matched inventory has missing or duplicate clinical IDs")
        if not source.is_file():
            raise ValueError("An original matched mask is no longer available")
        try:
            expected_bytes = int(row["bytes"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Original inventory lacks valid file sizes") from exc
        if source.stat().st_size != expected_bytes:
            raise ValueError("A mask size changed since the original inventory; rerun the full workflow")
        mapping[identifier] = str(source)
    if len(mapping) != expected_n:
        raise ValueError("Original matched inventory count differs from the report")
    return mapping


def distribution_counts(rows: list[dict], group: str) -> dict[tuple[int, str], int]:
    selected: dict[tuple[int, str], int] = {}
    for row in rows:
        if row["group"] != group:
            continue
        key = (int(row["month"]), str(row["state"]))
        if key in selected:
            raise ValueError("Original distributions contain duplicate cells")
        selected[key] = int(row["n"])
    if set(selected) != {(month, state) for month in MONTHS for state in STATES}:
        raise ValueError("Original distributions must contain all six visits and all states")
    return selected


def component_ids(report: dict) -> tuple[str, str]:
    identifiers = []
    for name in ("Basal Ganglia", "Internal Capsule"):
        matches = [row["roi_id"] for row in report["roi"]["roi_rows"] if row["roi_name"] == name]
        if len(matches) != 1:
            raise ValueError("The source run must contain one Basal Ganglia and one Internal Capsule ROI")
        identifiers.append(matches[0])
    return identifiers[0], identifiers[1]


def targeted_atlas(atlas: dict, bg_id: str, ic_id: str) -> dict:
    return union_atlas(atlas, [bg_id, ic_id], MERGED_ID, MERGED_NAME)


def union_atlas(atlas: dict, component_roi_ids: list[str], union_id: str, union_name: str) -> dict:
    """选取组成区并增加二值OR通道；不改变原图谱或原label_id含义。"""
    if atlas.get("atlas_format") != "independent_binary_channels":
        raise ValueError("The supplement requires the v2 independent binary ROI atlas")
    by_id = {row["roi_id"]: row for row in atlas["labels"]}
    if not component_roi_ids or len(set(component_roi_ids)) != len(component_roi_ids):
        raise ValueError("Union requires distinct component ROI IDs")
    if union_id in by_id or any(key not in by_id for key in component_roi_ids):
        raise ValueError("Component ROI IDs are absent from the frozen atlas")
    components = [(key, by_id[key]["roi_name"], atlas["private_channel_masks"][..., by_id[key]["label_id"]-1])
                  for key in component_roi_ids]
    union = np.logical_or.reduce([row[2] for row in components])
    components.append((union_id, union_name, union))
    masks = np.stack([row[2] for row in components], axis=3)
    labels = []
    for index, (roi_id, name, values) in enumerate(components, 1):
        count = int(np.count_nonzero(values))
        labels.append({"label_id": index, "roi_id": roi_id, "roi_name": name,
                       "roi_voxels": count, "roi_volume_ml": count*atlas["voxel_volume_ml"]})
    return {**atlas, "private_channel_masks": masks, "private_data": union, "labels": labels}


def _atlas_file(report: dict, key: str, config_key: str, repo: Path) -> Path:
    recorded = Path(report["atlas"][key])
    fallback = Path(report["config"][config_key])
    if not fallback.is_absolute():
        fallback = repo/fallback
    for candidate in (recorded, fallback):
        if candidate.is_file():
            return candidate
    raise ValueError("The original frozen atlas or label file is unavailable")


def load_reproduction_context(run_dir: Path) -> tuple:
    """只读加载原运行的派生状态与匹配清单，验证冻结图谱。返回值仅限进程内。"""
    run_dir = run_dir.resolve()
    report_path = run_dir/"share/ppt_data.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("status") not in {"COMPLETED", "COMPLETED_WITH_FLAGS"}:
        raise ValueError("Only a completed original run can be supplemented")
    cfg = report["config"]
    if (cfg.get("death_definition") != "mrs_equal_6_only_at_all_visits"
            or cfg.get("carry_forward") is not False
            or cfg.get("primary_membership") != "any_positive_overlap"):
        raise ValueError("The original run uses incompatible death or ROI definitions")
    original_n = int(report["roi"]["valid_n"])
    states_path = run_dir/"private/derived_states.csv"
    inventory_path = run_dir/"private/mask_inventory.csv"
    exclusions_path = run_dir/"private/mask_exclusions.csv"
    states = load_original_states(states_path)
    if len(states) != original_n:
        raise ValueError("Original state IDs do not match the valid cohort size")
    mapping = original_mapping(inventory_path, int(report["roi"]["input_n"]))
    exclusions_rows = read_csv(exclusions_path)
    old_exclusions = {row["patient_id"]: row["reason"] for row in exclusions_rows}
    if len(old_exclusions) != len(exclusions_rows) or len(old_exclusions) != report["roi"]["excluded_n"]:
        raise ValueError("Original exclusions contain duplicates or inconsistent counts")
    if set(states) & set(old_exclusions) or set(mapping) != set(states) | set(old_exclusions):
        raise ValueError("Original matched IDs do not partition into valid and excluded subjects")
    repo = Path(__file__).resolve().parent
    atlas = load_atlas(_atlas_file(report, "atlas_path", "atlas_path", repo),
                       _atlas_file(report, "labels_path", "atlas_labels", repo),
                       report["atlas"].get("atlas_format", "discrete_labels"))
    if (atlas["sha256"] != report["atlas"]["sha256"]
            or atlas["labels_sha256"] != report["atlas"]["labels_sha256"]):
        raise ValueError("Frozen atlas checksums differ from the original run")
    return report, states, mapping, old_exclusions, atlas, report_path, states_path, inventory_path


def verify_reproduction(report: dict, states: dict, old_exclusions: dict, roi: dict,
                        component_roi_ids: list[str]) -> tuple[dict, dict, dict, dict]:
    """复现原有效集合、排除原因和全部组成区两种阈值的完整mRS边际。"""
    if roi["private_valid_ids"] != set(states):
        raise ValueError("Rechecked valid patient set differs from the original run; rerun the full workflow")
    if roi["private_exclusions"] != old_exclusions:
        raise ValueError("Rechecked exclusion decisions differ from the original run")
    primary = aggregate_groups(states, roi["private_any_memberships"])
    sensitivity = aggregate_groups(states, roi["private_memberships"])
    original_rows = {row["roi_id"]: row for row in report["roi"]["roi_rows"]}
    current_rows = {row["roi_id"]: row for row in roi["roi_rows"]}
    for key in component_roi_ids:
        if current_rows[key]["any_overlap_n"] != original_rows[key]["any_overlap_n"]:
            raise ValueError("Component any-overlap count differs from the original report")
        if current_rows[key]["n"] != original_rows[key]["n"]:
            raise ValueError("Component sensitivity count differs from the original report")
    for version, current in (("primary", primary), ("sensitivity", sensitivity)):
        for key in ("overall", *component_roi_ids):
            if distribution_counts(current["distributions"], key) != distribution_counts(report[version]["distributions"], key):
                raise ValueError("Overall or component six-visit mRS distributions do not reproduce the source run")
    return primary, sensitivity, original_rows, current_rows


def summarize(run_dir: Path, progress=None) -> dict:
    """全部复现检查通过后才返回分享结果，异常中不包含患者标识。"""
    report, states, mapping, old_exclusions, atlas, report_path, states_path, inventory_path = load_reproduction_context(run_dir)
    cfg = report["config"]
    original_n = int(report["roi"]["valid_n"])
    bg_id, ic_id = component_ids(report)
    restricted = targeted_atlas(atlas, bg_id, ic_id)
    roi = analyse_masks(mapping, restricted, float(cfg["sensitivity_min_overlap_ml"]), progress=progress,
                        allow_unknown_units_if_exact_grid=cfg["allow_unknown_units_if_exact_grid"])
    primary, sensitivity, original_rows, current_rows = verify_reproduction(report, states, old_exclusions, roi, [bg_id, ic_id])
    bg = roi["private_any_memberships"][bg_id]
    ic = roi["private_any_memberships"][ic_id]
    union = roi["private_any_memberships"][MERGED_ID]
    intersection = bg & ic
    if union != bg | ic or len(union) != len(bg)+len(ic)-len(intersection):
        raise RuntimeError("Patient union inclusion-exclusion validation failed")
    merged = {**current_rows[MERGED_ID], "label_id": None, "roi_name_zh": MERGED_NAME_ZH,
              "component_roi_ids": [bg_id, ic_id],
              "component_label_ids": [original_rows[bg_id]["label_id"], original_rows[ic_id]["label_id"]],
              "min_overlap_ml": float(cfg["sensitivity_min_overlap_ml"]),
              "component_rule": "binary_voxel_union; one patient counted at most once"}
    ranked = [{**row} for row in report["roi"]["roi_rows"] if row["roi_id"] not in {bg_id, ic_id}]+[merged]
    ranked.sort(key=lambda row: (-row["any_overlap_n"], row["roi_name"]))
    for rank, row in enumerate(ranked, 1):
        row["rank"] = rank
    row_lookup = {row["roi_id"]: row for row in ranked}
    old_outcomes = [row for row in report["primary"]["outcome_rates"]
                    if row["group"] not in {"overall", bg_id, ic_id} and row["definition"] == "mrs_2_6"]
    merged_outcomes = [row for row in primary["outcome_rates"] if row["group"] == MERGED_ID]
    combined_outcomes = old_outcomes+[row for row in merged_outcomes if row["definition"] == "mrs_2_6"]
    ranked_outcomes = [{**row, "rank": row_lookup[row["group"]]["rank"],
                       "roi_name": row_lookup[row["group"]]["roi_name"],
                       "lesion_n": row_lookup[row["group"]]["any_overlap_n"],
                       "lesion_denominator_n": original_n,
                       "lesion_proportion": row_lookup[row["group"]]["any_overlap_n"]/original_n}
                      for row in combined_outcomes]
    ranked_outcomes.sort(key=lambda row: (row["month"], row["rank"]))
    expected_pairs = {(roi_id, month) for roi_id in row_lookup for month in MONTHS}
    actual_pairs = {(row["group"], row["month"]) for row in ranked_outcomes}
    if len(ranked_outcomes) != len(expected_pairs) or actual_pairs != expected_pairs:
        raise ValueError("The source report lacks complete ROI outcome rows")
    return {"status": "COMPLETED", "schema_version": "bgic_union_v1", "source_run_id": report["run_id"],
            "source_report_sha256": digest(report_path), "source_derived_states_sha256": digest(states_path),
            "source_inventory_sha256": digest(inventory_path), "cohort_n": original_n,
            "atlas_sha256": atlas["sha256"], "labels_sha256": atlas["labels_sha256"],
            "months": list(MONTHS), "merged_roi": merged,
            "inclusion_exclusion": {"basal_ganglia_n": len(bg), "internal_capsule_n": len(ic),
                "intersection_n": len(intersection), "union_n": len(union),
                "formula": "union_n = basal_ganglia_n + internal_capsule_n - intersection_n"},
            "merged_distributions": [row for row in primary["distributions"] if row["group"] == MERGED_ID],
            "merged_recovery_distributions": [row for row in primary["recovery_distributions"] if row["group"] == MERGED_ID],
            "merged_outcome_rates": merged_outcomes,
            "merged_sensitivity_distributions": [row for row in sensitivity["distributions"] if row["group"] == MERGED_ID],
            "ranked_rois": ranked, "ranked_outcome_rates": ranked_outcomes,
            "validation": {"all_passed": True, "valid_patient_set_reproduced": True,
                "exclusions_reproduced": True, "component_counts_reproduced": True,
                "component_and_overall_mrs_cells_checked": 2*3*len(MONTHS)*len(STATES),
                "patient_union_verified": True, "source_mask_sizes_reproduced": True,
                "original_mask_sha256_available": False},
            "definitions": {"membership": "any positive overlap with the voxel union of Basal Ganglia and Internal Capsule",
                "sensitivity": ">=0.1mL default, applied to union voxels once; not the union of separate threshold memberships",
                "outcome": "mRS2-6 numerator / valid observed mRS0-6 at each visit; mRS6 is death; no carry-forward",
                "ranking": "descending any-overlap patient count; ties resolved by ROI name; common cohort denominator",
                "unchanged_roi_values": "Other ROI aggregates are retained from the original run; their masks are not re-counted here",
                "reproducibility_limit": "Original per-mask SHA256 fingerprints were not recorded. File sizes, exact valid ID set, exclusions and component marginals reproduce; byte-level image immutability is not established."}}


def display_pages(result: dict) -> list[tuple[str, list[str]]]:
    merged = result["merged_roi"]
    rates = [row for row in result["merged_outcome_rates"] if row["definition"] == "mrs_2_6"]
    overlap = result["inclusion_exclusion"]
    first = [f"状态: {result['status']}; 来源运行: {result['source_run_id']}",
             f"有效总人数N={result['cohort_n']}；基底节-内囊={merged['any_overlap_n']} ({100*merged['any_overlap_proportion']:.4f}%)",
             f"基底节{overlap['basal_ganglia_n']} + 内囊{overlap['internal_capsule_n']} - 两区均累及{overlap['intersection_n']} = 并集{overlap['union_n']}",
             "分母=该合并ROI各次有效mRS0-6；分子=mRS2-6；死亡仍为mRS6。",
             "Month   Outcome_n/valid_n   Outcome%    mRS6_n"]
    for row in rates:
        value = f"{100*row['proportion']:.4f}" if row["proportion"] is not None else "NA"
        first.append(f"{row['month']:>5}   {row['numerator_n']}/{row['denominator_n']:<12} {value:>9} {row['death6_n']:>9}")
    first += ["复现检查：原有效ID集合、排除原因、总体与两组成区六次mRS频数均一致。",
              "原运行未记录每张掩膜SHA256，不能据此声称原始影像逐字节未变。"]
    second = ["按任意正重叠人数降序编号，共同分母N="+str(result["cohort_n"]),
              "Rank ROI                                Any_n      Any%"]
    for row in result["ranked_rois"]:
        second.append(f"{row['rank']:>4} {row['roi_name']:<34} {row['any_overlap_n']:>7} {100*row['any_overlap_proportion']:>9.4f}")
    second += ["同一患者可累及多个ROI；本表各ROI人数和比例不可相加。",
               "share/bgic_supplement.json和ranked_outcome_rates.csv含全部排序后的六时点数据。"]
    return [("基底节-内囊去重并集与六时点结局", first), ("合并后脑区编号与排序", second)]


def execute(run_dir: Path, out: Path | None = None) -> tuple[dict, Path]:
    run_dir = run_dir.expanduser().resolve()
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = out or run_dir.parent/("bgic_"+timestamp)
    output.mkdir(parents=True, exist_ok=False)
    share = output/"share"
    share.mkdir()
    try:
        result = summarize(run_dir, progress=lambda n,total: print(f"补统计掩膜 {n}/{total}", flush=True))
    except Exception as exc:
        private = output/"private"
        private.mkdir()
        import traceback
        (private/"failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
        result = {"status": "FAILED", "error_type": type(exc).__name__,
                  "detail": "Original run reproduction failed; see private/failure.txt. No merged counts released."}
        write_json(share/"bgic_supplement.json", result)
        write_json(output/"status.json", result)
        return result, output
    result["supplement_run_id"] = timestamp
    write_json(share/"bgic_supplement.json", result)
    for key in ("ranked_rois", "ranked_outcome_rates", "merged_distributions", "merged_recovery_distributions",
                "merged_outcome_rates", "merged_sensitivity_distributions"):
        csv_rows(share/(key+".csv"), result[key])
    chunks = []
    pages = display_pages(result)
    for index, (title, lines) in enumerate(pages, 1):
        chunk = f"===== PAGE {index}/{len(pages)} | {title} =====\n"+"\n".join(lines)+"\n"
        (share/f"page_{index:02d}.txt").write_text(chunk, encoding="utf-8")
        chunks.append(chunk)
    (share/"screenshots.txt").write_text("\n".join(chunks), encoding="utf-8")
    write_json(share/"manifest.json", {path.name: {"bytes": path.stat().st_size, "sha256": digest(path)}
                                       for path in sorted(share.iterdir())})
    with zipfile.ZipFile(output/"share.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(share.iterdir()):
            archive.write(path, arcname=path.name)
    with zipfile.ZipFile(output/"share.zip") as archive:
        if archive.testzip() is not None:
            raise RuntimeError("Supplement share ZIP integrity failed")
    write_json(output/"status.json", {"status": result["status"], "source_run_id": result["source_run_id"],
                                      "supplement_run_id": timestamp})
    return result, output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True, help="既有完整运行目录，必须保留share与private两个子目录")
    parser.add_argument("--out", type=Path, help="新的输出目录，已存在时拒绝覆盖")
    parser.add_argument("--no-pause", action="store_true")
    args = parser.parse_args()
    result, output = execute(args.run_dir, args.out)
    print("状态:", result["status"])
    print("输出:", output.resolve())
    if result["status"] != "COMPLETED":
        print("复现核对未通过，未发布合并统计。见该输出目录private/failure.txt。")
        return 2
    print("请下载:", (output/"share.zip").resolve())
    pages = display_pages(result)
    for i, (title, lines) in enumerate(pages, 1):
        print(f"\n===== PAGE {i}/{len(pages)} | {title} =====\n"+"\n".join(lines), flush=True)
        if i != len(pages) and not args.no_pause and sys.stdin.isatty():
            input("截图后按Enter继续：")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

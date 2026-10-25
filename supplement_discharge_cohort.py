#!/usr/bin/env python3
"""从原13区运行筛选出院致残性缺血性卒中队列，并重新汇总总体和四区联合人群。"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from itertools import combinations
import json
from pathlib import Path
import sys
import textwrap
import zipfile

import numpy as np

from discharge_selection import (apply_discharge_death_conflicts, load_discharge_records,
                                 select_discharge_cohort)
from mrs_summary import MONTHS, STATES, aggregate_groups
from roi_summary import analyse_masks
from run_summary import csv_rows, digest, public, write_json
from supplement_bgic import load_reproduction_context, verify_reproduction
from supplement_deep_union import UNION_ID, UNION_NAME, UNION_NAME_ZH, component_ids


def append_union_atlas(atlas: dict, components: list[str]) -> dict:
    """原13通道完整保留，追加四区二值OR；原通道编号不发生变化。"""
    if atlas.get("atlas_format") != "independent_binary_channels" or len(atlas["labels"]) != 13:
        raise ValueError("Discharge supplement requires the original independent 13-ROI atlas")
    by_id = {row["roi_id"]: row for row in atlas["labels"]}
    if len(components) != 4 or len(set(components)) != 4 or any(key not in by_id for key in components):
        raise ValueError("Four distinct component ROIs are required")
    if UNION_ID in by_id:
        raise ValueError("Reserved union ID already exists in original atlas")
    original = atlas["private_channel_masks"]
    union = np.logical_or.reduce([original[..., by_id[key]["label_id"]-1] for key in components])
    count = int(np.count_nonzero(union))
    labels = [dict(row) for row in atlas["labels"]] + [{
        "label_id": 14, "roi_id": UNION_ID, "roi_name": UNION_NAME,
        "roi_voxels": count, "roi_volume_ml": count*atlas["voxel_volume_ml"],
    }]
    return {**atlas, "private_channel_masks": np.concatenate((original, union[..., None]), axis=3),
            "labels": labels, "private_data": np.logical_or.reduce(original, axis=3)}


def group_curve(summary: dict, group: str) -> list[dict]:
    """曲线只以当次有效mRS为分母；零分母输出NA，不能显示为0%。"""
    counts = {(row["month"], row["state"]): row["n"]
              for row in summary["distributions"] if row["group"] == group}
    result = []
    for month in MONTHS:
        values = [counts[(month, str(i))] for i in range(7)]
        valid = sum(values)
        groups = [sum(values[:2]), sum(values[2:4]), sum(values[4:])]
        total = sum(counts[(month, state)] for state in STATES)
        result.append({"month": month, "total_n": total, "valid_n": valid,
            "good_recovery_n": groups[0], "functional_disability_n": groups[1],
            "functional_loss_n": groups[2],
            "good_recovery_proportion": groups[0]/valid if valid else None,
            "functional_disability_proportion": groups[1]/valid if valid else None,
            "functional_loss_proportion": groups[2]/valid if valid else None,
            "disability_death_n": groups[1]+groups[2],
            "disability_death_proportion": (groups[1]+groups[2])/valid if valid else None,
            "death6_n": values[6], "missing_n": counts[(month, "missing")],
            "invalid_n": counts[(month, "invalid")], "conflict_n": counts[(month, "conflict")]})
        if sum(groups) != valid:
            raise RuntimeError("Recovery groups do not partition valid mRS")
    return result


def overlap_audit(selected_ids: set[str], members: dict[str, set[str]], components: list[str]) -> dict:
    union_ids = set().union(*(members[key] for key in components))
    if members[UNION_ID] != union_ids:
        raise RuntimeError("Voxel union differs from the deduplicated patient set union")
    counts = Counter(identifier for key in components for identifier in members[key])
    histogram = {str(n): sum(counts.get(identifier, 0) == n for identifier in selected_ids) for n in range(5)}
    intersections, inclusion_exclusion = [], 0
    for size in range(1, 5):
        for keys in combinations(components, size):
            count = len(set.intersection(*(members[key] for key in keys)))
            inclusion_exclusion += count if size % 2 else -count
            if size > 1:
                intersections.append({"component_roi_ids": list(keys), "n": count})
    total_memberships = sum(len(members[key]) for key in components)
    if (inclusion_exclusion != len(union_ids)
            or sum(histogram.values()) != len(selected_ids)
            or sum(histogram[str(n)] for n in range(1, 5)) != len(union_ids)
            or sum(n*histogram[str(n)] for n in range(1, 5)) != total_memberships):
        raise RuntimeError("Four-region union conservation failed")
    return {"union_n": len(union_ids), "sum_component_memberships_n": total_memberships,
            "overcount_if_added_n": total_memberships-len(union_ids),
            "membership_multiplicity_histogram": histogram,
            "multiple_components_n": sum(histogram[str(n)] for n in range(2, 5)),
            "four_way_inclusion_exclusion_n": inclusion_exclusion, "intersections": intersections}


def summarize(run_dir: Path, progress=None) -> dict:
    report, states, mapping, exclusions, atlas, report_path, states_path, inventory_path = load_reproduction_context(run_dir)
    if len(report["roi"]["roi_rows"]) != 13:
        raise ValueError("Use the original completed 13-ROI run, not a union supplement")
    cfg = report["config"]
    records, source_audit = load_discharge_records(report)
    selected, selection_audit = select_discharge_cohort(states, records)
    if not selected:
        raise ValueError("No eligible discharge-disabling ischemic stroke patients in valid mask cohort")
    # SHA已验证；这里另核对原有效人群的诊断频数，避免旧队列对应错误。
    diagnosis_counts = {"infarction_1": sum(records.get(i, {}).get("D_DIAG") == "1" for i in states),
                        "TIA_2": sum(records.get(i, {}).get("D_DIAG") == "2" for i in states)}
    diagnosis_counts["other_or_missing"] = len(states)-sum(diagnosis_counts.values())
    if report.get("diagnosis_counts") is not None and diagnosis_counts != report["diagnosis_counts"]:
        raise ValueError("Original diagnosis counts no longer reproduce")
    components = component_ids(report)
    extended = append_union_atlas(atlas, components)
    roi = analyse_masks(mapping, extended, float(cfg["sensitivity_min_overlap_ml"]), progress=progress,
                        allow_unknown_units_if_exact_grid=cfg["allow_unknown_units_if_exact_grid"])
    ids13 = [row["roi_id"] for row in report["roi"]["roi_rows"]]
    verify_reproduction(report, states, exclusions, roi, ids13)
    selected_states, conflict_audit = apply_discharge_death_conflicts(
        {identifier: dict(states[identifier]) for identifier in selected}, records)
    members = {key: value & selected for key, value in roi["private_any_memberships"].items()}
    sensitive_members = {key: value & selected for key, value in roi["private_memberships"].items()}
    primary = aggregate_groups(selected_states, members)
    sensitivity = aggregate_groups(selected_states, sensitive_members)
    overlap = overlap_audit(selected, members, components)
    original_order = sorted(report["roi"]["roi_rows"], key=lambda row: (-row["any_overlap_n"], row["roi_name"]))
    n = len(selected)
    ranked = []
    for rank, original in enumerate(original_order, 1):
        key = original["roi_id"]
        # 仅保留结构标识，绝不沿用旧人群的病灶体积或覆盖统计。
        ranked.append({"rank": rank, "display_rank": rank, "roi_id": key,
            "label_id": original["label_id"], "roi_name": original["roi_name"],
            "any_overlap_n": len(members[key]), "any_overlap_proportion": len(members[key])/n,
            "n": len(sensitive_members[key]), "proportion": len(sensitive_members[key])/n,
            "denominator_n": n, "source_full_cohort_any_overlap_n": original["any_overlap_n"]})
    by_id = {row["roi_id"]: row for row in ranked}
    ranked_outcomes = [{**row, "rank": by_id[row["group"]]["rank"],
        "roi_name": by_id[row["group"]]["roi_name"],
        "lesion_n": by_id[row["group"]]["any_overlap_n"], "lesion_denominator_n": n,
        "lesion_proportion": by_id[row["group"]]["any_overlap_proportion"]}
        for row in primary["outcome_rates"] if row["group"] in by_id and row["definition"] == "mrs_2_6"]
    ranked_outcomes.sort(key=lambda row: (row["month"], row["rank"]))
    expected = {(key, month) for key in ids13 for month in MONTHS}
    if len(ranked_outcomes) != len(expected) or {(row["group"], row["month"]) for row in ranked_outcomes} != expected:
        raise RuntimeError("Incomplete 13-ROI six-visit outcome table")
    union = {"roi_id": UNION_ID, "roi_name": UNION_NAME, "roi_name_zh": UNION_NAME_ZH,
        "label_id": None, "component_roi_ids": components,
        "component_label_ids": [by_id[key]["label_id"] for key in components],
        "any_overlap_n": len(members[UNION_ID]), "any_overlap_proportion": len(members[UNION_ID])/n,
        "n": len(sensitive_members[UNION_ID]), "proportion": len(sensitive_members[UNION_ID])/n,
        "denominator_n": n, "min_overlap_ml": float(cfg["sensitivity_min_overlap_ml"])}
    conflict_n = sum(row["n"] for row in primary["distributions"] if row["group"] == "overall" and row["state"] == "conflict")
    return {"status": "COMPLETED_WITH_FLAGS" if conflict_n else "COMPLETED",
        "schema_version": "discharge_cohort_v1", "source_run_id": report["run_id"],
        "source_report_sha256": digest(report_path), "source_derived_states_sha256": digest(states_path),
        "source_inventory_sha256": digest(inventory_path), "atlas_sha256": atlas["sha256"],
        "labels_sha256": atlas["labels_sha256"], "source_full_cohort_n": len(states),
        "cohort_n": n, "months": list(MONTHS), "source_audit": source_audit,
        "selection_audit": selection_audit, "discharge_death_conflict_audit": conflict_audit,
        "original_diagnosis_counts": diagnosis_counts,
        "ranked_rois_13": ranked, "ranked_outcome_rates_13": ranked_outcomes,
        "union_roi": union, "overlap_audit": overlap,
        "components": [{"roi_id": key, "roi_name": by_id[key]["roi_name"],
                        "any_overlap_n": len(members[key])} for key in components],
        "overall_curve": group_curve(primary, "overall"), "union_curve": group_curve(primary, UNION_ID),
        "primary": primary, "sensitivity": sensitivity,
        "validation": {"all_passed": True, "valid_patient_set_reproduced": True,
            "exclusions_reproduced": True, "all_13_roi_counts_reproduced": True,
            "component_and_overall_mrs_cells_checked": 2*14*len(MONTHS)*len(STATES),
            "original_display_order_preserved": True, "patient_union_verified": True,
            "original_mask_sha256_available": False, "source_mask_sizes_reproduced": True,
            "source_diagnosis_counts_reproduced": report.get("diagnosis_counts") is not None},
        "definitions": {
            "cohort": "Original valid nonempty mask cohort AND D_DIAG=1 AND valid D_MRS in 2,3,4,5,6",
            "D_MRS_6": "Included under the literal >=2 rule and reported separately; no survival condition added",
            "discharge_death_conflict": "D_MRS6 followed by observed mRS0-5 is conflict at that visit; original files untouched; missing not filled",
            "curve_groups": "mRS0-1 good recovery, mRS2-3 functional disability, mRS4-6 functional loss; death only mRS6",
            "curve_denominator": "Valid observed mRS0-6 at each visit in the selected cohort; no imputation",
            "outcome": "mRS2-6 / valid observed mRS0-6 within each ROI at that visit",
            "roi_denominator": "All selected patients with valid nonempty masks, including those with no atlas overlap",
            "display_order": "Original full-cohort involvement ranking retained; NOT subgroup prevalence ranking",
            "union": "Exact patient union of Corona Radiata, Basal Ganglia, Internal Capsule and Centrum Semiovale",
            "sensitivity": "Threshold applied to binary union voxels; not union of separate threshold memberships",
            "early_followup": "If source mRS lacks 6 at 3/12 months, zero observed6 does not establish zero true deaths",
            "reproducibility_limit": "Original per-mask SHA256 absent; sizes, valid IDs, exclusions and all13 ROI marginals reproduced; byte identity not established"},
        "private_selection_rows": [{"patient_id": identifier, "D_DIAG": records.get(identifier, {}).get("D_DIAG", "source_absent"),
            "D_MRS": records.get(identifier, {}).get("D_MRS", "source_absent"), "selected": identifier in selected}
            for identifier in sorted(states)],
        "private_derived_states": [{"patient_id": identifier, **{f"m{month}": value for month, value in visits.items()}}
            for identifier, visits in sorted(selected_states.items())]}


def percent(value: float | None) -> str:
    return "NA" if value is None else f"{100*value:.4f}"


def curve_lines(rows: list[dict]) -> list[str]:
    lines = ["Month ValidN   0-1_n  2-3_n  4-6_n     0-1%     2-3%     4-6%   mRS6_n"]
    for row in rows:
        lines.append(f"{row['month']:>5} {row['valid_n']:>6} {row['good_recovery_n']:>7} "
            f"{row['functional_disability_n']:>6} {row['functional_loss_n']:>6} "
            + " ".join(f"{percent(row[key]):>8}" for key in ("good_recovery_proportion", "functional_disability_proportion", "functional_loss_proportion"))
            + f" {row['death6_n']:>8}")
    lines += ["分母=当次有效mRS0-6；分母>0时未舍入三类比例合计100%，零分母为NA。",
              "不填补缺失；mRS6为死亡；曲线高度必须对应实际比例。"]
    return lines


def review_pages(result: dict) -> list[tuple[str, list[str]]]:
    if not result.get("status", "").startswith("COMPLETED"):
        return [("统计未完成", ["状态: "+str(result.get("status")),
            "原因: "+str(result.get("safe_error", "复现或统计未通过")),
            "请在工作站查看输出目录private/failure.txt；未发布患者统计。"])]
    flow = result["selection_audit"]
    excluded = flow["mutually_exclusive_exclusions"]
    source = result["source_audit"]
    conflicts = result["discharge_death_conflict_audit"]
    pages = [("筛选定义及来源核对", [
        f"状态: {result['status']}; 原运行: {result['source_run_id']}",
        "原有效非空掩膜患者 ∩ D_DIAG=1 ∩ 有效D_MRS=2–6。",
        f"原有效总N={result['source_full_cohort_n']}; 本次入组N={result['cohort_n']}",
        "D_MRS=6按>=2原定义纳入并单列；不额外添加出院存活限制。",
        f"原有效人群中：早期记录={flow['early_record_n']}；D_DIAG=1人数={flow['ischemic_n']}。",
        "排除（互斥）：无早期记录={no_early_record}；TIA={diagnosis_tia_2}；诊断缺失={diagnosis_missing}；诊断非法={diagnosis_invalid}；其他诊断={diagnosis_other}".format(**excluded),
        "缺血性卒中中排除：D_MRS缺失={ischemic_discharge_mrs_missing}；非法={ischemic_discharge_mrs_invalid}；0–1分={ischemic_discharge_mrs_0_1}".format(**excluded),
        "入组出院mRS 2/3/4/5/6人数："+" / ".join(str(flow["included_discharge_mrs_counts"][str(k)]) for k in range(2, 7)),
        f"出院6后新增冲突：{conflicts['new_conflict_patient_n']}人 / {conflicts['new_conflict_visit_n']}访视；每访视="
            +" / ".join(str(conflicts["new_conflict_visit_counts"][str(m)]) for m in MONTHS),
        f"SAS记录数={source.get('rows')}; ID存储={source.get('id_storage')}; 原SAS的SHA256一致={source.get('original_sha256_matched')}",
        *[f"字段 {key} <- {info['actual']}; label={info.get('label')}; type={info.get('type')}"
          for key, info in source.get("schema", {}).items()],
        "原总体与13区两阈值六访视mRS 1680格全部复现。",
        "原始掩膜无逐张SHA256，仅复现大小、有效集合、排除与边际统计。"])]
    lines = ["13区沿用原全体人群编号和顺序；下列比例为本次筛选队列，不按新比例重排。",
             f"共同分母N={result['cohort_n']}; 主定义any>0; 敏感性体积阈值见配置。",
             "Rank ROI                        Any_n      Any%   Sens_n"]
    for row in result["ranked_rois_13"]:
        lines.append(f"{row['rank']:>4} {row['roi_name']:<25} {row['any_overlap_n']:>6} {percent(row['any_overlap_proportion']):>9} {row['n']:>8}")
    pages.append(("13个ROI累及人数及比例", lines))
    pages.append(("筛选队列总体功能分组曲线", curve_lines(result["overall_curve"])))
    union = result["union_roi"]
    pages.append(("四区联合去重及曲线", [
        "放射冠＋基底节＋内囊＋半卵圆中心，患者并集，每人一次。",
        f"联合N={union['any_overlap_n']} ({percent(union['any_overlap_proportion'])}%); 筛选队列N={result['cohort_n']}",
        "各组成区: "+", ".join(f"{r['roi_name']}={r['any_overlap_n']}" for r in result["components"]),
        "累及0/1/2/3/4区人数: "+" / ".join(str(result["overlap_audit"]["membership_multiplicity_histogram"][str(n)]) for n in range(5)),
        *curve_lines(result["union_curve"])]))
    for month in MONTHS:
        rows = [row for row in result["ranked_outcome_rates_13"] if row["month"] == month]
        lines = ["残疾/死亡率=mRS2–6 / 该ROI当次有效mRS0–6；mRS6为死亡。",
                 "Rank ROI                      Lesion_n  Lesion%   Outcome_n/valid_n  Outcome%"]
        for row in rows:
            lines.append(f"{row['rank']:>4} {row['roi_name']:<24} {row['lesion_n']:>6} {percent(row['lesion_proportion']):>9} "
                f"{row['numerator_n']:>6}/{row['denominator_n']:<6} {percent(row['proportion']):>10}")
        lines += ["累及比例分母为筛选队列总N；多个ROI患者可以重叠，不能相加。",
                  "3/12月若源字段无6分，不能将零计数解释为没有真实死亡。"]
        pages.append((f"{month}个月13区残疾/死亡率", lines))
    # 长标签和小终端可分成额外页面，核心13行ROI表仍完整保留。
    sized_pages = []
    for title, lines in pages:
        wrapped = [part for line in lines for part in (textwrap.wrap(line, 110) or [""])]
        sized_pages.extend((title, wrapped[start:start+24]) for start in range(0, len(wrapped), 24))
    return sized_pages


def save_review(result: dict, share: Path) -> None:
    pages = review_pages(result)
    contents = []
    for index, (title, lines) in enumerate(pages, 1):
        value = f"===== PAGE {index}/{len(pages)} | {title} =====\n"+"\n".join(lines)+"\n"
        (share/f"page_{index:02d}.txt").write_text(value, encoding="utf8")
        contents.append(value)
    (share/"screenshots.txt").write_text("\n".join(contents), encoding="utf8")


def execute(run_dir: Path, out: Path | None = None) -> tuple[dict, Path]:
    run_dir = run_dir.expanduser().resolve()
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = out or run_dir.parent/("discharge_"+timestamp)
    output.mkdir(parents=True, exist_ok=False)
    share, private = output/"share", output/"private"
    share.mkdir(); private.mkdir()
    try:
        result = summarize(run_dir, progress=lambda n,total:print(f"出院队列补统计 {n}/{total}", flush=True))
    except Exception as exc:
        import traceback
        (private/"failure.txt").write_text(traceback.format_exc(), encoding="utf8")
        # 只公开已知固定错误文本，外部I/O异常可能含患者路径，保留在private。
        safe_errors = {
            "Original early SAS is unavailable": "原早期SAS文件不可用，请核对原配置路径。",
            "Original early SAS SHA256 fingerprint is missing or ambiguous": "原报告缺少唯一的早期SAS指纹。",
            "Early SAS differs from the original run fingerprint": "早期SAS与原运行的SHA256不一致，未混用新旧数据。",
            "Required fields code_n, D_DIAG and D_MRS must each exist exactly once": "code_n、D_DIAG或D_MRS字段缺失或重复；未猜测替代字段。",
            "Early SAS changed while discharge fields were being read": "读取期间早期SAS发生变化。",
            "Early SAS has missing, invalid or duplicate normalized patient IDs": "早期SAS存在缺失、非法或重复ID。",
            "Early SAS row count differs from the original audit": "早期SAS记录数与原审计不一致。",
            "Early SAS ID storage differs from the original audit": "早期SAS的ID存储类型与原审计不一致。",
            "Original diagnosis counts no longer reproduce": "原有效人群的诊断人数未能复现。",
            "No eligible discharge-disabling ischemic stroke patients in valid mask cohort": "原有效掩膜人群中没有满足D_DIAG=1且D_MRS=2–6的患者。",
        }
        result = {"status": "FAILED", "error_type": type(exc).__name__,
                  "safe_error": safe_errors.get(str(exc), "原运行复现或统计未通过；具体原因保留在工作站private/failure.txt。"),
                  "detail": "Source/selection/reproduction failed; inspect private/failure.txt. No cohort statistics released."}
        write_json(share/"discharge_cohort.json", result)
        write_json(output/"status.json", result)
        return result, output
    csv_rows(private/"selection.csv", result.pop("private_selection_rows"))
    csv_rows(private/"derived_states.csv", result.pop("private_derived_states"))
    result["supplement_run_id"] = timestamp
    result = public(result)
    write_json(share/"discharge_cohort.json", result)
    for key in ("overall_curve", "union_curve", "ranked_rois_13", "ranked_outcome_rates_13"):
        csv_rows(share/(key+".csv"), result[key])
    for version in ("primary", "sensitivity"):
        for key in ("distributions", "recovery_distributions", "outcome_rates", "transitions", "change"):
            csv_rows(share/f"{version}_{key}.csv", result[version][key])
    save_review(result, share)
    write_json(share/"manifest.json", {p.name: {"bytes": p.stat().st_size, "sha256": digest(p)} for p in sorted(share.iterdir())})
    with zipfile.ZipFile(output/"share.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(share.iterdir()):
            archive.write(path, arcname=path.name)
    with zipfile.ZipFile(output/"share.zip") as archive:
        if archive.testzip() is not None:
            raise RuntimeError("Discharge cohort ZIP integrity failed")
    write_json(output/"status.json", {"status": result["status"], "source_run_id": result["source_run_id"], "supplement_run_id": timestamp})
    return result, output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--run-dir", type=Path, help="原13区运行目录（需保留private），不是四区补算目录")
    source.add_argument("--review", type=Path, help="只重放已有share/discharge_cohort.json")
    parser.add_argument("--out", type=Path, help="新目录，拒绝覆盖")
    parser.add_argument("--page", type=int)
    parser.add_argument("--no-pause", action="store_true")
    args = parser.parse_args()
    if args.review:
        if args.out:
            parser.error("--out cannot be combined with --review")
        result = json.loads(args.review.read_text(encoding="utf8"))
    else:
        result, output = execute(args.run_dir, args.out)
        print("状态:", result["status"]); print("输出:", output.resolve())
        if result["status"].startswith("COMPLETED"):
            print("请下载:", (output/"share.zip").resolve())
    pages = review_pages(result)
    if args.page is not None and not 1 <= args.page <= len(pages):
        parser.error(f"--page must be 1..{len(pages)}")
    chosen = [(args.page-1, pages[args.page-1])] if args.page else list(enumerate(pages))
    for position, (index, (title, lines)) in enumerate(chosen):
        print(f"\n===== PAGE {index+1}/{len(pages)} | {title} =====\n"+"\n".join(lines), flush=True)
        if not args.no_pause and sys.stdin.isatty() and position < len(chosen)-1:
            input("截图后按Enter继续：")
    return 0 if result["status"].startswith("COMPLETED") else 2


if __name__ == "__main__":
    raise SystemExit(main())

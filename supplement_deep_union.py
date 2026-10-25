#!/usr/bin/env python3
"""补算放射冠、基底节、内囊、半卵圆中心的去重联合人群，保留13区独立显示。"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
import sys
import zipfile

from mrs_summary import MONTHS, STATES
from roi_summary import analyse_masks
from run_summary import csv_rows, digest, write_json
from supplement_bgic import load_reproduction_context, union_atlas, verify_reproduction

COMPONENT_NAMES = ("Corona Radiata", "Basal Ganglia", "Internal Capsule", "Centrum Semiovale")
UNION_ID = "corona_basal_capsule_semiovale_union"
UNION_NAME = "Corona Radiata-Basal Ganglia-Internal Capsule-Centrum Semiovale union"
UNION_NAME_ZH = "放射冠、基底节、内囊、半卵圆中心联合人群"


def component_ids(report: dict) -> list[str]:
    identifiers = []
    for name in COMPONENT_NAMES:
        found = [row["roi_id"] for row in report["roi"]["roi_rows"] if row["roi_name"] == name]
        if len(found) != 1:
            raise ValueError("Source run must contain all four distinct union component ROIs")
        identifiers.append(found[0])
    if len(set(identifiers)) != 4:
        raise ValueError("Union component ROI IDs must be distinct")
    return identifiers


def summarize(run_dir: Path, progress=None) -> dict:
    report, states, mapping, old_exclusions, atlas, report_path, states_path, inventory_path = load_reproduction_context(run_dir)
    cfg = report["config"]
    cohort_n = int(report["roi"]["valid_n"])
    if len(report["roi"]["roi_rows"]) != 13:
        raise ValueError("This supplement expects the original independent 13-ROI run")
    components = component_ids(report)
    selected_atlas = union_atlas(atlas, components, UNION_ID, UNION_NAME)
    roi = analyse_masks(mapping, selected_atlas, float(cfg["sensitivity_min_overlap_ml"]), progress=progress,
                        allow_unknown_units_if_exact_grid=cfg["allow_unknown_units_if_exact_grid"])
    primary, sensitivity, original_rows, current_rows = verify_reproduction(report, states, old_exclusions, roi, components)
    members = roi["private_any_memberships"]
    union_ids = set().union(*(members[key] for key in components))
    if union_ids != members[UNION_ID]:
        raise RuntimeError("Voxel union membership differs from the exact patient set union")
    occurrence = Counter(identifier for key in components for identifier in members[key])
    histogram = {str(n): sum(occurrence.get(identifier, 0) == n for identifier in states) for n in range(5)}
    sum_memberships = sum(len(members[key]) for key in components)
    if (sum(histogram[str(n)] for n in range(1,5)) != len(union_ids)
            or sum(n*histogram[str(n)] for n in range(1,5)) != sum_memberships
            or sum(histogram.values()) != cohort_n):
        raise RuntimeError("Patient union histogram conservation failed")
    intersections = []
    inclusion_exclusion_n = 0
    for size in range(1,5):
        for keys in combinations(components, size):
            count = len(set.intersection(*(members[key] for key in keys)))
            inclusion_exclusion_n += count if size % 2 else -count
            if size > 1:
                intersections.append({"component_roi_ids": list(keys), "n": count})
    if inclusion_exclusion_n != len(union_ids):
        raise RuntimeError("Four-component inclusion-exclusion failed")
    union_row = {**current_rows[UNION_ID], "label_id": None, "roi_name_zh": UNION_NAME_ZH,
                 "component_roi_ids": components,
                 "component_label_ids": [original_rows[key]["label_id"] for key in components],
                 "min_overlap_ml": float(cfg["sensitivity_min_overlap_ml"])}
    union_distributions = [row for row in primary["distributions"] if row["group"] == UNION_ID]
    union_recovery = [row for row in primary["recovery_distributions"] if row["group"] == UNION_ID]
    union_outcomes = [row for row in primary["outcome_rates"] if row["group"] == UNION_ID]
    distribution_lookup = {(row["month"], row["state"]): row["n"] for row in union_distributions}
    curve = []
    for month in MONTHS:
        counts = [distribution_lookup[(month, str(i))] for i in range(7)]
        n = sum(counts)
        groups = [sum(counts[:2]), sum(counts[2:4]), sum(counts[4:])]
        if sum(groups) != n:
            raise RuntimeError("Union recovery groups do not partition the valid visit denominator")
        curve.append({"month":month, "total_n":len(union_ids), "valid_n":n,
                      "good_recovery_n":groups[0], "functional_disability_n":groups[1], "functional_loss_n":groups[2],
                      "good_recovery_proportion":groups[0]/n if n else None,
                      "functional_disability_proportion":groups[1]/n if n else None,
                      "functional_loss_proportion":groups[2]/n if n else None,
                      "disability_death_n":groups[1]+groups[2],
                      "disability_death_proportion":(groups[1]+groups[2])/n if n else None,
                      "death6_n":counts[6]})
    # 联合人群单独出图，原13区既不删除也不替换。
    ranked = [{**row} for row in report["roi"]["roi_rows"]]
    ranked.sort(key=lambda row:(-row["any_overlap_n"],row["roi_name"]))
    for rank,row in enumerate(ranked,1):
        row["rank"] = rank
    by_id = {row["roi_id"]:row for row in ranked}
    ranked_outcomes = [{**row, "rank":by_id[row["group"]]["rank"], "roi_name":by_id[row["group"]]["roi_name"],
                       "lesion_n":by_id[row["group"]]["any_overlap_n"], "lesion_denominator_n":cohort_n,
                       "lesion_proportion":by_id[row["group"]]["any_overlap_n"]/cohort_n}
                      for row in report["primary"]["outcome_rates"]
                      if row["group"] != "overall" and row["definition"] == "mrs_2_6"]
    ranked_outcomes.sort(key=lambda row:(row["month"],row["rank"]))
    expected = {(key,month) for key in by_id for month in MONTHS}
    if len(ranked_outcomes) != len(expected) or {(row["group"],row["month"]) for row in ranked_outcomes} != expected:
        raise ValueError("Original independent ROI outcome table is incomplete")
    return {"status":"COMPLETED", "schema_version":"deep_union_v1", "source_run_id":report["run_id"],
            "source_report_sha256":digest(report_path), "source_derived_states_sha256":digest(states_path),
            "source_inventory_sha256":digest(inventory_path), "atlas_sha256":atlas["sha256"],
            "labels_sha256":atlas["labels_sha256"], "cohort_n":cohort_n, "months":list(MONTHS),
            "union_roi":union_row,
            "components":[{"roi_id":key,"label_id":original_rows[key]["label_id"],
                           "roi_name":original_rows[key]["roi_name"],"any_overlap_n":len(members[key])}
                          for key in components],
            "overlap_audit":{"union_n":len(union_ids), "sum_component_memberships_n":sum_memberships,
                "overcount_if_added_n":sum_memberships-len(union_ids),
                "membership_multiplicity_histogram":histogram,
                "multiple_components_n":sum(histogram[str(n)] for n in range(2,5)),
                "intersections":intersections, "four_way_inclusion_exclusion_n":inclusion_exclusion_n},
            "union_distributions":union_distributions, "union_recovery_distributions":union_recovery,
            "union_outcome_rates":union_outcomes, "union_curve":curve,
            "union_sensitivity_distributions":[row for row in sensitivity["distributions"] if row["group"]==UNION_ID],
            "union_sensitivity_recovery_distributions":[row for row in sensitivity["recovery_distributions"] if row["group"]==UNION_ID],
            "ranked_rois_13":ranked, "ranked_outcome_rates_13":ranked_outcomes,
            "validation":{"all_passed":True,"valid_patient_set_reproduced":True,"exclusions_reproduced":True,
                "component_counts_reproduced":True,"component_and_overall_mrs_cells_checked":2*5*len(MONTHS)*len(STATES),
                "patient_union_verified":True,"source_mask_sizes_reproduced":True,"original_mask_sha256_available":False,
                "all_original_roi_rows_preserved":True,"original_roi_count":len(ranked)},
            "definitions":{"membership":"Any positive overlap with any of the four component ROIs; each patient counted once",
                "curve_groups":"mRS0-1 good recovery, mRS2-3 functional disability, mRS4-6 functional loss; mRS6 remains death",
                "curve_denominator":"Valid observed mRS0-6 within the union cohort at each visit; no imputation or carry-forward",
                "roi_display":"Original independent ROIs remain unchanged and sorted by any-overlap counts; the union is a separate cohort, not a replacement ROI",
                "sensitivity":"Volume threshold applied to binary union voxels, counting shared voxels only once",
                "reproducibility_limit":"Original per-mask SHA256 fingerprints were not recorded. File sizes, exact valid IDs, exclusions and all four component marginals reproduce; byte-level image immutability is not established."}}


def display_page(result: dict) -> str:
    union = result["union_roi"]
    lines = ["===== PAGE 1/1 | 四脑区去重联合人群及六次结局 =====",
        "组成区：放射冠 + 基底节 + 内囊 + 半卵圆中心（患者并集，不相加）",
        f"有效总N={result['cohort_n']}；联合人群N={union['any_overlap_n']} ({100*union['any_overlap_proportion']:.4f}%)",
        "原13脑区保留独立统计；本联合人群仅另加一张曲线。",
        "组成区人数："+", ".join(f"{r['roi_name']}={r['any_overlap_n']}" for r in result["components"]),
        "累及1/2/3/4个组成区的人数："+" / ".join(str(result["overlap_audit"]["membership_multiplicity_histogram"][str(k)]) for k in range(1,5)),
        "Month ValidN  0-1_n   2-3_n   4-6_n    0-1%    2-3%    4-6%   2-6_n    2-6%"]
    for r in result["union_curve"]:
        proportions = [r[key] for key in ("good_recovery_proportion","functional_disability_proportion","functional_loss_proportion","disability_death_proportion")]
        values = [f"{100*x:.4f}" if x is not None else "NA" for x in proportions]
        lines.append(f"{r['month']:>5} {r['valid_n']:>6} {r['good_recovery_n']:>6} {r['functional_disability_n']:>7} {r['functional_loss_n']:>7} "
                     +" ".join(f"{v:>7}" for v in values[:3])+f" {r['disability_death_n']:>7} {values[3]:>8}")
    lines += ["分母为该次有效mRS0-6；有效n>0时三组未舍入比例合计100%。死亡仍为mRS6。",
        "复现核对：原有效ID集合、排除原因、总体和四组成区全部六次mRS频数一致。",
        "原运行无逐张掩膜SHA256，不能据此声称影像逐字节未变。",
        "share/deep_union_supplement.json含精确计数、原13区排序表与独立联合曲线数据。"]
    return "\n".join(lines)+"\n"


def execute(run_dir: Path, out: Path | None = None) -> tuple[dict,Path]:
    run_dir = run_dir.expanduser().resolve()
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = out or run_dir.parent/("deep_union_"+timestamp)
    output.mkdir(parents=True,exist_ok=False)
    share = output/"share";share.mkdir()
    try:
        result = summarize(run_dir,progress=lambda n,total:print(f"联合人群补统计 {n}/{total}",flush=True))
    except Exception as exc:
        import traceback
        private=output/"private";private.mkdir()
        (private/"failure.txt").write_text(traceback.format_exc(),encoding="utf8")
        result={"status":"FAILED","error_type":type(exc).__name__,
                "detail":"Original run reproduction failed; see private/failure.txt. No union counts released."}
        write_json(share/"deep_union_supplement.json",result)
        write_json(output/"status.json",result)
        return result,output
    result["supplement_run_id"] = timestamp
    write_json(share/"deep_union_supplement.json",result)
    for table in ("union_distributions","union_recovery_distributions","union_outcome_rates","union_curve",
                  "union_sensitivity_distributions","union_sensitivity_recovery_distributions", "ranked_rois_13","ranked_outcome_rates_13"):
        csv_rows(share/(table+".csv"),result[table])
    page = display_page(result)
    (share/"page_01.txt").write_text(page,encoding="utf8")
    (share/"screenshots.txt").write_text(page,encoding="utf8")
    write_json(share/"manifest.json",{p.name:{"bytes":p.stat().st_size,"sha256":digest(p)} for p in sorted(share.iterdir())})
    with zipfile.ZipFile(output/"share.zip","x",zipfile.ZIP_DEFLATED) as archive:
        for p in sorted(share.iterdir()):archive.write(p,arcname=p.name)
    with zipfile.ZipFile(output/"share.zip") as archive:
        if archive.testzip() is not None:raise RuntimeError("Union share ZIP integrity failed")
    write_json(output/"status.json",{"status":result["status"],"source_run_id":result["source_run_id"],"supplement_run_id":timestamp})
    return result,output


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir",type=Path,required=True)
    parser.add_argument("--out",type=Path,help="新的输出目录，拒绝覆盖已有目录")
    args=parser.parse_args()
    result,output=execute(args.run_dir,args.out)
    print("状态:",result["status"]);print("输出:",output.resolve())
    if result["status"]!="COMPLETED":
        print("复现未通过，未发布联合人群统计。请查看该输出目录private/failure.txt。")
        return 2
    print("请下载:",(output/"share.zip").resolve())
    print(display_page(result),flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""仅从SAS汇总全部缺血性卒中与出院致残亚组；不读取或核查任何影像。"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import textwrap
import zipfile

from discharge_selection import apply_discharge_death_conflicts
from inputs import load_clinical
from mrs_summary import MONTHS, MRS_FIELDS, RAW_STATES, STATES, aggregate_groups, derive_states, parse_mrs

DISABLED = "discharge_disabling"


class ClinicalInputsError(ValueError):
    def __init__(self, audit: dict):
        super().__init__("Clinical fields or IDs failed validation")
        self.audit = audit


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            value.update(block)
    return value.hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf8")


def csv_rows(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        if not rows:
            return
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows(rows)


def anonymous_audit(value):
    """临床读取审计只保留字段、类型、计数；源路径留在工作站。"""
    if isinstance(value, dict):
        return {key: anonymous_audit(item) for key, item in value.items()
                if key not in {"path", "source_files"} and not str(key).startswith("private_")}
    if isinstance(value, list):
        return [anonymous_audit(item) for item in value]
    return value


def curve(summary: dict, group: str) -> list[dict]:
    cells = {(r["month"], r["state"]): r["n"] for r in summary["distributions"] if r["group"] == group}
    rows = []
    for month in MONTHS:
        counts = [cells[month, str(i)] for i in range(7)]
        valid = sum(counts)
        grouped = (sum(counts[:2]), sum(counts[2:4]), sum(counts[4:]))
        assert sum(grouped) == valid
        rows.append({"month": month, "valid_n": valid,
            "total_n": sum(cells[month, state] for state in STATES),
            "good_n": grouped[0], "disability_n": grouped[1], "loss_n": grouped[2],
            "good_percent": grouped[0]*100/valid if valid else None,
            "disability_percent": grouped[1]*100/valid if valid else None,
            "loss_percent": grouped[2]*100/valid if valid else None,
            "death6_n": counts[6], "missing_n": cells[month, "missing"],
            "invalid_n": cells[month, "invalid"], "conflict_n": cells[month, "conflict"]})
    return rows


def summarize(config: dict, base: Path) -> dict:
    # 只复制临床字段配置；lesion_root、atlas与原影像报告完全不参与本入口。
    source_cfg = {key: config.get(key) for key in ("sas_encoding", "id_width")}
    configured = config.get("source_files", {})
    if not isinstance(configured, dict) or any(not configured.get(role) for role in ("early", "late")):
        raise ValueError("Both early and late SAS paths are required")
    paths = {}
    for role in ("early", "late"):
        path = Path(configured[role]).expanduser()
        paths[role] = (path if path.is_absolute() else base/path).resolve()
    source_cfg["source_files"] = {role: str(path) for role, path in paths.items()}
    fingerprints = {role: {"bytes": path.stat().st_size, "sha256": digest(path)} for role, path in paths.items()}
    records, audit = load_clinical(source_cfg, include_discharge_mrs=True)
    if not audit["clinical_allowed"]:
        raise ClinicalInputsError(anonymous_audit(audit))
    if len({info.get("id_storage") for info in audit["sources"].values()}) != 1:
        raise ValueError("SAS ID storage types differ between sources")
    for role, path in paths.items():
        if fingerprints[role] != {"bytes": path.stat().st_size, "sha256": digest(path)}:
            raise ValueError("SAS changed while clinical records were being read")
    canonical = {identifier: {"D_DIAG": parse_mrs(row.get("D_DIAG")), "D_MRS": parse_mrs(row.get("D_MRS"))}
                 for identifier, row in records.items()}
    diagnoses = Counter({key: 0 for key in ("ischemic_1", "TIA_2", "diagnosis_missing", "diagnosis_invalid", "diagnosis_other", "no_early_record")})
    for identifier, row in records.items():
        diagnosis = canonical[identifier]["D_DIAG"]
        if "D_DIAG" not in row:
            category = "no_early_record"
        else:
            category = {"1": "ischemic_1", "2": "TIA_2", "missing": "diagnosis_missing", "invalid": "diagnosis_invalid"}.get(diagnosis, "diagnosis_other")
        diagnoses[category] += 1
    ischemic = {identifier: row for identifier, row in records.items() if canonical[identifier]["D_DIAG"] == "1"}
    if not ischemic:
        raise ValueError("No patients meet D_DIAG=1 in the SAS records")
    states, mrs_audit = derive_states(ischemic)
    disabled_ids = {identifier for identifier in ischemic if canonical[identifier]["D_MRS"] in ("2", "3", "4", "5", "6")}
    discharge6 = {identifier: visits for identifier, visits in states.items() if canonical[identifier]["D_MRS"] == "6"}
    corrected, discharge_conflicts = apply_discharge_death_conflicts(discharge6, canonical)
    states.update(corrected)
    summary = aggregate_groups(states, {DISABLED: disabled_ids})
    distributions = {"overall": {state: 0 for state in RAW_STATES}, "disabled": {state: 0 for state in RAW_STATES}}
    for identifier in ischemic:
        value = canonical[identifier]["D_MRS"]
        distributions["overall"][value] += 1
        if identifier in disabled_ids:
            distributions["disabled"][value] += 1
    overall_n, disabled_n = len(ischemic), len(disabled_ids)
    exclusions = {"discharge_mrs_0_1": sum(distributions["overall"][str(i)] for i in (0, 1)),
                  "discharge_mrs_missing": distributions["overall"]["missing"],
                  "discharge_mrs_invalid": distributions["overall"]["invalid"]}
    if (sum(diagnoses.values()) != len(records) or diagnoses["ischemic_1"] != overall_n
            or sum(exclusions.values())+disabled_n != overall_n):
        raise RuntimeError("Clinical cohort flow conservation failed")
    flags = []
    if not disabled_n:
        flags.append("No discharge-disabling patients; subgroup curve has no valid percentages")
    if any(row["state"] == "conflict" and row["n"] for row in summary["distributions"]):
        flags.append("Alive mRS after previous mRS6 is flagged as conflict at that visit")
    if any(row["valid_n"] == 0 for row in curve(summary, "overall")+curve(summary, DISABLED)):
        flags.append("Some visit denominators are zero; percentages are null, not zero")
    return {"schema_version": "clinical_trajectories_v1", "status": "COMPLETED_WITH_FLAGS" if flags else "COMPLETED",
        "issues": flags, "months": list(MONTHS), "overall_n": overall_n, "disabled_n": disabled_n,
        "clinical_audit": anonymous_audit(audit), "source_fingerprints": fingerprints,
        "selection": {"clinical_union_n": len(records), "diagnosis_counts": dict(diagnoses),
            "overall_n": overall_n, "disabled_n": disabled_n, "disabled_exclusions": exclusions,
            "flow_conservation_passed": True},
        "discharge_mrs_counts": distributions,
        "followup_only_mrs_audit_before_discharge6_conflicts": mrs_audit,
        "discharge_death_conflict_audit": discharge_conflicts,
        "overall_curve": curve(summary, "overall"), "disabled_curve": curve(summary, DISABLED),
        "primary_distributions": [row for row in summary["distributions"] if row["group"] == "overall"],
        "disabled_distributions": [row for row in summary["distributions"] if row["group"] == DISABLED],
        "recovery_distributions": summary["recovery_distributions"],
        "outcome_rates": summary["outcome_rates"], "transitions": summary["transitions"],
        "validation": {**{k:v for k,v in summary["validation"].items()
            if k not in {"roi_membership_assignments_n", "patients_in_multiple_roi_groups_n"}},
            "disabled_subgroup_membership_n": disabled_n,
            "no_imaging_filter": True, "source_unchanged_during_read": True},
        "definitions": {
            "overall": "All unique SAS patients with D_DIAG=1; discharge mRS is not an inclusion condition",
            "disabled": "Subset of overall with valid integer D_MRS in2..6; includes6 under the literal>=2 rule",
            "imaging": "No mask/atlas existence, validity, nonempty lesion or matching requirement; no imaging input read",
            "brain_map": "Any separately displayed lesion map has its own imaging cohort and denominator; not this SAS cohort",
            "groups": "good=mRS0..1, disability=mRS2..3, loss=mRS4..6; mRS6 means death",
            "percent_units": "Percent0..100; exact grouped counts / valid observed mRS0..6 at that visit *100",
            "missing": "No complete-case restriction or imputation; missing/invalid/conflict excluded only from visit denominators",
            "discharge6": "Observed0..5 after discharge6 becomes conflict at that visit; missing remains missing",
            "early_death_limit": "Source fields without6 at3/12 months cannot establish zero true mortality",
        },
        "private_source_files": source_cfg["source_files"],
        "private_selection_rows": [{"patient_id": identifier, **canonical[identifier],
            "overall_selected": identifier in ischemic, "disabled_selected": identifier in disabled_ids,
            **{field: records[identifier].get(field) for field in MRS_FIELDS.values()}}
            for identifier in sorted(records)],
        "private_derived_states": [{"patient_id": identifier, "disabled_selected": identifier in disabled_ids,
            **{f"m{month}": value for month, value in visits.items()}}
            for identifier, visits in sorted(states.items())]}


def review_pages(result: dict) -> list[tuple[str, list[str]]]:
    if not result.get("status", "").startswith("COMPLETED"):
        lines = ["状态: "+str(result.get("status")), "原因: "+result.get("safe_error", "临床输入或统计未通过。")]
        if result.get("clinical_audit"):
            lines.append(json.dumps(result["clinical_audit"], ensure_ascii=False))
        lines.append("详细错误留在工作站private/failure.txt。未发布轨迹统计。")
        return [("尚未完成", lines)]
    audit = result["clinical_audit"]
    flow = result["selection"]
    counts = result["discharge_mrs_counts"]
    lines = [f"状态: {result['status']}; 本次运行: {result.get('run_id', 'not_saved')}",
        "仅SAS记录：全人群D_DIAG=1；出院致残亚组另要求有效D_MRS=2–6。",
        "完全不读取或核查病灶/图谱，不以是否有影像或病灶非空筛选。",
        f"两SAS并集N={flow['clinical_union_n']}; 缺血性卒中N={result['overall_n']}; 出院致残N={result['disabled_n']}",
        "诊断分布: "+json.dumps(flow["diagnosis_counts"], ensure_ascii=False),
        "致残亚组排除: "+json.dumps(flow["disabled_exclusions"], ensure_ascii=False),
        "缺血性卒中出院mRS0/1/2/3/4/5/6/缺失/非法: "+" / ".join(str(counts['overall'][s]) for s in RAW_STATES),
        "致残亚组出院mRS2/3/4/5/6: "+" / ".join(str(counts['disabled'][str(s)]) for s in range(2,7)),
        "连接范围: "+", ".join(f"{key}={audit.get(key)}" for key in ("early_n", "late_n", "early_only_n", "late_only_n", "joined_n")),
        "SAS读取前后SHA256一致；仅读取ID、诊断、出院及六次随访mRS。",
        "出院6后新增冲突访视: "+str(result['discharge_death_conflict_audit']['new_conflict_visit_counts']),
        "图中若另配病灶累及脑图，其影像队列分母独立标注，不与本次SAS队列混用。"]
    for source, details in audit.get("sources", {}).items():
        for field, schema in details.get("schema", {}).items():
            if field in ("D_DIAG", "D_MRS"):
                lines.append(f"{source} {field} <- {schema['actual']}; label={schema.get('label')}")
    pages = [("SAS队列及筛选核对", lines)]
    for title, key, n in (("全部缺血性卒中恢复分组", "overall_curve", result["overall_n"]),
                          ("出院致残性缺血性卒中恢复分组", "disabled_curve", result["disabled_n"])):
        lines = [f"固定队列N={n}；比例以该次有效mRS0–6为分母。",
            "Month ValidN  0-1_n  2-3_n  4-6_n     0-1%     2-3%     4-6%  mRS6_n Miss Invalid Conflict"]
        for row in result[key]:
            percentages = ["NA" if row[k] is None else f"{row[k]:.4f}" for k in ("good_percent", "disability_percent", "loss_percent")]
            lines.append(f"{row['month']:>5} {row['valid_n']:>6} {row['good_n']:>6} {row['disability_n']:>6} {row['loss_n']:>6} "
                +" ".join(f"{p:>8}" for p in percentages)+f" {row['death6_n']:>7} {row['missing_n']:>4} {row['invalid_n']:>7} {row['conflict_n']:>8}")
        lines += ["三组为mRS0–1、2–3、4–6；有效n>0时未舍入比例合计100%，零分母为NA。",
            "PPT不绘制缺失组，不填补结局。死亡仅由mRS6定义。",
            "3/12月源字段无6分时，不代表真实死亡为零。百分比曲线不代表个体轨迹。"]
        pages.append((title, lines))
    output = []
    for title, lines in pages:
        wrapped = [part for line in lines for part in (textwrap.wrap(line, 120) or [""])]
        output.extend((title, wrapped[start:start+24]) for start in range(0,len(wrapped),24))
    return output


def save_review(result: dict, directory: Path) -> None:
    pages = review_pages(result)
    texts = []
    for index,(title,lines) in enumerate(pages,1):
        value=f"===== PAGE {index}/{len(pages)} | {title} =====\n"+"\n".join(lines)+"\n"
        (directory/f"page_{index:02d}.txt").write_text(value,encoding="utf8")
        texts.append(value)
    (directory/"screenshots.txt").write_text("\n".join(texts),encoding="utf8")


def execute(config_path: Path, out: Path | None = None) -> tuple[dict, Path]:
    config_path=config_path.expanduser().resolve()
    config=json.loads(config_path.read_text(encoding="utf8"))
    timestamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output=out or config_path.parent/"outputs"/("clinical_trajectories_"+timestamp)
    output.mkdir(parents=True,exist_ok=False)
    share,private=output/"share",output/"private"
    share.mkdir();private.mkdir()
    try:
        result=summarize(config,config_path.parent)
    except Exception as exc:
        import traceback
        (private/"failure.txt").write_text(traceback.format_exc(),encoding="utf8")
        messages={"Both early and late SAS paths are required":"配置需包含early和late两个SAS路径。",
                  "SAS ID storage types differ between sources":"两SAS的ID存储类型不一致，未猜测跨来源ID转换。",
                  "SAS changed while clinical records were being read":"读取期间SAS文件发生变化，未发布结果。",
                  "No patients meet D_DIAG=1 in the SAS records":"SAS中没有满足D_DIAG=1的患者。"}
        result={"status":"INPUTS_REQUIRED" if isinstance(exc,ClinicalInputsError) else "FAILED",
                "error_type":type(exc).__name__,"safe_error":messages.get(str(exc),"临床字段、ID或输入检查未通过。"),"run_id":timestamp}
        if isinstance(exc,ClinicalInputsError):result["clinical_audit"]=exc.audit
        write_json(share/"clinical_trajectories.json",result);write_json(output/"status.json",result)
        save_review(result,share)
        return result,output
    write_json(private/"source_files.json",result.pop("private_source_files"))
    csv_rows(private/"selection.csv",result.pop("private_selection_rows"))
    csv_rows(private/"derived_states.csv",result.pop("private_derived_states"))
    result["run_id"]=timestamp
    write_json(share/"clinical_trajectories.json",result)
    for key in ("overall_curve","disabled_curve","primary_distributions","disabled_distributions","recovery_distributions","outcome_rates","transitions"):
        csv_rows(share/(key+".csv"),result[key])
    save_review(result,share)
    write_json(share/"manifest.json",{p.name:{"bytes":p.stat().st_size,"sha256":digest(p)} for p in sorted(share.iterdir())})
    with zipfile.ZipFile(output/"share.zip","x",zipfile.ZIP_DEFLATED) as z:
        for file in sorted(share.iterdir()):z.write(file,arcname=file.name)
    with zipfile.ZipFile(output/"share.zip") as z:
        if z.testzip() is not None:raise RuntimeError("Clinical share ZIP integrity failed")
    write_json(output/"status.json",{"status":result["status"],"run_id":timestamp,"overall_n":result["overall_n"],"disabled_n":result["disabled_n"]})
    return result,output


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=Path(__file__).parent/"config.v2.local.json")
    parser.add_argument("--out",type=Path)
    parser.add_argument("--review",type=Path)
    parser.add_argument("--page",type=int)
    parser.add_argument("--no-pause",action="store_true")
    args=parser.parse_args()
    if args.review:
        if args.out:parser.error("--out不能与--review同时使用")
        result=json.loads(args.review.read_text(encoding="utf8"))
    else:
        if not args.config.is_file():parser.error("配置文件不存在，请用--config指定包含两个SAS路径的配置。")
        result,output=execute(args.config,args.out)
        print("状态:",result["status"]);print("输出:",output.resolve())
        if result["status"].startswith("COMPLETED"):print("请下载:",(output/"share.zip").resolve())
    pages=review_pages(result)
    if args.page is not None and not 1<=args.page<=len(pages):parser.error(f"--page应为1..{len(pages)}")
    chosen=[(args.page-1,pages[args.page-1])] if args.page else list(enumerate(pages))
    for position,(index,(title,lines)) in enumerate(chosen):
        print(f"\n===== PAGE {index+1}/{len(pages)} | {title} =====\n"+"\n".join(lines),flush=True)
        if not args.no_pause and sys.stdin.isatty() and position<len(chosen)-1:input("截图后按Enter继续：")
    return 0 if result["status"].startswith("COMPLETED") else 2


if __name__=="__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""只读提取急性梗死ROI与mRS群体汇总；患者级清单不进入分享包。"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shlex
import sys
import textwrap
import zipfile

from inputs import load_clinical, scan_masks
from mrs_summary import MONTHS, STATES, derive_states, aggregate_groups
from roi_summary import load_atlas, analyse_masks


def public(value):
    if isinstance(value, dict):
        return {key: public(val) for key, val in value.items() if not str(key).startswith("private_")}
    if isinstance(value, (tuple, list, set)):
        return [public(val) for val in value]
    return value


def write_json(path: Path, data) -> None:
    path.write_text(json.dumps(public(data), ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def csv_rows(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        if not rows:
            handle.write("")
            return
        keys = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def source_fingerprint(path: str) -> dict:
    p = Path(path)
    info = p.stat()
    return {"path": str(p), "bytes": info.st_size, "mtime_ns": info.st_mtime_ns, "sha256": digest(p)}


def compact(data) -> str:
    return json.dumps(data, ensure_ascii=False, separators=(",", ":"))


def mrs_table(rows: list[dict], group: str) -> list[str]:
    selected = {(r["month"], r["state"]): r for r in rows if r["group"] == group}
    header = "Month " + " ".join(f"{state:>6}" for state in STATES[:7]) + "   Miss Invalid Conflict"
    result = [header]
    for month in MONTHS:
        result.append(f"{month:>5} " + " ".join(f"{selected[(month, s)]['n']:>6}" for s in STATES[:7])
                      + " " + " ".join(f"{selected[(month, s)]['n']:>7}" for s in STATES[7:]))
    return result


def review_pages(report: dict) -> list[tuple[str, list[str]]]:
    sections = [("运行与定义", [
        "status: " + report["status"], "run: " + report["run_id"],
        "病灶目录: " + report["config"]["lesion_root"],
        "死亡仅由同次mRS=6定义；不读取独立死亡标志。",
        "缺失不前填，包括6分之后的缺失；6分后0-5仅该访视标冲突。",
        "3/12月源字段没有6分时，图中死亡为0不代表没有真实死亡。",
        "当前为描述性展示，不运行回归、插补、配准或分割。",
        "issues: " + compact(report.get("issues", [])),
        "software: " + compact(report["software"]),
    ])]
    clinical = report.get("clinical", {})
    sections.append(("临床来源与连接", [
        *[f"{key}: {compact(value)}" for key, value in clinical.items() if key not in {"sources", "source_files"}],
        *[f"{source}: {path}" for source, path in clinical.get("source_files", {}).items()],
        *[f"{source} ID: {compact({k:v for k,v in info.items() if k not in {'schema','path','required_fields'}})}"
          for source, info in clinical.get("sources", {}).items()],
    ]))
    inventory = report.get("mask_inventory", {})
    sections.append(("掩膜文件与ID匹配", [f"{key}: {compact(value)}" for key, value in inventory.items()]))
    if "roi" not in report:
        sections.append(("尚未统计", ["输入问题已保存。先核对以上页，不会把未完成结果当0。",
                                         "若为未知文件名，请发匿名命名格式页；无需发送患者ID。"] ))
    else:
        roi = report["roi"]
        sections.append(("有效掩膜与空间", [
            *[f"{key}: {compact(roi[key])}" for key in ("input_n", "valid_n", "excluded_n", "exclusion_counts", "orientation_reordered_n", "units_assumed_mm_n")],
            "未记录单位但与毫米图谱网格完全相符者按mm解释，人数如上。",
            "网格相容不等于验证了解剖配准。病例排除明细只留在private目录。",
            "诊断分布（不据此排除）: " + compact(report.get("diagnosis_counts", {})),
            "覆盖: " + compact(roi["coverage_summary"]),
        ]))
        threshold=report["config"]["sensitivity_min_overlap_ml"]
        lines = [f"主定义：any>0；敏感性：重叠>={threshold:g}mL；两者分母相同。",
                 f"ROI                  any_n     any_%   >={threshold:g}mL_n   N"]
        for row in sorted(roi["roi_rows"], key=lambda r: (-r["any_overlap_n"], r["label_id"])):
            percent = "NA" if row["any_overlap_proportion"] is None else f"{100*row['any_overlap_proportion']:.2f}"
            lines.append(f"{row['roi_name']:<21} {row['any_overlap_n']:>6} {percent:>9} {row['n']:>11} {row['denominator_n']:>6}")
        lines += ["同一患者可属于多个ROI，不把脑区人数或比例相加。",
                  "图谱为MNI9＋HO脑干的自定义展示分区，双侧合并；未覆盖比例另报。"]
        sections.append(("脑区累及比例", lines))
    if "primary" in report:
        primary = report["primary"]
        names = {r["roi_id"]: r["roi_name"] for r in report["roi"]["roi_rows"]}
        for group in primary["groups"]:
            title = "全体有效匹配患者" if group["group"] == "overall" else names[group["group"]]
            sections.append((f"{title} mRS分布 N={group['n']}", [
                *mrs_table(primary["distributions"], group["group"]),
                "每一行的全部状态合计为固定N；Miss=缺失，Invalid=非法，Conflict=6后活着评分。",
            ]))
        audit = report["mrs_audit"]
        sections.append(("转移核对与数据限制", [
            "人数/转移边际验证: " + compact(primary["validation"]),
            "原始6后0-5患者数: " + str(audit["death_to_alive_patient_n"]),
            "原始6后0-5访视数: " + compact(audit["death_to_alive_visit_counts"]),
            "全部相邻逐人转移计数已保存到CSV和ppt_data.json，未由边际分布拼接。",
            "无法观察的死亡仍可能藏在缺失中，不能把mRS6比例解释成完整死亡率。",
            "未覆盖ROI的有效病灶仍在总体分母中。",
        ]))
    sections.append(("交付文件", [
        "截图：发送全部PAGE页面用于核对。",
        "精确PPT数值：优先下载并发送本次share.zip，只有群体汇总，不含患者ID。",
        "share/ppt_data.json：总体及ROI mRS分布、转移与阈值敏感性版本。",
        "share/*.csv：可直接查看的脑区统计、mRS与转移表。",
        "private/包含受试者ID和病灶路径，仅留在工作站，不在share.zip内。",
    ]))
    pages = []
    for title, lines in sections:
        wrapped = [part for line in lines for part in (textwrap.wrap(line, 110) or [""])]
        for start in range(0, len(wrapped), 24):
            pages.append((title, wrapped[start:start+24]))
    return pages


def save_review(report: dict, share: Path) -> None:
    pages = review_pages(report)
    all_text = []
    for index, (title, lines) in enumerate(pages, 1):
        value = f"===== PAGE {index}/{len(pages)} | {title} =====\n" + "\n".join(lines) + "\n"
        (share / f"page_{index:02d}.txt").write_text(value, encoding="utf-8")
        all_text.append(value)
    (share / "screenshots.txt").write_text("\n".join(all_text), encoding="utf-8")


def display(report: dict, page: int | None, pause: bool) -> None:
    pages = review_pages(report)
    if page is not None and not 1 <= page <= len(pages):
        raise ValueError(f"page must be 1..{len(pages)}")
    selected = [(page-1,pages[page-1])] if page else list(enumerate(pages))
    for i, (index, (title, lines)) in enumerate(selected):
        print(f"\n===== PAGE {index+1}/{len(pages)} | {title} =====\n" + "\n".join(lines), flush=True)
        if pause and i < len(selected)-1:
            input("截图后按Enter继续：")


def execute(config_path: Path, output: Path | None = None) -> tuple[dict, Path]:
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    base = config_path.resolve().parent
    if cfg.get("death_definition") != "mrs_equal_6_only_at_all_visits" or cfg.get("carry_forward") is not False:
        raise ValueError("This command supports only the confirmed mRS6-only, no carry-forward definition")
    if cfg.get("primary_membership") != "any_positive_overlap" or cfg.get("diagnosis_filter") is not None:
        raise ValueError("Unsupported cohort or primary ROI rule; do not silently alter the requested population")
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = output if output is not None else base / cfg.get("output_root", "outputs") / run_id
    out.mkdir(parents=True, exist_ok=False)
    share, private = out / "share", out / "private"
    share.mkdir(); private.mkdir()
    versions = {"python":sys.version.split()[0]}
    for name in ("numpy", "nibabel", "pyreadstat", "pandas"):
        versions[name] = importlib.metadata.version(name)
    report = {"status":"RUNNING", "run_id":run_id, "config":cfg, "software":versions, "issues":[]}
    write_json(out/"status.json", {"status":"RUNNING", "stage":"clinical", "run_id":run_id})
    try:
        print("1/4 读取指定的两个SAS，仅提取ID、mRS与诊断...", flush=True)
        records, clinical = load_clinical(cfg)
        report["clinical"] = clinical
        if not clinical["clinical_allowed"]:
            report["status"] = "INPUTS_REQUIRED"
        elif cfg["id_mode"] == "numeric_code_n" and not clinical["numeric_id_compatible"]:
            report["status"] = "INPUTS_REQUIRED"
            report["issues"].append("numeric_filename_rule_requires_numeric_code_n_in_both_sources")
        else:
            report["source_fingerprints"] = [source_fingerprint(p) for p in cfg["source_files"].values()]
            print("2/4 扫描全体掩膜并精确连接临床ID...", flush=True)
            mapping, inventory_audit, inventory = scan_masks(cfg["lesion_root"],set(records),cfg["id_regex"],cfg["id_mode"])
            report["mask_inventory"] = inventory_audit
            csv_rows(private/"mask_inventory.csv",inventory)
            if not inventory_audit["mapping_allowed"]:
                report["status"] = "INPUTS_REQUIRED"
            else:
                def resolve(value):
                    path=Path(value)
                    return path if path.is_absolute() else base/path
                atlas = load_atlas(resolve(cfg["atlas_path"]), resolve(cfg["atlas_labels"]))
                provenance=json.loads(resolve(cfg["atlas_provenance"]).read_text(encoding="utf-8"))
                if atlas["sha256"] != provenance["atlas_sha256"] or atlas["labels_sha256"] != provenance["labels_sha256"]:
                    raise ValueError("Frozen atlas or label checksum mismatch")
                report["atlas"] = {**public(atlas),"provenance":provenance}
                print(f"3/4 读取{len(mapping)}张已匹配掩膜，统计ROI重叠...",flush=True)
                roi = analyse_masks(mapping,atlas,float(cfg["sensitivity_min_overlap_ml"]),
                    progress=lambda n,total:print(f"  {n}/{total}",flush=True),
                    allow_unknown_units_if_exact_grid=cfg["allow_unknown_units_if_exact_grid"])
                report["roi"] = public(roi)
                csv_rows(private/"mask_exclusions.csv",[{"patient_id":i,"path":mapping[i],"reason":r}
                                                       for i,r in roi["private_exclusions"].items()])
                if not roi["valid_n"]:
                    report["status"] = "NO_VALID_MASKS"
                else:
                    print("4/4 汇总总体和所有ROI的六次mRS及相邻转移...",flush=True)
                    valid_records={i:records[i] for i in roi["private_valid_ids"]}
                    states,audit=derive_states(valid_records)
                    primary=aggregate_groups(states,roi["private_any_memberships"])
                    sensitivity=aggregate_groups(states,roi["private_memberships"])
                    diagnoses=Counter(str(row.get("D_DIAG")) for row in valid_records.values())
                    # 诊断原始编码只统计1/2/缺失/其他，避免任意字符串进入分享资料。
                    report["diagnosis_counts"]={"infarction_1":sum(v for k,v in diagnoses.items() if k in {"1","1.0"}),
                        "TIA_2":sum(v for k,v in diagnoses.items() if k in {"2","2.0"}),
                        "other_or_missing":sum(v for k,v in diagnoses.items() if k not in {"1","1.0","2","2.0"})}
                    report.update(primary=primary,sensitivity=sensitivity,mrs_audit=audit)
                    csv_rows(private/"derived_states.csv",[{"patient_id":i,**{f"m{m}":s for m,s in visits.items()}} for i,visits in sorted(states.items())])
                    report["status"]="COMPLETED_WITH_FLAGS" if roi["excluded_n"] or audit["has_conflicts"] else "COMPLETED"
                    if roi["excluded_n"]:
                        report["issues"].append("Some matched masks excluded for explicit image-format/space/empty reasons; denominator is valid subset")
                    if audit["has_conflicts"]:
                        report["issues"].append("Alive mRS after prior mRS6 is retained as conflict, not a recovery transition")
                    for title, result in (("primary",primary),("sensitivity",sensitivity)):
                        for table in ("distributions","transitions","state3_distributions","state3_transitions","change"):
                            csv_rows(share/f"{title}_{table}.csv",result[table])
                csv_rows(share/"roi_summary.csv",roi["roi_rows"])
    except Exception as exc:
        report["status"]="FAILED"
        report["issues"].append({"error_type":type(exc).__name__,"detail":"See private/failure.txt on workstation"})
        import traceback
        (private/"failure.txt").write_text(traceback.format_exc(),encoding="utf-8")
    report=public(report)
    write_json(share/"ppt_data.json",report)
    write_json(out/"status.json",{"status":report["status"],"run_id":run_id,"issues":report["issues"]})
    save_review(report,share)
    files=sorted(p for p in share.iterdir() if p.is_file())
    write_json(share/"manifest.json",{p.name:{"sha256":digest(p),"bytes":p.stat().st_size} for p in files})
    with zipfile.ZipFile(out/"share.zip","x",zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(share.iterdir()):
            archive.write(path,arcname=path.name)
    with zipfile.ZipFile(out/"share.zip") as archive:
        if archive.testzip() is not None:
            raise RuntimeError("Share archive integrity failure")
    return report,out


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=Path(__file__).parent/"config.local.json")
    parser.add_argument("--out",type=Path)
    parser.add_argument("--review",type=Path,help="重放已有share/ppt_data.json，不重读原数据")
    parser.add_argument("--page",type=int)
    parser.add_argument("--no-pause",action="store_true")
    args=parser.parse_args()
    if args.review:
        report=json.loads(args.review.read_text(encoding="utf-8"))
    else:
        if not args.config.is_file():
            parser.error("本地配置不存在；先运行 configure_workstation.py --data-root <项目父目录>。")
        report,out=execute(args.config,args.out)
        print(f"\n状态: {report['status']}\n输出: {out.resolve()}\n请下载: {(out/'share.zip').resolve()}",flush=True)
        print("重看第1页:\npython "+shlex.quote(str(Path(__file__).resolve()))+" --review "+shlex.quote(str((out/'share/ppt_data.json').resolve()))+" --page 1",flush=True)
    display(report,args.page,not args.no_pause and sys.stdin.isatty())
    return 0 if report["status"].startswith("COMPLETED") else 2


if __name__=="__main__":
    raise SystemExit(main())

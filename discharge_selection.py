"""从冻结早期SAS读取出院筛选字段；仅将匿名汇总交给分享报告。"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from inputs import _id
from mrs_summary import MONTHS, STATES, parse_mrs
from run_summary import digest

REQUIRED_FIELDS = ("code_n", "D_DIAG", "D_MRS")
DISCHARGE_INCLUDED = ("2", "3", "4", "5", "6")
PARSED_VALUES = tuple(str(n) for n in range(7)) + ("missing", "invalid")


def load_discharge_records(report: dict) -> tuple[dict[str, dict[str, str]], dict]:
    """仅从原早期SAS读取三个字段；SHA256、字段或ID检查失败均阻断。

    返回records含患者ID，只能用于工作站内存连接或private文件。
    返回audit不含患者ID、输入路径或任意患者原始字段值。
    """
    cfg = report["config"]
    width = cfg.get("id_width")
    if width is not None and (isinstance(width, bool) or not isinstance(width, int) or width < 1):
        raise ValueError("Invalid original id_width")
    source = Path(cfg["source_files"]["early"]).expanduser().resolve()
    if not source.is_file():
        raise ValueError("Original early SAS is unavailable")
    fingerprints = [row for row in report.get("source_fingerprints", [])
                    if Path(row.get("path", "")).expanduser().resolve() == source]
    if len(fingerprints) != 1 or not fingerprints[0].get("sha256"):
        raise ValueError("Original early SAS SHA256 fingerprint is missing or ambiguous")
    frozen = fingerprints[0]
    before_sha = digest(source)
    before_bytes = source.stat().st_size
    if before_sha != frozen["sha256"] or before_bytes != frozen.get("bytes"):
        raise ValueError("Early SAS differs from the original run fingerprint")
    original_source = report["clinical"]["sources"]["early"]
    expected_storage = original_source.get("id_storage")
    if expected_storage not in ("numeric", "character"):
        raise ValueError("Original early SAS ID storage audit is missing")

    import pyreadstat

    _, meta = pyreadstat.read_sas7bdat(str(source), metadataonly=True,
                                     encoding=cfg.get("sas_encoding"))
    names = list(meta.column_names)
    spelling = {name.casefold(): name for name in names}
    frequencies = Counter(name.casefold() for name in names)
    if any(frequencies[field.casefold()] != 1 for field in REQUIRED_FIELDS):
        raise ValueError("Required fields code_n, D_DIAG and D_MRS must each exist exactly once")
    actual = {field: spelling[field.casefold()] for field in REQUIRED_FIELDS}
    types = getattr(meta, "readstat_variable_types", {}) or {}
    labels = getattr(meta, "column_names_to_labels", {}) or {}
    frame, _ = pyreadstat.read_sas7bdat(str(source), usecols=list(actual.values()),
                                      encoding=cfg.get("sas_encoding"),
                                      disable_datetime_conversion=True, user_missing=True)
    after_sha = digest(source)
    if after_sha != before_sha or source.stat().st_size != before_bytes:
        raise ValueError("Early SAS changed while discharge fields were being read")
    if not len(frame):
        raise ValueError("Early SAS contains no records")
    if len(frame) != original_source.get("rows"):
        raise ValueError("Early SAS row count differs from the original audit")
    dtype = str(frame[actual["code_n"]].dtype)
    id_type = types.get(actual["code_n"])
    numeric = id_type == "double" or (id_type is None and dtype.startswith(
        ("float", "int", "uint", "Float", "Int", "UInt")))
    storage = "numeric" if numeric else "character"
    if storage != expected_storage:
        raise ValueError("Early SAS ID storage differs from the original audit")
    columns = {field: frame[name].tolist() for field, name in actual.items()}
    ids: list[str] = []
    reasons: Counter = Counter()
    for raw in columns["code_n"]:
        identifier, reason = _id(raw, numeric, width)
        if reason or identifier is None:
            reasons[reason or "invalid"] += 1
        else:
            ids.append(identifier)
    counts = Counter(ids)
    duplicates = sum(n > 1 for n in counts.values())
    if reasons or duplicates:
        raise ValueError("Early SAS has missing, invalid or duplicate normalized patient IDs")
    records = {identifier: {field: parse_mrs(columns[field][index])
                            for field in ("D_DIAG", "D_MRS")}
               for index, identifier in enumerate(ids)}
    audit = {
        "source_role": "early", "source_sha256": before_sha, "source_bytes": before_bytes,
        "original_sha256_matched": True, "unchanged_during_read": True,
        "required_fields": list(REQUIRED_FIELDS),
        "schema": {field: {"actual": name, "type": types.get(name), "label": labels.get(name)}
                   for field, name in actual.items()},
        "rows": len(frame), "id_storage": storage, "id_source_type": id_type,
        "id_pandas_dtype": dtype, "id_unique_n": len(counts), "id_missing_n": 0,
        "id_invalid_n": 0, "id_duplicate_n": 0,
        "id_outer_whitespace_n": sum(identifier != identifier.strip() for identifier in ids),
        "independent_death_fields_read": [],
    }
    return records, audit


def _checked_record(row: dict) -> tuple[str, str]:
    diagnosis, discharge = row.get("D_DIAG"), row.get("D_MRS")
    if diagnosis not in PARSED_VALUES or discharge not in PARSED_VALUES:
        raise ValueError("Discharge records must contain canonical D_DIAG and D_MRS values")
    return diagnosis, discharge


def select_discharge_cohort(valid_states: dict, records: dict) -> tuple[set[str], dict]:
    """在原有效影像人群内筛选D_DIAG=1且D_MRS为2–6，不依赖随访完整性。"""
    excluded = {key: 0 for key in (
        "no_early_record", "diagnosis_tia_2", "diagnosis_missing", "diagnosis_invalid",
        "diagnosis_other", "ischemic_discharge_mrs_missing", "ischemic_discharge_mrs_invalid",
        "ischemic_discharge_mrs_0_1")}
    included_counts = {score: 0 for score in DISCHARGE_INCLUDED}
    diagnosis_counts = {key: 0 for key in ("infarction_1", "TIA_2", "missing", "invalid", "other", "no_early_record")}
    selected: set[str] = set()
    ischemic_n = 0
    for identifier in valid_states:
        if identifier not in records:
            excluded["no_early_record"] += 1
            diagnosis_counts["no_early_record"] += 1
            continue
        diagnosis, discharge = _checked_record(records[identifier])
        if diagnosis != "1":
            key = "tia_2" if diagnosis == "2" else diagnosis if diagnosis in ("missing", "invalid") else "other"
            excluded["diagnosis_" + key] += 1
            diagnosis_counts["TIA_2" if diagnosis == "2" else key] += 1
            continue
        diagnosis_counts["infarction_1"] += 1
        ischemic_n += 1
        if discharge in DISCHARGE_INCLUDED:
            selected.add(identifier)
            included_counts[discharge] += 1
        else:
            key = discharge if discharge in ("missing", "invalid") else "0_1"
            excluded["ischemic_discharge_mrs_" + key] += 1
    conserved = (len(selected) + sum(excluded.values()) == len(valid_states)
                 and sum(included_counts.values()) == len(selected)
                 and sum(diagnosis_counts.values()) == len(valid_states))
    if not conserved:
        raise RuntimeError("Discharge cohort flow conservation failed")
    return selected, {
        "original_valid_n": len(valid_states), "early_record_n": len(valid_states) - excluded["no_early_record"],
        "ischemic_n": ischemic_n, "included_n": len(selected),
        "mutually_exclusive_exclusions": excluded, "diagnosis_counts": diagnosis_counts,
        "included_discharge_mrs_counts": included_counts,
        "discharge_mrs6_included_n": included_counts["6"],
        "flow_conservation_passed": True,
        "rule": "D_DIAG=1 and integer D_MRS in 2..6 within the original valid nonempty mask cohort",
        "followup_completeness_required": False,
    }


def apply_discharge_death_conflicts(states: dict, records: dict) -> tuple[dict, dict]:
    """出院mRS6后的0–5仅标记当次冲突；不填补缺失，不改变原文件。"""
    updated: dict[str, dict[int, str]] = {}
    by_month = {str(month): 0 for month in MONTHS}
    patients = 0
    discharge_deaths = 0
    for identifier, visits in states.items():
        if set(visits) != set(MONTHS) or any(state not in STATES for state in visits.values()):
            raise ValueError("Each selected patient must have six validated follow-up states")
        if identifier not in records:
            raise ValueError("A selected patient lacks a discharge record")
        diagnosis, discharge = _checked_record(records[identifier])
        if diagnosis != "1" or discharge not in DISCHARGE_INCLUDED:
            raise ValueError("Discharge conflict checks require the selected ischemic cohort")
        copy = dict(visits)
        patient_changed = False
        if discharge == "6":
            discharge_deaths += 1
            for month in MONTHS:
                if copy[month] in ("0", "1", "2", "3", "4", "5"):
                    copy[month] = "conflict"
                    by_month[str(month)] += 1
                    patient_changed = True
        patients += patient_changed
        updated[identifier] = copy
    return updated, {
        "discharge_mrs6_n": discharge_deaths, "new_conflict_patient_n": patients,
        "new_conflict_visit_n": sum(by_month.values()), "new_conflict_visit_counts": by_month,
        "missing_after_discharge6": "retain_missing", "carry_forward": False,
        "rule": "Observed follow-up mRS0-5 after discharge mRS6 becomes conflict; all other states are retained",
    }

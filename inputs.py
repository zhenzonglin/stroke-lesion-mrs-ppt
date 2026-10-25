"""固定来源的mRS读取与全体病灶目录清单；不读取任何死亡标志字段。"""
from __future__ import annotations

from collections import Counter
import math
import os
from pathlib import Path
import re
import stat
from typing import Any


SOURCE_FIELDS = {
    "early": ("code_n", "F3_MRS", "F12_MRS", "D_DIAG"),
    "late": ("code_n", "m24_mrs", "m36_mrs", "m48_mrs", "m60_mrs"),
}
DEFAULT_ID_REGEX = r"^(?:sub-)?(?P<id>[0-9]+)(?:$|[_-])"


def _id(value: Any, numeric: bool, width: int | None) -> tuple[str | None, str | None]:
    """数字ID只接受可精确表示的非负整数；字符ID保留大小写及前导零。"""
    if value is None or (isinstance(value, str) and not value.strip()):
        return None, "missing"
    if numeric:
        if isinstance(value, str) and re.fullmatch(r"\.?[A-Z_]|\.", value.strip()):
            return None, "missing"
        try:
            number = float(value)
        except (TypeError, ValueError, OverflowError):
            return None, "invalid"
        if math.isnan(number):
            return None, "missing"
        if not math.isfinite(number) or number < 0 or not number.is_integer() or number > 2**53 - 1:
            return None, "invalid"
        text = str(int(number))
        return (text.zfill(width) if width else text), None
    if isinstance(value, str):
        # 前后空格保留，不偷偷修复患者主键。
        return value, None
    try:
        if math.isnan(value):
            return None, "missing"
    except (TypeError, ValueError):
        pass
    return None, "invalid"


def load_clinical(config: dict, *, include_discharge_mrs: bool = False) -> tuple[dict[str, dict], dict]:
    """读取两个SAS并保留ID并集；显式请求时才从early读取出院D_MRS。"""
    required_by_source = dict(SOURCE_FIELDS)
    if include_discharge_mrs:
        required_by_source["early"] = (*required_by_source["early"], "D_MRS")
    audit = {"clinical_allowed": False, "numeric_id_compatible": False,
             "source_files": {}, "sources": {}, "errors": [],
             "id_rule": "exact_character_or_integer_numeric_explicit_width",
             "death_fields_read": [], "join": "outer_union",
             "early_n": 0, "late_n": 0, "early_only_n": 0, "late_only_n": 0,
             "early_missing_late_n": 0}
    width = config.get("id_width")
    if width is not None and (isinstance(width, bool) or not isinstance(width, int) or width < 1):
        audit["errors"].append({"stage": "config", "error": "invalid_id_width"})
        return {}, audit
    source_files = config.get("source_files", {})
    if not isinstance(source_files, dict) or any(not source_files.get(key) for key in SOURCE_FIELDS):
        audit["errors"].append({"stage": "config", "error": "source_files_requires_early_and_late"})
        return {}, audit
    try:
        import pyreadstat
    except ImportError:
        audit["errors"].append({"stage": "dependency", "error": "pyreadstat_unavailable"})
        return {}, audit
    payloads = {}
    for source, required in required_by_source.items():
        path = Path(source_files[source]).expanduser().resolve()
        audit["source_files"][source] = str(path)
        info = {"path": str(path), "required_fields": list(required), "schema": {},
                "rows": 0, "status": "not_read"}
        audit["sources"][source] = info
        try:
            if not path.is_file():
                raise FileNotFoundError("configured_source_missing")
            _, meta = pyreadstat.read_sas7bdat(str(path), metadataonly=True,
                                             encoding=config.get("sas_encoding"))
            names = list(meta.column_names)
            lower_counts = Counter(name.casefold() for name in names)
            spelling = {name.casefold(): name for name in names}
            missing = [field for field in required if field.casefold() not in spelling]
            duplicate_columns = [field for field in required if lower_counts[field.casefold()] > 1]
            if missing or duplicate_columns:
                info.update(status="field_blocked", missing_fields=missing,
                            duplicate_columns=duplicate_columns)
                audit["errors"].append({"source": source, "stage": "fields", "error": "required_fields_not_unique"})
                continue
            types = getattr(meta, "readstat_variable_types", {}) or {}
            labels = getattr(meta, "column_names_to_labels", {}) or {}
            actual = {field: spelling[field.casefold()] for field in required}
            info["schema"] = {field: {"actual": name, "type": types.get(name), "label": labels.get(name)}
                              for field, name in actual.items()}
            frame, _ = pyreadstat.read_sas7bdat(str(path), usecols=list(actual.values()),
                                              encoding=config.get("sas_encoding"),
                                              disable_datetime_conversion=True, user_missing=True)
            dtype = str(frame[actual["code_n"]].dtype)
            id_type = types.get(actual["code_n"])
            numeric = (id_type == "double" or (id_type is None and dtype.startswith(
                ("float", "int", "uint", "Float", "Int", "UInt"))))
            columns = {field: frame[name].tolist() for field, name in actual.items()}
            ids, reasons = [], Counter()
            for value in columns["code_n"]:
                identifier, reason = _id(value, numeric, width)
                ids.append(identifier)
                if reason:
                    reasons[reason] += 1
            counts = Counter(value for value in ids if value is not None)
            duplicates = sum(count > 1 for count in counts.values())
            info.update(rows=len(frame), id_storage="numeric" if numeric else "character",
                        id_source_type=id_type, id_pandas_dtype=dtype,
                        id_unique_n=len(counts), id_missing_n=reasons["missing"],
                        id_invalid_n=reasons["invalid"], id_duplicate_n=duplicates,
                        id_outer_whitespace_n=sum(value != value.strip() for value in ids if value is not None),
                        status="read")
            if not len(frame):
                info["status"] = "empty_source"
                audit["errors"].append({"source": source, "stage": "read", "error": "empty_clinical_source"})
                continue
            if reasons or duplicates:
                info["status"] = "id_blocked"
                audit["errors"].append({"source": source, "stage": "id", "error": "invalid_or_duplicate_ids"})
                continue
            payloads[source] = {identifier: {field: values[index] for field, values in columns.items()
                                              if field != "code_n"}
                                for index, identifier in enumerate(ids)}
        except Exception as exc:
            info["status"] = "read_error"
            audit["errors"].append({"source": source, "stage": "read", "error": type(exc).__name__})
    audit["numeric_id_compatible"] = (len(audit["sources"]) == 2 and all(
        source.get("id_storage") == "numeric" for source in audit["sources"].values()))
    if audit["errors"] or set(payloads) != set(SOURCE_FIELDS):
        return {}, audit
    early, late = payloads["early"], payloads["late"]
    # 任一来源存在的患者均保留；未在另一来源出现时，其随访字段保持缺席。
    records = {identifier: {**early.get(identifier, {}), **late.get(identifier, {})}
               for identifier in sorted(early.keys() | late.keys())}
    audit.update(clinical_allowed=True, early_n=len(early), late_n=len(late),
                 early_only_n=len(early.keys() - late.keys()),
                 late_only_n=len(late.keys() - early.keys()),
                 early_missing_late_n=len(early.keys() - late.keys()), joined_n=len(records))
    return records, audit


def _anonymous_pattern(name: str) -> str:
    """仅保留ASCII技术词及._-分隔符；姓名和所有非ASCII字符均匿名化。"""
    known = {"sub", "lesion", "mask", "space", "mni", "desc", "acute", "thr",
             "flair", "t", "w", "label", "to", "mm", "nii", "gz"}
    parts = []
    for token in re.findall(r"[A-Za-z]+|[0-9]+|[_.-]+|[^A-Za-z0-9_.-]+", name):
        if re.fullmatch(r"[_.-]+", token):
            parts.append(token)
        elif re.fullmatch(r"[0-9]+", token):
            parts.append("<n>")
        elif token.lower() in known:
            parts.append(token.lower())
        else:
            parts.append("<text>")
    return re.sub(r"(?:<text>)+", "<text>", "".join(parts))


def scan_masks(root: str, clinical_ids: set[str], id_regex: str = DEFAULT_ID_REGEX,
               id_mode: str = "numeric_code_n") -> tuple[dict[str, str], dict, list[dict]]:
    """扫描所有NIfTI；不允许以静默漏扫或猜测ID构成临床分析子集。

    调用端只有在load_clinical的numeric_id_compatible为True时才可使用numeric模式。
    mapping使用临床主表的规范ID作为键；private_inventory只能保存在工作站。
    """
    directory = Path(root).expanduser().resolve()
    audit = {"root": str(directory), "id_regex": id_regex, "id_mode": id_mode,
             "mapping_allowed": False, "all_file_n": 0, "nifti_n": 0,
             "recognized_n": 0, "unrecognized_n": 0, "matched_n": 0,
             "unmatched_n": 0, "duplicate_id_n": 0, "empty_file_n": 0,
             "symlink_directory_n": 0, "filename_patterns": [], "errors": []}
    inventory: list[dict] = []
    mapping: dict[str, str] = {}
    if not directory.is_dir():
        audit["errors"].append("lesion_root_missing")
        return mapping, audit, inventory
    if id_mode not in {"numeric_code_n", "character"}:
        audit["errors"].append("unknown_id_mode")
        return mapping, audit, inventory
    try:
        expression = re.compile(id_regex)
        if "id" not in expression.groupindex:
            raise ValueError("regex_requires_named_id_group")
    except (re.error, ValueError):
        audit["errors"].append("invalid_id_regex")
        return mapping, audit, inventory
    clinical_lookup = {}
    for identifier in clinical_ids:
        if id_mode == "numeric_code_n":
            if not isinstance(identifier, str) or not re.fullmatch(r"[0-9]+", identifier):
                audit["errors"].append("clinical_ids_incompatible_with_numeric_mode")
                return mapping, audit, inventory
            canonical = str(int(identifier))
        else:
            canonical = identifier
        if canonical in clinical_lookup:
            audit["errors"].append("clinical_id_collision_after_numeric_normalization")
            return mapping, audit, inventory
        clinical_lookup[canonical] = identifier
    groups: dict[str, list[dict]] = {}
    patterns = Counter()

    def scan_error(exc: OSError) -> None:
        raise exc

    try:
        for folder, directories, files in os.walk(directory, onerror=scan_error, followlinks=False):
            directories.sort()
            symlinks = [name for name in directories if (Path(folder) / name).is_symlink()]
            audit["symlink_directory_n"] += len(symlinks)
            for filename in sorted(files):
                path = Path(folder) / filename
                if not filename.casefold().endswith((".nii", ".nii.gz")):
                    audit["all_file_n"] += 1
                    continue
                info = path.stat()
                if not stat.S_ISREG(info.st_mode):
                    continue
                audit["all_file_n"] += 1
                audit["nifti_n"] += 1
                stem = filename[:-7] if filename.casefold().endswith(".nii.gz") else filename[:-4]
                match = expression.match(stem)
                row = {"id": None, "raw_id": None, "clinical_id": None,
                       "path": str(path.absolute()), "status": "unrecognized_name", "bytes": info.st_size}
                inventory.append(row)
                if info.st_size == 0:
                    audit["empty_file_n"] += 1
                    row["status"] = "empty_file"
                if not match or match.group("id") is None or not match.group("id"):
                    audit["unrecognized_n"] += 1
                    patterns[_anonymous_pattern(filename)] += 1
                    continue
                raw_id = match.group("id")
                if id_mode == "numeric_code_n" and not re.fullmatch(r"[0-9]+", raw_id):
                    audit["unrecognized_n"] += 1
                    patterns[_anonymous_pattern(filename)] += 1
                    continue
                identifier = str(int(raw_id)) if id_mode == "numeric_code_n" else raw_id
                audit["recognized_n"] += 1
                row.update(id=identifier, raw_id=raw_id, clinical_id=clinical_lookup.get(identifier))
                row["status"] = ("empty_file" if info.st_size == 0 else
                                 "matched" if identifier in clinical_lookup else "no_clinical_match")
                groups.setdefault(identifier, []).append(row)
    except OSError as exc:
        audit["errors"].append("directory_scan_" + type(exc).__name__)
    for identifier, rows in groups.items():
        if len(rows) > 1:
            audit["duplicate_id_n"] += 1
            for row in rows:
                row["status"] = "duplicate_id"
        elif rows[0]["clinical_id"] is not None and rows[0]["status"] == "matched":
            mapping[rows[0]["clinical_id"]] = rows[0]["path"]
    audit["matched_n"] = len(mapping)
    audit["unmatched_n"] = sum(row["clinical_id"] is None and row["id"] is not None for row in inventory)
    audit["filename_patterns"] = [{"pattern": key, "n": count} for key, count in patterns.most_common(10)]
    for field, error in (("duplicate_id_n", "duplicate_mask_ids"),
                         ("unrecognized_n", "unrecognized_nifti_names"),
                         ("empty_file_n", "empty_nifti_files"),
                         ("symlink_directory_n", "symlink_directories_not_scanned")):
        if audit[field]:
            audit["errors"].append(error)
    if not audit["nifti_n"]:
        audit["errors"].append("no_nifti_files")
    if not mapping:
        audit["errors"].append("no_clinical_matches")
    audit["mapping_allowed"] = not audit["errors"]
    return mapping, audit, inventory

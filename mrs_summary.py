"""仅按原始mRS定义的固定人群随访汇总；不读取任何独立死亡字段。"""
from __future__ import annotations

from collections import Counter
import math
import re
from typing import Any

MONTHS = (3, 12, 24, 36, 48, 60)
MRS_FIELDS = {3: "F3_MRS", 12: "F12_MRS", **{m: f"m{m}_mrs" for m in MONTHS[2:]}}
STATES = tuple(str(value) for value in range(7)) + ("missing", "invalid", "conflict")
RAW_STATES = STATES[:-1]
STATE3 = ("independent", "dependent", "death", "missing", "invalid", "conflict")
CHANGE = ("improved", "stable", "worsened", "to_death", "remaining_death",
          "death_after_missing", "missing_after_death", "missing", "invalid", "conflict")
_SPECIAL_MISSING = re.compile(r"^(?:|\.|\.?[A-Z_])$", re.IGNORECASE)


def parse_mrs(value: Any) -> str:
    """接受整数0–6；缺失、非法值分开，不输出任意原始文本。"""
    if value is None:
        return "missing"
    if isinstance(value, str) and _SPECIAL_MISSING.fullmatch(value.strip()):
        return "missing"
    if isinstance(value, bool):
        return "invalid"
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return "invalid"
    if math.isnan(number):
        return "missing"
    if not math.isfinite(number) or not number.is_integer() or not 0 <= number <= 6:
        return "invalid"
    return str(int(number))


def derive_states(records: dict[str, dict]) -> tuple[dict[str, dict[int, str]], dict]:
    """mRS=6为唯一死亡证据；6后0–5标conflict，其他时点保留原观测。

    所有D_DEATH/F*_DEATH字段均被忽略。死亡后缺失仍为missing，不前填；
    先前的6不回溯改写。返回states含ID，仅供内存连接，audit不含ID。
    """
    raw_counts = {month: Counter({state: 0 for state in RAW_STATES}) for month in MONTHS}
    counts = {month: Counter({state: 0 for state in STATES}) for month in MONTHS}
    absent = Counter({month: 0 for month in MONTHS})
    conflict_counts = Counter({month: 0 for month in MONTHS})
    first_death = Counter({month: 0 for month in MONTHS})
    conflict_patients = 0
    states = {}
    for identifier, row in records.items():
        if not isinstance(row, dict):
            raise TypeError("Each clinical record must be a field-value dictionary")
        known_dead, patient_conflict = False, False
        derived = {}
        for month, field in MRS_FIELDS.items():
            absent[month] += field not in row
            raw = parse_mrs(row.get(field))
            raw_counts[month][raw] += 1
            state = raw
            if known_dead and raw in tuple(str(score) for score in range(6)):
                state = "conflict"
                conflict_counts[month] += 1
                patient_conflict = True
            if raw == "6" and not known_dead:
                first_death[month] += 1
                known_dead = True
            derived[month] = state
            counts[month][state] += 1
        states[identifier] = derived
        conflict_patients += patient_conflict
    audit = {
        "n_subjects": len(records), "months": list(MONTHS),
        "fields": {str(month): field for month, field in MRS_FIELDS.items()},
        "rules": {
            "death_evidence": "mrs_equal_6_only_at_all_visits",
            "independent_death_fields_used": [],
            "carry_forward": False,
            "missing_after_mrs6": "retain_missing",
            "conflict": "observed_0_to_5_after_any_prior_6_becomes_conflict_at_that_visit",
            "earlier_mrs6": "retain_original_6",
            "cohort_exclusion": "none_for_missing_invalid_or_conflicting_outcomes",
        },
        "raw_mrs_counts": {str(month): dict(raw_counts[month]) for month in MONTHS},
        "derived_counts": {str(month): dict(counts[month]) for month in MONTHS},
        "source_row_or_field_absent_n": {str(month): absent[month] for month in MONTHS},
        "death_to_alive_patient_n": int(conflict_patients),
        "death_to_alive_visit_n": int(sum(conflict_counts.values())),
        "death_to_alive_visit_counts": {str(month): conflict_counts[month] for month in MONTHS},
        "first_observed_mrs6_counts": {str(month): first_death[month] for month in MONTHS},
        "has_conflicts": bool(conflict_patients),
        "raw_and_derived_counts_preserve_n": all(
            sum(raw_counts[month].values()) == sum(counts[month].values()) == len(records)
            for month in MONTHS),
    }
    return states, audit


def _state3(state: str) -> str:
    if state in ("0", "1", "2"):
        return "independent"
    if state in ("3", "4", "5"):
        return "dependent"
    return "death" if state == "6" else state


def _change(source: str, target: str) -> str:
    if "conflict" in (source, target):
        return "conflict"
    if "invalid" in (source, target):
        return "invalid"
    if source == "missing" and target == "6":
        return "death_after_missing"
    if source == "6" and target == "missing":
        return "missing_after_death"
    if "missing" in (source, target):
        return "missing"
    if target == "6":
        return "remaining_death" if source == "6" else "to_death"
    if source == "6":
        raise ValueError("Derived states must not contain an unflagged death-to-alive transition")
    delta = int(target) - int(source)
    return "improved" if delta < 0 else "worsened" if delta > 0 else "stable"


def _ratio(n: int, denominator: int) -> float | None:
    return n / denominator if denominator else None


def aggregate_groups(states: dict[str, dict[int, str]], groups: dict[str, set[str]]) -> dict:
    """输出聚合表，自动加overall；ROI可重叠，未知ID/错误状态不静默丢弃。

    distributions始终含全部10种状态，三状态表含6种。transitions仅保留
    非零边；其分母为固定组总N，另给起始状态人数与条件比例。change的
    改善/稳定/恶化只比较相邻两次均0–5者，死亡单独分类。
    """
    all_ids = set(states)
    for visit_states in states.values():
        if set(visit_states) != set(MONTHS) or any(state not in STATES for state in visit_states.values()):
            raise ValueError("Each subject must have exactly six validated derived visit states")
        seen_death = False
        for month in MONTHS:
            state = visit_states[month]
            if seen_death and state in ("0", "1", "2", "3", "4", "5"):
                raise ValueError("Unflagged alive score after an earlier mRS6")
            seen_death |= state == "6"
    memberships = {name: set(ids) for name, ids in groups.items()}
    if any(not isinstance(name, str) or not name for name in memberships):
        raise ValueError("Group names must be nonempty strings")
    if any(ids - all_ids for ids in memberships.values()):
        raise ValueError("Group contains IDs absent from derived states")
    if "overall" in memberships and memberships["overall"] != all_ids:
        raise ValueError("Reserved group overall must contain all derived subjects")
    memberships = {"overall": all_ids, **memberships}
    result = {"groups": [], "distributions": [], "transitions": [],
              "state3_distributions": [], "state3_transitions": [], "change": [],
              "validation": {"all_passed": True, "death_to_alive_transitions_n": 0,
                             "fixed_denominator": True, "marginal_checks_n": 0},
              "definitions": {
                  "proportion": "n divided by fixed group denominator, from 0 to 1",
                  "proportion_from_state": "transition n divided by source-state n",
                  "proportion_alive_pairs": "only improved/stable/worsened, denominator adjacent pairs both mRS 0-5",
                  "missing_after_death": "mRS6 followed by missing; no death carry-forward",
                  "death_after_missing": "missing followed by mRS6; death time is not inferred",
                  "roi_overlap": "One participant may enter multiple ROIs; do not sum ROI Ns or proportions",
              }}
    for group, ids in memberships.items():
        n_group = len(ids)
        result["groups"].append({"group": group, "n": n_group})
        counters = {month: Counter(states[identifier][month] for identifier in ids) for month in MONTHS}
        counters3 = {month: Counter(_state3(states[identifier][month]) for identifier in ids) for month in MONTHS}
        for month in MONTHS:
            for state in STATES:
                n = counters[month][state]
                result["distributions"].append({"group": group, "month": month, "state": state,
                    "n": n, "denominator": n_group, "proportion": _ratio(n, n_group)})
            for state in STATE3:
                n = counters3[month][state]
                result["state3_distributions"].append({"group": group, "month": month, "state": state,
                    "n": n, "denominator": n_group, "proportion": _ratio(n, n_group)})
            if sum(counters[month].values()) != n_group or sum(counters3[month].values()) != n_group:
                raise RuntimeError("Visit counts do not preserve the group denominator")
        for previous, current in zip(MONTHS[:-1], MONTHS[1:]):
            pairs = Counter((states[identifier][previous], states[identifier][current]) for identifier in ids)
            pairs3 = Counter((_state3(source), _state3(target)) for source, target in pairs.elements())
            change = Counter(_change(source, target) for source, target in pairs.elements())
            alive_pairs = sum(change[name] for name in ("improved", "stable", "worsened"))
            if sum(pairs.values()) != n_group or sum(pairs3.values()) != n_group or sum(change.values()) != n_group:
                raise RuntimeError("Transition counts do not preserve the group denominator")
            for pair_counts, known_states, before, after, destination in (
                (pairs, STATES, counters[previous], counters[current], "transitions"),
                (pairs3, STATE3, counters3[previous], counters3[current], "state3_transitions"),
            ):
                for state in known_states:
                    if (sum(n for (source, _), n in pair_counts.items() if source == state) != before[state]
                            or sum(n for (_, target), n in pair_counts.items() if target == state) != after[state]):
                        raise RuntimeError("Transition marginal differs from the visit distribution")
                    result["validation"]["marginal_checks_n"] += 2
                for (source, target), n in sorted(pair_counts.items()):
                    result[destination].append({"group": group, "from_month": previous, "to_month": current,
                        "from_state": source, "to_state": target, "n": n, "denominator": n_group,
                        "proportion": _ratio(n, n_group), "from_state_n": before[source],
                        "proportion_from_state": _ratio(n, before[source])})
            for category in CHANGE:
                n = change[category]
                result["change"].append({"group": group, "from_month": previous, "to_month": current,
                    "change": category, "n": n, "denominator": n_group, "proportion": _ratio(n, n_group),
                    "alive_pair_denominator": alive_pairs,
                    "proportion_alive_pairs": _ratio(n, alive_pairs) if category in ("improved", "stable", "worsened") else None})
    roi_memberships = [ids for name, ids in memberships.items() if name != "overall"]
    occurrence = Counter(identifier for ids in roi_memberships for identifier in ids)
    result["validation"].update(groups_n=len(memberships), overall_n=len(all_ids),
        roi_membership_assignments_n=sum(map(len, roi_memberships)),
        patients_in_multiple_roi_groups_n=sum(n > 1 for n in occurrence.values()))
    return result

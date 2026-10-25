"""合成数据验证，不含真实患者资料。"""
import copy
import json
import math
import unittest

from mrs_summary import MONTHS, MRS_FIELDS, STATES, aggregate_groups, derive_states, parse_mrs


def record(scores):
    return {MRS_FIELDS[month]: score for month, score in zip(MONTHS, scores)}


class MRSummaryTest(unittest.TestCase):
    def test_death_flags_ignored_at_all_visits(self):
        rows = {"private_a": {**record([0, 1, 2, 3, 4, 5]), "D_DEATH": 1,
                            "F3_DEATH": 1, "F6_DEATH": 1, "F12_DEATH": 1, "F2Y_DEATH": 1}}
        states, audit = derive_states(rows)
        self.assertEqual(list(states["private_a"].values()), ["0", "1", "2", "3", "4", "5"])
        self.assertEqual(audit["death_to_alive_patient_n"], 0)
        self.assertEqual(audit["rules"]["independent_death_fields_used"], [])

    def test_only_later_alive_visits_become_conflict(self):
        states, audit = derive_states({"private_a": record([2, 6, None, 1, 6, 0])})
        self.assertEqual(list(states["private_a"].values()), ["2", "6", "missing", "conflict", "6", "conflict"])
        self.assertEqual(audit["death_to_alive_patient_n"], 1)
        self.assertEqual(audit["death_to_alive_visit_n"], 2)
        self.assertEqual(audit["raw_mrs_counts"]["36"]["1"], 1)
        self.assertEqual(audit["derived_counts"]["36"]["conflict"], 1)

    def test_mrs6_is_accepted_at_3_and_12_months(self):
        states, _ = derive_states({"a": record([6, None, 6, None, 6, None])})
        self.assertEqual(list(states["a"].values()), ["6", "missing", "6", "missing", "6", "missing"])

    def test_missing_never_carried_forward_or_backward(self):
        states, _ = derive_states({"a": record([None, 2, None, 6, None, None])})
        self.assertEqual(list(states["a"].values()), ["missing", "2", "missing", "6", "missing", "missing"])

    def test_special_missing_invalid_and_source_unchanged(self):
        rows = {"private_a": record([".A", float("nan"), 1.5, 99, "private_text", 0])}
        before = copy.deepcopy(rows)
        states, audit = derive_states(rows)
        self.assertEqual(list(states["private_a"].values()), ["missing", "missing", "invalid", "invalid", "invalid", "0"])
        self.assertTrue(math.isnan(rows["private_a"]["F12_MRS"]))
        self.assertEqual(rows["private_a"]["m48_mrs"], before["private_a"]["m48_mrs"])
        self.assertNotIn("private_a", json.dumps(audit))
        self.assertNotIn("private_text", json.dumps(audit))
        self.assertEqual(parse_mrs(True), "invalid")
        self.assertEqual(parse_mrs(".z"), "missing")

    def test_fixed_denominator_overlap_and_marginals(self):
        states, _ = derive_states({"private_a": record([3, 2, None, 6, None, 0]),
                                  "private_b": record([1, 2, 2, 3, 4, 4]),
                                  "private_c": record([0, 0, 0, 0, 0, 0])})
        summary = aggregate_groups(states, {"roi_a": {"private_a", "private_b"}, "roi_b": {"private_b", "private_c"}, "empty": set()})
        self.assertEqual(summary["validation"]["patients_in_multiple_roi_groups_n"], 1)
        self.assertTrue(summary["validation"]["all_passed"])
        for group, n in [("overall", 3), ("roi_a", 2), ("roi_b", 2), ("empty", 0)]:
            for month in MONTHS:
                rows = [row for row in summary["distributions"] if row["group"] == group and row["month"] == month]
                self.assertEqual(len(rows), len(STATES))
                self.assertEqual(sum(row["n"] for row in rows), n)
                self.assertTrue(all(row["denominator"] == n for row in rows))
            for previous, current in zip(MONTHS[:-1], MONTHS[1:]):
                rows = [row for row in summary["transitions"] if row["group"] == group and row["from_month"] == previous and row["to_month"] == current]
                self.assertEqual(sum(row["n"] for row in rows), n)
        self.assertNotIn("private_a", json.dumps(summary))
        self.assertTrue(all(row["proportion"] is None for row in summary["distributions"] if row["group"] == "empty"))

    def test_change_denominator_excludes_death_missing_invalid_conflict(self):
        states, _ = derive_states({"a": record([3, 1, 1, 1, 1, 1]), "b": record([1, 2, 2, 2, 2, 2]),
                                  "c": record([0, 0, 0, 0, 0, 0]), "d": record([2, 6, None, 6, 6, 6]),
                                  "e": record([None, 6, None, None, None, None])})
        summary = aggregate_groups(states, {})
        first = {row["change"]: row for row in summary["change"] if row["from_month"] == 3}
        self.assertEqual(first["improved"]["alive_pair_denominator"], 3)
        self.assertAlmostEqual(first["improved"]["proportion_alive_pairs"], 1 / 3)
        self.assertEqual(first["to_death"]["n"], 1)
        self.assertEqual(first["death_after_missing"]["n"], 1)
        self.assertIsNone(first["to_death"]["proportion_alive_pairs"])

    def test_unknown_group_ids_or_unflagged_resurrection_rejected(self):
        states, _ = derive_states({"a": record([0, 0, 0, 0, 0, 0])})
        with self.assertRaises(ValueError):
            aggregate_groups(states, {"roi": {"not_in_states"}})
        with self.assertRaises(ValueError):
            aggregate_groups(states, {"overall": set()})
        states["a"][3] = "6"
        states["a"][12] = "missing"
        with self.assertRaises(ValueError):
            aggregate_groups(states, {})

    def test_empty_cohort_is_explicit(self):
        states, audit = derive_states({})
        summary = aggregate_groups(states, {})
        self.assertEqual(audit["n_subjects"], 0)
        self.assertEqual(summary["groups"], [{"group": "overall", "n": 0}])
        self.assertEqual(len(summary["distributions"]), 60)
        self.assertEqual(summary["transitions"], [])

    def test_requested_recovery_groups_and_all_outcome_definitions_use_valid_observations(self):
        records = {str(score): record([score]*6) for score in range(7)}
        records["missing"] = record([None]*6)
        records["invalid"] = record([99]*6)
        records["conflict"] = record([6, 2, 2, 2, 2, 2])
        states, _ = derive_states(records)
        result = aggregate_groups(states, {"roi": set(records), "empty": set()})
        recovered = {row["category"]: row for row in result["recovery_distributions"]
                     if row["group"] == "roi" and row["month"] == 12}
        self.assertEqual([recovered[key]["n"] for key in ("good_recovery", "functional_disability", "functional_loss")], [2, 2, 3])
        self.assertTrue(all(row["valid_n"] == row["denominator"] == 7 for row in recovered.values()))
        self.assertTrue(all(row["total_n"] == 10 for row in recovered.values()))
        rates = {row["definition"]: row for row in result["outcome_rates"]
                 if row["group"] == "roi" and row["month"] == 12}
        self.assertEqual([rates[key]["numerator_n"] for key in ("mrs_2_5", "mrs_2_6", "mrs_2_3")], [4, 5, 2])
        self.assertAlmostEqual(rates["mrs_2_5"]["proportion"], 4/7)  # 含死亡的7人分母，非存活者6人。
        self.assertEqual(rates["mrs_2_6"]["death6_n"], 1)
        self.assertEqual(rates["mrs_2_6"]["conflict_n"], 1)
        self.assertTrue(all(row["proportion"] is None for row in result["outcome_rates"] if row["group"] == "empty"))
        self.assertTrue(all(row["denominator_n"] == 0 for row in result["outcome_rates"] if row["group"] == "empty"))


if __name__ == "__main__":
    unittest.main()

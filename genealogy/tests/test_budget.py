"""Tests for budget-constrained planning and carrying matches through a merge."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from roots import budget as budget_mod
from roots.evidence import ResearchTask
from roots.matches import ingest as ingest_mod
from roots.store import Store
from tests import helpers


def _task(priority, repository, cost_band, record_set="", subject="X"):
    return ResearchTask(
        subject=subject, subject_xref=None, question="q",
        record_set=record_set, repository=repository,
        cost_band=cost_band, priority=priority,
    )


def _temp_store():
    fd, path = tempfile.mkstemp(suffix=".db")
    os.close(fd)
    os.unlink(path)
    return Store(path), path


class TestCostEstimates(unittest.TestCase):
    def test_free_is_free(self):
        c = budget_mod.estimate_cost(_task(0.5, "irishgenealogy.ie / GRONI", "free"))
        self.assertEqual(c.amount, 0.0)
        self.assertEqual(c.kind, "free")

    def test_subscription_is_a_month_not_a_record(self):
        c = budget_mod.estimate_cost(
            _task(0.5, "FindMyPast / Ancestry / FamilySearch", "subscription")
        )
        self.assertEqual(c.kind, "subscription")
        self.assertGreater(c.amount, 5)
        self.assertLess(c.amount, 30)

    def test_old_gro_births_use_the_cheap_digital_image(self):
        c = budget_mod.estimate_cost(
            _task(0.5, "GRO / FreeBMD", "low", "gro birth index and pdf, around 1901")
        )
        self.assertEqual(c.amount, budget_mod.GRO_DIGITAL_IMAGE)

    def test_marriages_have_no_cheap_option(self):
        c = budget_mod.estimate_cost(
            _task(0.5, "GRO / FreeBMD", "low",
                  "gro marriage index and certificate, around 1923")
        )
        self.assertEqual(c.amount, 12.50)
        self.assertIn("no cheaper", c.note)

    def test_scottish_records_are_cheaper_than_english(self):
        scots = budget_mod.estimate_cost(
            _task(0.5, "ScotlandsPeople", "low", "statutory death register"))
        english = budget_mod.estimate_cost(
            _task(0.5, "GRO / FreeBMD", "low", "gro marriage index and certificate"))
        self.assertLess(scots.amount, english.amount)


class TestBudgetPlanning(unittest.TestCase):
    def test_free_tasks_are_always_taken(self):
        tasks = [_task(0.1, "National Archives of Ireland", "free")]
        plan = budget_mod.plan_within_budget(tasks, 0.0)
        self.assertEqual(len(plan.selected), 1)
        self.assertEqual(plan.spent, 0.0)

    def test_one_subscription_fee_covers_every_subscription_task(self):
        tasks = [
            _task(0.5, "FindMyPast / Ancestry / FamilySearch", "subscription")
            for _ in range(10)
        ]
        plan = budget_mod.plan_within_budget(tasks, 20.0)
        self.assertEqual(len(plan.selected), 10)
        self.assertIsNotNone(plan.subscription_taken)
        self.assertLess(
            plan.spent, 15.0,
            "ten subscription lookups cost one month between them, not ten",
        )

    def test_budget_is_never_exceeded(self):
        tasks = [
            _task(0.9, "GRO / FreeBMD", "low", "gro marriage index and certificate")
            for _ in range(20)
        ]
        plan = budget_mod.plan_within_budget(tasks, 30.0)
        self.assertLessEqual(plan.spent, 30.0)
        self.assertEqual(len(plan.selected), 2, "£30 buys two £12.50 certificates")
        self.assertEqual(len(plan.deferred), 18)

    def test_cheap_high_value_records_outrank_expensive_ones(self):
        cheap = _task(0.70, "GRO / FreeBMD", "low", "gro birth index, around 1880",
                      subject="cheap")
        dear = _task(0.75, "GRO / FreeBMD", "low",
                     "gro marriage index and certificate, around 1900", subject="dear")
        plan = budget_mod.plan_within_budget([dear, cheap], 4.0)
        chosen = [t.subject for t, _c in plan.selected]
        self.assertEqual(chosen, ["cheap"])

    def test_a_subscription_is_skipped_when_it_is_poor_value(self):
        # One low-priority subscription task against several cheap, valuable
        # records: the money belongs on the records.
        tasks = [_task(0.05, "FindMyPast / Ancestry / FamilySearch", "subscription")]
        tasks += [
            _task(0.9, "GRO / FreeBMD", "low", f"gro birth index, around 18{60 + i}")
            for i in range(3)
        ]
        plan = budget_mod.plan_within_budget(tasks, 10.0)
        self.assertIsNone(plan.subscription_taken)
        self.assertEqual(len(plan.selected), 3)

    def test_zero_budget_still_reports_what_is_deferred(self):
        tasks = [_task(0.9, "GRO / FreeBMD", "low", "gro marriage certificate")]
        plan = budget_mod.plan_within_budget(tasks, 0.0)
        self.assertEqual(plan.selected, [])
        self.assertEqual(len(plan.deferred), 1)


class TestLedger(unittest.TestCase):
    def setUp(self):
        self.store, self.path = _temp_store()

    def tearDown(self):
        helpers.cleanup(self.store, self.path)

    def test_spending_accumulates_against_a_cap(self):
        from roots import evidence as ev

        budget_mod.set_cap(self.store, 50.0)
        ev.add_source(self.store, "cert", cost=12.5)
        ev.add_source(self.store, "image", cost=3.0)
        led = budget_mod.ledger(self.store)
        self.assertAlmostEqual(led.spent, 15.5)
        self.assertAlmostEqual(led.remaining, 34.5)
        self.assertEqual(led.records, 2)

    def test_no_cap_means_no_remaining(self):
        led = budget_mod.ledger(self.store)
        self.assertIsNone(led.cap)
        self.assertIsNone(led.remaining)


class TestCarryMatches(unittest.TestCase):
    """Merging kits must not strand the match lists on the old ones."""

    def setUp(self):
        self.store, self.path = _temp_store()
        self.a = self.store.create_kit("a", "self", "23andMe", "37", "-")
        self.b = self.store.create_kit("b", "self", "Ancestry", "37", "-")
        self.merged = self.store.create_kit("m", "self", "merged", "37", "-")
        m1 = self.store.upsert_match(self.a, "23andMe", "x1", name="Cousin One",
                                     total_cm=800, side="paternal")
        m2 = self.store.upsert_match(self.a, "23andMe", "x2", name="Cousin Two",
                                     total_cm=300)
        self.store.replace_match_segments(m1, [("1", 1000, 90000000, 40.0, 5000)])
        self.store.add_shared_match(self.a, m1, m2, 60.0)
        self.store.upsert_match(self.b, "Ancestry", "y1", name="Cousin Three",
                                total_cm=150)
        self.store.commit()

    def tearDown(self):
        helpers.cleanup(self.store, self.path)

    def test_matches_segments_and_links_all_carry(self):
        stats = ingest_mod.carry_matches(self.store, [self.a, self.b], self.merged)
        self.assertEqual(stats["matches"], 3)
        self.assertEqual(stats["segments"], 1)
        self.assertEqual(stats["shared"], 1)
        self.assertEqual(len(self.store.matches(self.merged)), 3)

    def test_platform_is_preserved(self):
        ingest_mod.carry_matches(self.store, [self.a, self.b], self.merged)
        sources = {m["source"] for m in self.store.matches(self.merged)}
        self.assertEqual(
            sources, {"23andMe", "Ancestry"},
            "side inference compares lists within one platform, so the "
            "platform must survive the merge",
        )

    def test_side_assignments_survive(self):
        ingest_mod.carry_matches(self.store, [self.a, self.b], self.merged)
        one = [m for m in self.store.matches(self.merged) if m["name"] == "Cousin One"]
        self.assertEqual(one[0]["side"], "paternal")

    def test_carrying_twice_does_not_duplicate(self):
        ingest_mod.carry_matches(self.store, [self.a, self.b], self.merged)
        ingest_mod.carry_matches(self.store, [self.a, self.b], self.merged)
        self.assertEqual(len(self.store.matches(self.merged)), 3)
        segs = [s for m in self.store.matches(self.merged)
                for s in self.store.match_segments(m["id"])]
        self.assertEqual(len(segs), 1, "segments must not accumulate on re-run")


if __name__ == "__main__":
    unittest.main(verbosity=2)

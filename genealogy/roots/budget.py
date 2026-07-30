"""Spending limits: what to buy next, and what it will cost.

The research planner in `roots.evidence` produces more work than anyone wants
to pay for in one go. This module turns that list into a plan that fits an
actual budget, and keeps a ledger of what has been spent so far.

Two modelling points matter more than the arithmetic.

**Subscription costs are shared, not per-record.** A census lookup and a
parish register lookup and forty more all cost one month's subscription
between them. Charging each task separately -- which is what a naive
allocation does -- makes subscription research look ruinous and pushes you
toward buying certificates you did not need. The correct treatment is that
the first subscription task costs a month's fee and every subsequent one is
free, which is also the practical advice: decide what you need from a
subscription site, then take one month and do all of it.

**Free work should never be deferred.** The GRO birth index carries the
mother's maiden surname and costs nothing to search; the surviving Irish
censuses are free; FamilySearch is free. Any plan that spends money before
exhausting these is wrong regardless of budget, so free tasks are always
included and never counted against the limit.

On the prices below: `roots.evidence` deliberately keeps costs as bands
rather than figures, on the grounds that a stale number in code is worse than
no number. Fitting a budget needs actual money, so the table here is the one
place figures live. It is dated, sourced to docs/RECORD_SOURCES.md, and
overridable from a file -- and if it is more than a year old you should treat
its output as indicative and recheck.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .store import Store

#: When the figures below were last checked against providers.
#: See docs/RECORD_SOURCES.md for sourcing and caveats.
VERIFIED_ON = "2026-07-30"

#: Cost of one record, in pounds, keyed by repository.
PER_RECORD: Dict[str, float] = {
    "GRO": 12.50,               # full certificate; the safe assumption
    "GRO / FreeBMD": 12.50,
    "ScotlandsPeople": 1.50,    # 6 credits at 25p
    "GRONI": 2.50,
    "irishgenealogy.ie / GRONI": 0.0,
    "gov.uk Find a Will": 1.50,
    "National Archives of Ireland": 0.0,
    "GRO / FreeBMD / FamilySearch": 12.50,
}

#: One month of access, in pounds, for sites where the marginal record is free.
SUBSCRIPTION_MONTH: Dict[str, float] = {
    "FindMyPast / Ancestry / FamilySearch": 9.99,
    "FamilySearch / FindMyPast / county record office": 9.99,
    "Ancestry / FindMyPast / FamilySearch trees": 10.99,
    "default": 9.99,
}

#: Cheaper alternatives that exist for some GRO records, applied when the
#: task text shows the record is old enough to qualify.
GRO_DIGITAL_IMAGE = 3.00
GRO_PDF = 8.00


@dataclass
class CostEstimate:
    amount: float
    kind: str  # 'free' | 'per-item' | 'subscription'
    note: str = ""


def estimate_cost(task: Any, costs: Optional[Dict[str, Any]] = None) -> CostEstimate:
    """What one research task would cost to carry out."""
    table = costs or {}
    per_record = {**PER_RECORD, **table.get("per_record", {})}
    subs = {**SUBSCRIPTION_MONTH, **table.get("subscription_month", {})}

    band = getattr(task, "cost_band", "unknown")
    repo = getattr(task, "repository", "") or ""
    record_set = (getattr(task, "record_set", "") or "").lower()

    if band == "free":
        return CostEstimate(0.0, "free", "no charge")

    if band == "subscription":
        return CostEstimate(
            subs.get(repo, subs["default"]), "subscription",
            "one month's access, shared with every other subscription task",
        )

    if band == "low":
        amount = per_record.get(repo)
        if amount is None:
            amount = 5.00
        # GRO offers much cheaper formats for older records, and the planner
        # says which era it is aiming at.
        if repo.startswith("GRO") and "birth" in record_set:
            year = _year_in(record_set)
            if year and year <= 1934:
                return CostEstimate(
                    GRO_DIGITAL_IMAGE, "per-item",
                    "GRO digital image, available for older births",
                )
            if year and year >= 1984:
                return CostEstimate(GRO_PDF, "per-item", "GRO PDF")
        if repo.startswith("GRO") and "marriage" in record_set:
            return CostEstimate(
                per_record.get("GRO", 12.50), "per-item",
                "marriages have no cheaper GRO format at any date",
            )
        return CostEstimate(amount, "per-item", "")

    return CostEstimate(0.0, "free", "cost unknown; treated as free to consider")


def _year_in(text: str) -> Optional[int]:
    import re

    hits = re.findall(r"\b(1[6-9]\d{2}|20[0-2]\d)\b", text)
    return int(hits[0]) if hits else None


@dataclass
class BudgetPlan:
    budget: float
    selected: List[Tuple[Any, CostEstimate]] = field(default_factory=list)
    deferred: List[Tuple[Any, CostEstimate]] = field(default_factory=list)
    subscription_taken: Optional[Tuple[str, float]] = None
    spent: float = 0.0

    @property
    def free_count(self) -> int:
        return sum(1 for _t, c in self.selected if c.kind == "free")

    @property
    def remaining(self) -> float:
        return max(0.0, self.budget - self.spent)

    def summary(self) -> str:
        paid = len(self.selected) - self.free_count
        bits = [
            f"{len(self.selected)} tasks selected ({self.free_count} free, {paid} paid)",
            f"£{self.spent:.2f} of £{self.budget:.2f}",
        ]
        if self.deferred:
            bits.append(f"{len(self.deferred)} deferred")
        return ", ".join(bits)


def plan_within_budget(
    tasks: Sequence[Any],
    budget: float,
    costs: Optional[Dict[str, Any]] = None,
) -> BudgetPlan:
    """Choose the highest-yield set of tasks that fits a budget.

    Greedy by priority per pound, which is the standard knapsack heuristic
    and is more than accurate enough here given that the priorities are
    themselves estimates. Two refinements matter:

    Free tasks are taken unconditionally, since deferring them cannot save
    money and can only delay the answer.

    The subscription decision is evaluated as a whole rather than task by
    task. Two plans are costed -- one that takes a month's subscription and
    one that does not -- and the better is returned. That is the only way to
    get the right answer when a single fee unlocks many tasks.
    """
    plan_without = _greedy(tasks, budget, costs, take_subscription=False)
    plan_with = _greedy(tasks, budget, costs, take_subscription=True)
    best = max(
        (plan_without, plan_with),
        key=lambda p: (_yield(p), -p.spent),
    )
    best.budget = budget
    return best


def _yield(plan: BudgetPlan) -> float:
    return sum(getattr(t, "priority", 0.0) for t, _c in plan.selected)


def _greedy(
    tasks: Sequence[Any],
    budget: float,
    costs: Optional[Dict[str, Any]],
    take_subscription: bool,
) -> BudgetPlan:
    plan = BudgetPlan(budget=budget)
    priced = [(t, estimate_cost(t, costs)) for t in tasks]

    free = [(t, c) for t, c in priced if c.kind == "free"]
    subscription = [(t, c) for t, c in priced if c.kind == "subscription"]
    per_item = [(t, c) for t, c in priced if c.kind == "per-item"]

    plan.selected.extend(free)

    remaining = budget
    if take_subscription and subscription:
        repo, fee = _cheapest_subscription(subscription)
        if fee <= remaining:
            plan.subscription_taken = (repo, fee)
            plan.spent += fee
            remaining -= fee
            plan.selected.extend(subscription)
        else:
            plan.deferred.extend(subscription)
    else:
        plan.deferred.extend(subscription)

    per_item.sort(
        key=lambda tc: -(getattr(tc[0], "priority", 0.0) / max(tc[1].amount, 0.01))
    )
    for task, cost in per_item:
        if cost.amount <= remaining:
            plan.selected.append((task, cost))
            plan.spent += cost.amount
            remaining -= cost.amount
        else:
            plan.deferred.append((task, cost))
    return plan


def _cheapest_subscription(
    subscription: Sequence[Tuple[Any, CostEstimate]]
) -> Tuple[str, float]:
    best_repo, best_fee = "a subscription site", 1e9
    for task, cost in subscription:
        if cost.amount < best_fee:
            best_fee = cost.amount
            best_repo = getattr(task, "repository", "a subscription site")
    return best_repo, best_fee


# ---------------------------------------------------------------------------
# ledger


@dataclass
class Ledger:
    cap: Optional[float]
    spent: float
    records: int

    @property
    def remaining(self) -> Optional[float]:
        return None if self.cap is None else max(0.0, self.cap - self.spent)

    def summary(self) -> str:
        if self.cap is None:
            return f"£{self.spent:.2f} spent across {self.records} records (no cap set)"
        return (
            f"£{self.spent:.2f} of £{self.cap:.2f} spent across {self.records} "
            f"records; £{self.remaining:.2f} left"
        )


def ledger(store: Store) -> Ledger:
    """What has actually been spent, from the recorded sources."""
    row = store.db.execute(
        "SELECT COUNT(*) AS n, COALESCE(SUM(cost), 0) AS c FROM source"
    ).fetchone()
    cap = store.get_meta("budget_cap")
    return Ledger(
        cap=float(cap) if cap else None,
        spent=float(row["c"] or 0.0),
        records=int(row["n"]),
    )


def set_cap(store: Store, amount: Optional[float]) -> None:
    store.set_meta("budget_cap", "" if amount is None else str(amount))


def load_costs(path: str) -> Dict[str, Any]:
    """Load a cost override file, for when the built-in table has aged."""
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)

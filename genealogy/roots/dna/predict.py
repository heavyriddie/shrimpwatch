"""Rank candidate relationships for an observed amount of shared DNA.

The question a genealogist actually asks is "this person shares 212 cM with
me -- who could they be?", and the honest answer is a ranked list with
probabilities, not a single label.  Shared-cM ranges for adjacent
relationships overlap heavily; at 212 cM a second cousin, a first cousin
twice removed, and a half first cousin once removed are all entirely
plausible.

Three ingredients go into the ranking:

*Likelihood* -- how probable is this much sharing under each candidate
relationship.  Comes from the pedigree simulator, including the possibility
of no detectable sharing at all, which matters enormously for distant
relationships.

*Prior* -- how many relatives of that type a person actually has.  This is
the ingredient most tools omit, and it dominates the answer at low cM.  You
have roughly 7 first cousins but roughly 900 fourth cousins, so a 40 cM
match is far more likely to be a distant relative than an unusually low
close one.  The prior comes from a simple demographic model (average
children per couple, rate of half-sibships) that you can tune or switch off.

*Optional evidence* -- segment count, fully-identical sharing, and age
difference.  Fully-identical regions are decisive: they only appear at full
sibling or closer, so their presence or absence separates relationships that
total sharing alone cannot.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from ..genome import GeneticMap
from .simulate import CATALOGUE, Distribution, Relationship, cached_distribution

GENERATION_YEARS = 30.0


# ---------------------------------------------------------------------------
# priors


def expected_relative_count(
    rel: Relationship, kids: float = 2.5, half_rate: float = 0.25
) -> float:
    """Roughly how many relatives of this type a person has.

    ``kids`` is the average number of children per couple across the
    generations in question; ``half_rate`` is the average number of
    half-sibships an ancestor has by other partners.  The model is crude but
    its shape is right, and its shape is what matters: relative counts grow
    by roughly 5x per cousin degree, which is the force that pulls low-cM
    matches toward distant explanations.
    """
    if rel.special == "twin":
        return 0.004  # identical twinning rate, roughly
    if rel.special == "double":
        return 0.5

    def one_direction(up: int, down: int) -> float:
        if up == 0:
            # Direct line: ancestors above, descendants below.
            return 2.0 ** down + kids ** down
        couples = 2.0 ** (up - 1)
        if rel.ancestors == 1:
            # Half relationships hang off an ancestor's other partnerships.
            siblings = 2.0 * couples * half_rate
        else:
            siblings = couples * max(kids - 1.0, 0.0)
        return siblings * (kids ** (down - 1))

    total = one_direction(rel.up, rel.down)
    if rel.up != rel.down and rel.up > 0:
        total += one_direction(rel.down, rel.up)
    return max(total, 1e-6)


def coexistence_prior(rel: Relationship) -> float:
    """Penalise relationships that need an implausible age difference.

    Both people have to be alive, and old enough to test, at the same moment.
    A relationship whose two sides sit g generations apart implies an age
    difference near 30g years, so a "2x-great-grandparent" match implies
    about 120 years and effectively cannot appear in a database however many
    such ancestors you have on paper.

    This is deliberately gentle -- a 30-year gap is completely ordinary and a
    60-year gap is common -- and it is skipped entirely when a real age
    difference is known, since the observed value is better evidence than
    this stand-in for it.
    """
    gap = abs(rel.generation_gap)
    if gap == 0:
        return 1.0
    implied_years = 30.0 * gap
    return math.exp(-0.5 * (implied_years / 40.0) ** 2)


def age_prior(rel: Relationship, age_gap_years: Optional[float]) -> float:
    """Weight a relationship by how well it fits an observed age difference.

    ``age_gap_years`` is the match's birth year minus yours, so a positive
    number means the match is younger.  A relationship implying a generation
    gap of g predicts an age difference near g * 30 years; a first cousin who
    is 62 years younger than you is far more likely a first cousin twice
    removed.
    """
    if age_gap_years is None:
        return 1.0
    expected = rel.generation_gap * GENERATION_YEARS
    sd = 18.0  # generation lengths vary a lot; keep this forgiving
    z = (abs(age_gap_years) - abs(expected)) / sd
    return math.exp(-0.5 * z * z)


# ---------------------------------------------------------------------------
# prediction


@dataclass
class Candidate:
    rel: Relationship
    likelihood: float
    prior: float
    posterior: float = 0.0
    dist: Optional[Distribution] = None

    @property
    def name(self) -> str:
        return self.rel.name

    def range_text(self) -> str:
        if not self.dist:
            return ""
        lo, hi = self.dist.range(0.05, 0.95)
        return f"{lo:.0f}-{hi:.0f} cM (90% of simulated pairs)"


@dataclass
class Prediction:
    observed_cm: float
    candidates: List[Candidate] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def top(self, n: int = 8) -> List[Candidate]:
        return self.candidates[:n]

    def grouped_by_degree(self) -> Dict[int, float]:
        """Collapse to total posterior per meiotic degree of separation."""
        out: Dict[int, float] = {}
        for c in self.candidates:
            out[c.rel.meioses] = out.get(c.rel.meioses, 0.0) + c.posterior
        return out

    def credible_set(self, mass: float = 0.90) -> List[Candidate]:
        """Smallest set of relationships covering ``mass`` of the posterior."""
        acc, out = 0.0, []
        for c in self.candidates:
            out.append(c)
            acc += c.posterior
            if acc >= mass:
                break
        return out


def predict(
    store,
    observed_cm: float,
    *,
    gmap: Optional[GeneticMap] = None,
    min_cm: float = 7.0,
    iterations: int = 1500,
    segments: Optional[int] = None,
    ibd2_cm: Optional[float] = None,
    age_gap_years: Optional[float] = None,
    use_prior: bool = True,
    kids: float = 2.5,
    keys: Optional[Sequence[str]] = None,
    calibration: float = 1.0,
) -> Prediction:
    """Rank relationships that could explain ``observed_cm`` of shared DNA."""
    pred = Prediction(observed_cm=observed_cm)
    rels = [CATALOGUE[k] for k in keys] if keys else list(CATALOGUE.values())

    raw: List[Candidate] = []
    for rel in rels:
        dist = cached_distribution(
            store, rel, iterations=iterations, min_cm=min_cm, gmap=gmap,
            calibration=calibration,
        )
        lik = _likelihood(dist, observed_cm, segments, ibd2_cm)
        if lik <= 0:
            continue
        prior = 1.0
        if use_prior:
            prior = expected_relative_count(rel, kids=kids)
            if age_gap_years is None:
                prior *= coexistence_prior(rel)
        prior *= age_prior(rel, age_gap_years)
        raw.append(Candidate(rel=rel, likelihood=lik, prior=prior, dist=dist))

    total = sum(c.likelihood * c.prior for c in raw)
    if total > 0:
        for c in raw:
            c.posterior = c.likelihood * c.prior / total
    raw.sort(key=lambda c: c.posterior, reverse=True)
    pred.candidates = raw
    pred.notes = _notes(observed_cm, raw, ibd2_cm)
    return pred


def _likelihood(
    dist: Distribution,
    observed_cm: float,
    segments: Optional[int],
    ibd2_cm: Optional[float],
) -> float:
    if observed_cm <= 0:
        return dist.p_undetected
    # Mixture: the density over nonzero outcomes, scaled by the chance of
    # producing any detectable sharing at all.
    lik = dist.density(observed_cm) * (1.0 - dist.p_undetected)
    if lik <= 0:
        return 0.0

    if segments is not None and dist.segment_counts:
        # Empirical probability of that many segments, Laplace-smoothed so a
        # count we happened not to simulate does not zero out the candidate.
        n = len(dist.segment_counts)
        hits = sum(1 for s in dist.segment_counts if s == segments)
        lik *= (hits + 0.5) / (n + 0.5 * 60)

    if ibd2_cm is not None and dist.ibd2:
        # Fully-identical sharing is a hard discriminator: only full siblings
        # and closer produce it in quantity.
        n = len(dist.ibd2)
        typical = sum(dist.ibd2) / n
        if ibd2_cm > 100:
            frac = sum(1 for v in dist.ibd2 if v > 100) / n
            lik *= max(frac, 1e-4)
        elif typical > 100:
            frac = sum(1 for v in dist.ibd2 if v <= 100) / n
            lik *= max(frac, 1e-4)
    return lik


def _notes(observed_cm: float, cands: Sequence[Candidate], ibd2_cm: Optional[float]) -> List[str]:
    notes: List[str] = []
    if not cands:
        notes.append("no candidate relationship can explain this amount of sharing")
        return notes
    if observed_cm > 3300:
        notes.append(
            "sharing this high means parent/child or identical twin; nothing "
            "else reaches it"
        )
    elif 1300 <= observed_cm <= 2300 and ibd2_cm is None:
        notes.append(
            "this range is the classic ambiguity: half sibling, grandparent, "
            "and aunt/uncle all sit here. Total sharing cannot separate them -- "
            "use ages, the X chromosome, and which of them appears on your "
            "mother's match list"
        )
    if observed_cm and observed_cm < 30:
        notes.append(
            "below about 30 cM a substantial share of reported matches are not "
            "genealogically meaningful within a documentable number of "
            "generations, even when the segment is real"
        )
    top = cands[0]
    if len(cands) > 1 and cands[1].posterior > top.posterior * 0.6:
        notes.append(
            "the top candidates are close in probability; treat this as a set "
            "of possibilities rather than an identification"
        )
    return notes


# ---------------------------------------------------------------------------
# comparing an observation against a specific hypothesis


@dataclass
class FitResult:
    rel: Relationship
    observed_cm: float
    expected_cm: float
    percentile: float
    plausible: bool
    comment: str


def fit_hypothesis(
    store,
    rel: Relationship,
    observed_cm: float,
    *,
    gmap: Optional[GeneticMap] = None,
    min_cm: float = 7.0,
    iterations: int = 1500,
    calibration: float = 1.0,
) -> FitResult:
    """Test a specific proposed relationship against observed sharing.

    This is the workhorse for checking a paper trail against DNA: the tree
    says these two are third cousins, so is 340 cM believable?  (It is not --
    that is roughly the 99.9th percentile for third cousins, and points at
    either an error in the tree or an additional relationship line.)
    """
    dist = cached_distribution(
        store, rel, iterations=iterations, min_cm=min_cm, gmap=gmap, calibration=calibration
    )
    n = dist.n or 1
    below = sum(1 for s in dist.samples if s < observed_cm)
    pct = below / n
    lo, hi = dist.range(0.02, 0.98)
    plausible = lo <= observed_cm <= hi
    if observed_cm > hi:
        comment = (
            f"observed sharing is higher than {pct:.1%} of simulated "
            f"{rel.name} pairs -- the pair is likely more closely related than "
            "the tree says, or related through more than one line"
        )
    elif observed_cm < lo:
        comment = (
            f"observed sharing is lower than {1 - pct:.1%} of simulated "
            f"{rel.name} pairs -- possible for real relatives at this distance, "
            "but check the connection"
        )
    else:
        comment = f"consistent with {rel.name} (at the {pct:.0%} percentile)"
    return FitResult(
        rel=rel,
        observed_cm=observed_cm,
        expected_cm=dist.mean,
        percentile=pct,
        plausible=plausible,
        comment=comment,
    )

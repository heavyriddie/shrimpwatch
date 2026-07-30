"""Monte Carlo meiosis simulation over an explicit pedigree.

Everything this project claims about "how much DNA should a second cousin
once removed share" is computed here rather than looked up in a table.  The
model is deliberately concrete: build the actual pedigree, give each founder
two uniquely-labelled chromosome copies, simulate recombination down every
line, and see which stretches of chromosome the two target people ended up
carrying from the same founder haplotype.

Why simulate the whole pedigree instead of using a closed-form expectation?
Because the shortcut of "treat the two common ancestors as independent" is
wrong in a way that matters.  Consider an aunt and her nephew.  Both of the
nephew's shared grandparents reach him through a single parent -- his
father -- so at any given position his paternal chromosome traces to exactly
one of them; his maternal chromosome is unrelated to his aunt entirely.
Treating the two grandparents as independent contributors predicts about
1525 cM.  Simulating the pedigree gives 1742 cM, which is what aunt-nephew
pairs actually show.  The same subtlety affects every relationship that
descends through a shared couple.

Calibration against known averages, at a 7 cM detection threshold:

    parent/child        3545 cM   (whole genome, by construction)
    full siblings       ~2600 cM  (published average 2613)
    half sib / grandparent / aunt-uncle
                        ~1740 cM  (published average ~1750)
    first cousins       ~860 cM   (published average 866)
    second cousins      ~215 cM   (published average 229)

The close relationships land essentially on the published averages.  Cousin
levels come in a few percent low, and the gap widens with distance -- which
is expected, since crowd-sourced averages accumulate small false-positive
segments exactly where true sharing thins out.  `calibration` scales
simulated totals if you would rather match the published tables directly.

Recombination is modelled as a Poisson process at one crossover per Morgan
with no interference.  Interference makes crossovers more evenly spaced than
Poisson, which narrows the spread slightly without moving the mean; ignoring
it is the standard simplification and is conservative for our purposes.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..genome import AUTOSOMES, CHROM_CM, GeneticMap

Interval = Tuple[float, float, int]  # start cM, end cM, founder-haplotype label
Span = Tuple[float, float]

SIM_VERSION = 2  # bump to invalidate cached distributions


# ---------------------------------------------------------------------------
# relationship catalogue


@dataclass(frozen=True)
class Relationship:
    """A relationship expressed as a pedigree shape.

    ``up`` is the number of meioses from person A to the most recent common
    ancestor, ``down`` the number from person B.  ``ancestors`` is 2 when the
    pair descends from a couple (a "full" relationship) and 1 when they
    descend from a single shared ancestor with different partners (a "half"
    relationship).  ``up == 0`` means A *is* the common ancestor, i.e. a
    direct line.
    """

    key: str
    name: str
    up: int
    down: int
    ancestors: int = 2
    special: Optional[str] = None

    @property
    def meioses(self) -> int:
        return self.up + self.down

    @property
    def generation_gap(self) -> int:
        """Positive when B sits generations below A."""
        return self.down - self.up

    def label(self) -> str:
        return self.name


def _cousin_name(n: int, removed: int, half: bool) -> str:
    ord_map = {1: "first", 2: "second", 3: "third", 4: "fourth", 5: "fifth",
               6: "sixth", 7: "seventh", 8: "eighth"}
    base = f"{ord_map.get(n, str(n) + 'th')} cousin"
    if half:
        base = "half " + base
    if removed == 1:
        base += " once removed"
    elif removed == 2:
        base += " twice removed"
    elif removed > 2:
        base += f" {removed} times removed"
    return base


def catalogue(max_cousin: int = 6) -> List[Relationship]:
    """The relationship space we score candidates against."""
    rels: List[Relationship] = [
        Relationship("twin", "identical twin", 0, 0, 1, special="twin"),
        Relationship("parent", "parent/child", 0, 1, 1),
        Relationship("gparent", "grandparent/grandchild", 0, 2, 1),
        Relationship("ggparent", "great-grandparent/great-grandchild", 0, 3, 1),
        Relationship("gggparent", "2x-great-grandparent/grandchild", 0, 4, 1),
        Relationship("sibling", "full sibling", 1, 1, 2),
        Relationship("half_sibling", "half sibling", 1, 1, 1),
        Relationship("avuncular", "aunt/uncle - niece/nephew", 1, 2, 2),
        Relationship("half_avuncular", "half aunt/uncle - niece/nephew", 1, 2, 1),
        Relationship("gavuncular", "great-aunt/uncle - great-niece/nephew", 1, 3, 2),
        Relationship("half_gavuncular", "half great-aunt/uncle", 1, 3, 1),
        Relationship("ggavuncular", "2x-great-aunt/uncle", 1, 4, 2),
        Relationship("double_1c", "double first cousin", 2, 2, 4, special="double"),
    ]
    for n in range(1, max_cousin + 1):
        for removed in range(0, 4):
            up = n + 1
            down = n + 1 + removed
            rels.append(
                Relationship(
                    f"{n}c{removed}r" if removed else f"{n}c",
                    _cousin_name(n, removed, half=False),
                    up, down, 2,
                )
            )
            rels.append(
                Relationship(
                    f"h{n}c{removed}r" if removed else f"h{n}c",
                    _cousin_name(n, removed, half=True),
                    up, down, 1,
                )
            )
    return rels


CATALOGUE = {r.key: r for r in catalogue()}


def _greats(n: int) -> str:
    """The 'great-' prefix genealogists use, abbreviated past two."""
    if n <= 0:
        return ""
    if n == 1:
        return "great-"
    if n == 2:
        return "great-great-"
    return f"{n}x-great-"


def describe_shape(up: int, down: int, ancestors: int) -> str:
    """Name any pedigree shape, in the terms a genealogist would use.

    Three families of relationship, split by how far person A sits from the
    common ancestor.  A direct line (``up == 0``) gives parents and
    grandparents; one step off the line (``up == 1``) gives siblings, aunts
    and great-aunts; two or more steps gives cousins, whose degree is
    ``up - 1`` and whose removal is the difference between the two sides.
    """
    half = "half " if ancestors == 1 else ""
    if up == 0:
        if down == 1:
            return "parent/child"
        greats = _greats(down - 2)
        return f"{greats}grandparent/{greats}grandchild"
    if up == 1:
        if down == 1:
            return f"{half}sibling" if half else "full sibling"
        greats = _greats(down - 2)
        return f"{half}{greats}aunt/uncle - {greats}niece/nephew"
    return _cousin_name(up - 1, down - up, half=(ancestors == 1))


def relationship_for(up: int, down: int, ancestors: int) -> Relationship:
    """Name an arbitrary pedigree shape, whether or not it is in the catalogue."""
    if up > down:
        up, down = down, up
    # A direct-line relationship only ever runs through one ancestor: your
    # grandmother is your grandmother regardless of whom she married.
    if up == 0:
        ancestors = 1
    for rel in CATALOGUE.values():
        if rel.up == up and rel.down == down and rel.ancestors == ancestors and not rel.special:
            return rel
    prefix = "h" if ancestors == 1 else ""
    key = f"anc{down}" if up == 0 else f"{prefix}{up}_{down}"
    return Relationship(key, describe_shape(up, down, ancestors), up, down, ancestors)


# ---------------------------------------------------------------------------
# interval algebra


def _slice(intervals: Sequence[Interval], a: float, b: float) -> List[Interval]:
    out: List[Interval] = []
    for s, e, lab in intervals:
        if e <= a:
            continue
        if s >= b:
            break
        out.append((max(s, a), min(e, b), lab))
    return out


def _coalesce(intervals: Sequence[Interval]) -> List[Interval]:
    out: List[Interval] = []
    for s, e, lab in intervals:
        if e - s <= 0:
            continue
        if out and out[-1][2] == lab and abs(out[-1][1] - s) < 1e-9:
            out[-1] = (out[-1][0], e, lab)
        else:
            out.append((s, e, lab))
    return out


def _matching(a: Sequence[Interval], b: Sequence[Interval]) -> List[Span]:
    """Spans where two haplotypes carry the same founder label."""
    out: List[Span] = []
    i = j = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if s < e and a[i][2] == b[j][2]:
            if out and abs(out[-1][1] - s) < 1e-9:
                out[-1] = (out[-1][0], e)
            else:
                out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def _union(spans: Iterable[Span]) -> List[Span]:
    ordered = sorted(spans)
    out: List[Span] = []
    for s, e in ordered:
        if out and s <= out[-1][1] + 1e-9:
            if e > out[-1][1]:
                out[-1] = (out[-1][0], e)
        else:
            out.append((s, e))
    return out


def _intersect(a: Sequence[Span], b: Sequence[Span]) -> List[Span]:
    out: List[Span] = []
    i = j = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def _total(spans: Iterable[Span], min_cm: float = 0.0) -> float:
    return sum(e - s for s, e in spans if (e - s) >= min_cm)


# ---------------------------------------------------------------------------
# meiosis


class _Person:
    __slots__ = ("h1", "h2")

    def __init__(self, h1: List[Interval], h2: List[Interval]):
        self.h1 = h1
        self.h2 = h2


class _Pedigree:
    """Builds people on one chromosome, handing out fresh founder labels."""

    def __init__(self, length: float, rng: random.Random):
        self.length = length
        self.rng = rng
        self._next_label = 0

    def founder(self) -> _Person:
        a, b = self._next_label, self._next_label + 1
        self._next_label += 2
        return _Person([(0.0, self.length, a)], [(0.0, self.length, b)])

    def gamete(self, p: _Person) -> List[Interval]:
        """One recombined chromosome, ready to be passed to a child."""
        rng = self.rng
        length = self.length
        points: List[float] = []
        # Crossovers arrive as a Poisson process at 1 per Morgan (100 cM),
        # so the gaps between them are exponential with mean 100 cM.
        x = rng.expovariate(0.01)
        while x < length:
            points.append(x)
            x += rng.expovariate(0.01)
        bounds = [0.0] + points + [length]
        cur = rng.getrandbits(1)
        out: List[Interval] = []
        for i in range(len(bounds) - 1):
            src = p.h1 if cur == 0 else p.h2
            out.extend(_slice(src, bounds[i], bounds[i + 1]))
            cur ^= 1
        return _coalesce(out)

    def child(self, p1: _Person, p2: _Person) -> _Person:
        return _Person(self.gamete(p1), self.gamete(p2))

    def descend(self, p: _Person, generations: int) -> _Person:
        """Walk down `generations` meioses, marrying in unrelated founders."""
        for _ in range(generations):
            p = self.child(p, self.founder())
        return p


def _build_pair(ped: _Pedigree, rel: Relationship) -> Tuple[_Person, _Person]:
    if rel.special == "twin":
        p = ped.founder()
        return p, _Person(list(p.h1), list(p.h2))

    if rel.special == "double":
        # Two siblings marry two siblings; their children share all four
        # grandparents.  Genetically this lands close to a half sibling.
        f1, f2 = ped.founder(), ped.founder()
        g1, g2 = ped.founder(), ped.founder()
        c1, c2 = ped.child(f1, f2), ped.child(f1, f2)
        d1, d2 = ped.child(g1, g2), ped.child(g1, g2)
        return ped.child(c1, d1), ped.child(c2, d2)

    if rel.up == 0:
        # Direct line: A is the ancestor, B descends from A through partners
        # unrelated to the line.
        anc = ped.founder()
        head = ped.child(anc, ped.founder())
        return anc, ped.descend(head, rel.down - 1)

    if rel.ancestors == 2:
        f1, f2 = ped.founder(), ped.founder()
        head_a = ped.child(f1, f2)
        head_b = ped.child(f1, f2)
    else:
        shared = ped.founder()
        head_a = ped.child(shared, ped.founder())
        head_b = ped.child(shared, ped.founder())
    return ped.descend(head_a, rel.up - 1), ped.descend(head_b, rel.down - 1)


@dataclass
class Draw:
    """One simulated pair of people."""

    total_cm: float
    segments: int
    largest_cm: float
    ibd2_cm: float


def simulate_once(
    rel: Relationship,
    rng: random.Random,
    min_cm: float = 7.0,
    chrom_lengths: Optional[Dict[str, float]] = None,
) -> Draw:
    lengths = chrom_lengths or {c: CHROM_CM[c] for c in AUTOSOMES}
    total = 0.0
    n_seg = 0
    largest = 0.0
    ibd2 = 0.0
    for _chrom, length in lengths.items():
        ped = _Pedigree(length, rng)
        a, b = _build_pair(ped, rel)
        m11 = _matching(a.h1, b.h1)
        m12 = _matching(a.h1, b.h2)
        m21 = _matching(a.h2, b.h1)
        m22 = _matching(a.h2, b.h2)
        hir = _union(m11 + m12 + m21 + m22)
        kept = [(s, e) for s, e in hir if (e - s) >= min_cm]
        for s, e in kept:
            total += e - s
            n_seg += 1
            if e - s > largest:
                largest = e - s
        # Both chromosome copies matching: only happens for siblings and closer.
        both = _union(_intersect(m11, m22) + _intersect(m12, m21))
        ibd2 += _total(both, min_cm)
    return Draw(total_cm=total, segments=n_seg, largest_cm=largest, ibd2_cm=ibd2)


# ---------------------------------------------------------------------------
# distributions


@dataclass
class Distribution:
    key: str
    name: str
    samples: List[float] = field(default_factory=list)
    segment_counts: List[int] = field(default_factory=list)
    largest: List[float] = field(default_factory=list)
    ibd2: List[float] = field(default_factory=list)
    min_cm: float = 7.0

    @property
    def n(self) -> int:
        return len(self.samples)

    @property
    def mean(self) -> float:
        return sum(self.samples) / self.n if self.n else 0.0

    @property
    def stdev(self) -> float:
        if self.n < 2:
            return 0.0
        m = self.mean
        return math.sqrt(sum((x - m) ** 2 for x in self.samples) / (self.n - 1))

    @property
    def p_undetected(self) -> float:
        """Chance the pair shares nothing above the detection threshold.

        This is the number people forget: a real fourth cousin has a
        substantial chance of showing no shared DNA at all, so absence of a
        match is weak evidence against a distant relationship.
        """
        if not self.n:
            return 0.0
        return sum(1 for x in self.samples if x <= 0) / self.n

    def percentile(self, p: float) -> float:
        if not self.samples:
            return 0.0
        ordered = sorted(self.samples)
        k = (len(ordered) - 1) * p
        lo, hi = math.floor(k), math.ceil(k)
        if lo == hi:
            return ordered[int(k)]
        return ordered[lo] + (ordered[hi] - ordered[lo]) * (k - lo)

    def range(self, lo: float = 0.05, hi: float = 0.95) -> Tuple[float, float]:
        return self.percentile(lo), self.percentile(hi)

    def density(self, x: float) -> float:
        """Smoothed likelihood of observing ``x`` cM under this relationship.

        Gaussian kernel density with Silverman's bandwidth.  For
        relationships that frequently produce no detectable sharing, the
        point mass at zero is handled separately by the caller.
        """
        nz = [s for s in self.samples if s > 0]
        if not nz:
            return 0.0
        n = len(nz)
        m = sum(nz) / n
        var = sum((v - m) ** 2 for v in nz) / n if n > 1 else 1.0
        sd = math.sqrt(var) or 1.0
        h = 1.06 * sd * (n ** -0.2)
        # Parent/child and identical twins have no simulated spread at all --
        # they share the whole genome every time. Real reported totals for
        # them still vary by a percent or two through call rates, map
        # differences and segment-edge handling, so the bandwidth gets a
        # floor proportional to the mean. Without it those relationships
        # would have zero likelihood at any value but exactly 3545 cM.
        h = max(h, 0.02 * m, 2.0)
        acc = 0.0
        inv = 1.0 / h
        for v in nz:
            z = (x - v) * inv
            if -6.0 < z < 6.0:
                acc += math.exp(-0.5 * z * z)
        return acc * inv / (n * math.sqrt(2 * math.pi))

    def to_json(self) -> Dict:
        return {
            "key": self.key,
            "name": self.name,
            "min_cm": self.min_cm,
            "samples": [round(s, 2) for s in self.samples],
            "segment_counts": self.segment_counts,
            "largest": [round(s, 2) for s in self.largest],
            "ibd2": [round(s, 2) for s in self.ibd2],
        }

    @classmethod
    def from_json(cls, blob: Dict) -> "Distribution":
        return cls(
            key=blob["key"],
            name=blob["name"],
            samples=[float(s) for s in blob["samples"]],
            segment_counts=[int(s) for s in blob.get("segment_counts", [])],
            largest=[float(s) for s in blob.get("largest", [])],
            ibd2=[float(s) for s in blob.get("ibd2", [])],
            min_cm=float(blob.get("min_cm", 7.0)),
        )


def simulate(
    rel: Relationship,
    iterations: int = 2000,
    min_cm: float = 7.0,
    seed: Optional[int] = None,
    gmap: Optional[GeneticMap] = None,
    calibration: float = 1.0,
) -> Distribution:
    """Run the simulation and collect a distribution of outcomes."""
    rng = random.Random(seed if seed is not None else 0xC0FFEE + hash(rel.key) % 100003)
    lengths = (
        {c: gmap.chrom_cm(c) for c in AUTOSOMES if gmap.chrom_cm(c) > 0}
        if gmap
        else {c: CHROM_CM[c] for c in AUTOSOMES}
    )
    dist = Distribution(key=rel.key, name=rel.name, min_cm=min_cm)
    for _ in range(iterations):
        draw = simulate_once(rel, rng, min_cm=min_cm, chrom_lengths=lengths)
        dist.samples.append(draw.total_cm * calibration)
        dist.segment_counts.append(draw.segments)
        dist.largest.append(draw.largest_cm * calibration)
        dist.ibd2.append(draw.ibd2_cm * calibration)
    return dist


def cached_distribution(
    store,
    rel: Relationship,
    iterations: int = 2000,
    min_cm: float = 7.0,
    gmap: Optional[GeneticMap] = None,
    calibration: float = 1.0,
) -> Distribution:
    """Fetch a distribution from the project cache, simulating on a miss."""
    src = gmap.source if gmap else "builtin"
    key = f"v{SIM_VERSION}:{rel.key}:{iterations}:{min_cm}:{calibration}:{src}"
    blob = store.sim_get(key)
    if blob:
        return Distribution.from_json(blob)
    dist = simulate(rel, iterations, min_cm, gmap=gmap, calibration=calibration)
    store.sim_put(key, dist.to_json())
    return dist

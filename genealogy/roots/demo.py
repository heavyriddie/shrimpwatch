"""A synthetic family that is genetically real, for tests and for trust.

Two problems this solves.

The first is testing.  Every other module in this project takes real people's
DNA as input, and real DNA is exactly the thing you cannot commit to a
repository or paste into a bug report.  Without a fixture the test suite can
only check that functions run, not that they are *right*: does the segment
detector actually recover the segments that were inherited, does the
relationship predictor actually name the relationship that exists.  So this
module builds a family from the ground up -- founder haplotypes, meiosis with
crossovers, transmission down an explicit pedigree -- and writes out both the
raw files and the ground truth.  A test can then assert against a known
answer instead of against a golden output nobody can justify.

The second is trust.  Before pointing a tool at your own kit and believing
what it says about a stranger who shares 214 cM with you, you want to watch it
get a known answer right.  `generate` writes a complete worked example -- two
vendors' raw files for you and for your mother, a match list, a segment list,
shared-match data, and a GEDCOM -- so the whole pipeline can be run end to end
against numbers whose true values are printed in `truth.json`.

Why simulate rather than fabricate.  It would be far less code to draw
genotypes at random and write plausible-looking cM totals beside them, and the
result would be worthless: parent and child would violate Mendelian
inheritance at a third of all sites, and the "segments" would not correspond
to anything in the genotypes, so no detector could ever recover them.  Here
the numbers in the CSVs are *measurements of the simulation*.  The segments in
``segments_self.csv`` are the stretches SELF and the match genuinely inherited
from the same founder chromosome, the cM totals are those segments added up
above the detection threshold, and the raw files are the genotypes those same
haplotypes carry.  The files are internally consistent because they all
describe one underlying event.

The pedigree mirrors the situation the toolkit was written for: a tested
person, a tested mother, a documented maternal line, and an undocumented
paternal side known only through DNA matches.  The father is deliberately
*not* tested.  The paternal cluster -- a first cousin, two second cousins who
are also first cousins to each other, and a half first cousin -- is what makes
the demo worth clustering, and `maternal_tree.ged` deliberately stops at the
mother, because filling in the other half is the problem the tools exist to
attack.

Modelling choices worth knowing about:

* Meiosis is a Poisson process at one crossover per Morgan, no interference,
  on the uniform-rate genetic map -- the same model as `roots.dna.simulate`,
  whose interval algebra this module reuses rather than reimplements.  With a
  uniform map, cM and bp are proportional within a chromosome, so a crossover
  drawn in cM lands unambiguously between two SNPs.
* Allele frequencies come from a folded Beta(0.5, 0.5): the arcsine density
  is the neutral site-frequency spectrum's usual stand-in, and folding it onto
  [0, 0.5] gives the minor allele the low-frequency skew real arrays show.  It
  is not a real spectrum and no allele frequency here means anything; what
  matters is that heterozygosity comes out near the real ~25-30%, because that
  is what determines how informative a site is for IBD detection.
* Genotyping errors flip one allele of a call.  That is the failure mode that
  matters, because it is the only way a true parent and child can end up as
  opposite homozygotes -- which is precisely the signal every downstream
  mismatch-bunching rule is written to survive.  A test that never sees one is
  not testing anything.
* The two vendor panels overlap only partially, by construction, so merging
  two kits from the same person actually has work to do.
* ``Side`` in the match lists is the *true* side taken from the pedigree, not
  a guess.  In real use it is inferred by bucketing matches against a tested
  parent; here it is the answer that inference is supposed to arrive at.
"""

from __future__ import annotations

import bisect
import csv
import json
import math
import os
import random
from array import array
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from .dna.simulate import (
    Interval,
    Span,
    _coalesce,
    _matching,
    _slice,
    _union,
    relationship_for,
)
from .genome import AUTOSOMES, CHROM_BP, CHROM_CM, GeneticMap

# The interval helpers above are private to `simulate`, but they are the
# tested implementation of exactly this project's IBD algebra; duplicating
# them here would mean two implementations that could drift apart.

ALLELES = "ACGT"
CHROM_ORDER: Tuple[str, ...] = AUTOSOMES + ("X", "Y", "MT")

#: Ancestry's numeric chromosome coding for the non-numeric chromosomes.
_ANCESTRY_CHROM: Dict[str, str] = {"X": "23", "Y": "24", "MT": "26"}

_MONTHS = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN",
           "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

#: Y and MT do not recombine, so they need no genetic length; provenance
#: intervals on them are placeholders and never reach a cM report.
_NONRECOMBINING_CM = 1.0


# ---------------------------------------------------------------------------
# configuration


@dataclass
class DemoConfig:
    """Every knob the generator has.

    ``snp_count`` is the size of the *master* panel across the 22 autosomes
    plus X, split between chromosomes in proportion to physical length.  Each
    vendor kit gets a subset of it (see ``panel_overlap``), so an individual
    file is smaller than this number.  Y and MT panels are sized separately
    because their real SNP counts bear no relation to their physical length.
    """

    seed: int = 20260730
    build: str = "37"

    snp_count: int = 640_000
    y_snps: int = 1_400
    mt_snps: int = 2_000

    error_rate: float = 0.002
    no_call_rate: float = 0.005
    #: Fraction of one vendor kit's positions that also appear in the other.
    panel_overlap: float = 0.65

    min_segment_cm: float = 7.0
    #: Minimum SNPs in a reported segment.  None means "500 at the default
    #: density, scaled with ``snp_count``".  The threshold exists to reject
    #: long stretches that clear 7 cM only because the array barely covers
    #: them, so it is a statement about array density, not about biology --
    #: pinning it at 500 while dropping snp_count for a fast test would throw
    #: away every genuine cousin segment.
    min_segment_snps: Optional[int] = None

    include_sibling_kit: bool = True

    filenames: Dict[str, str] = field(
        default_factory=lambda: {
            "self_23andme": "self_23andme.txt",
            "self_ancestry": "self_ancestry.txt",
            "mother_23andme": "mother_23andme.txt",
            "mother_ancestry": "mother_ancestry.txt",
            "sibling_23andme": "sibling_23andme.txt",
            "matches_self": "matches_self.csv",
            "matches_mother": "matches_mother.csv",
            "segments_self": "segments_self.csv",
            "icw_self": "icw_self.csv",
            "maternal_tree": "maternal_tree.ged",
            "truth": "truth.json",
        }
    )

    @property
    def segment_snp_threshold(self) -> int:
        if self.min_segment_snps is not None:
            return self.min_segment_snps
        return max(20, int(round(500 * self.snp_count / 640_000)))


# ---------------------------------------------------------------------------
# people


@dataclass
class _Genome:
    """One person's chromosomes.

    ``hap1`` carries what came from the father, ``hap2`` what came from the
    mother; each is a per-chromosome ``array('B')`` of 0 (reference allele) or
    1 (alternate), indexed the same way as the master panel.  A male has no
    paternal X (``hap1['X']`` is None) and no maternal Y.  ``prov1``/``prov2``
    are the matching founder-haplotype provenance in centimorgans, in the same
    ``(start, end, label)`` form `roots.dna.simulate` uses -- this is what
    makes shared segments computable rather than guessable.
    """

    hap1: Dict[str, Optional[array]] = field(default_factory=dict)
    hap2: Dict[str, Optional[array]] = field(default_factory=dict)
    prov1: Dict[str, List[Interval]] = field(default_factory=dict)
    prov2: Dict[str, List[Interval]] = field(default_factory=dict)


@dataclass
class DemoPerson:
    """A simulated person.

    ``genome`` is excluded from repr and from every serialised form; it is
    tens of megabytes and means nothing outside this module.
    """

    label: str
    given: str
    surname: str
    sex: str
    birth_year: int
    birth_place: str
    father: Optional[str] = None
    mother: Optional[str] = None
    tested: bool = False

    birth_day: int = 1
    birth_month: int = 1

    relationship_to_self: str = ""
    relationship_to_mother: str = ""
    side: str = "Unknown"

    genome: Optional[_Genome] = field(default=None, repr=False, compare=False)

    @property
    def name(self) -> str:
        return f"{self.given} {self.surname}"

    @property
    def gedcom_name(self) -> str:
        return f"{self.given} /{self.surname}/"

    @property
    def birth_date(self) -> str:
        return f"{self.birth_day} {_MONTHS[self.birth_month - 1]} {self.birth_year}"

    def public(self) -> Dict[str, object]:
        return {
            "label": self.label,
            "name": self.name,
            "sex": self.sex,
            "birth_year": self.birth_year,
            "birth_place": self.birth_place,
            "father": self.father,
            "mother": self.mother,
            "tested": self.tested,
            "relationship_to_self": self.relationship_to_self,
            "relationship_to_mother": self.relationship_to_mother,
            "side": self.side,
        }


# label, given, surname, sex, birth year, birth place, father, mother, tested
#
# Order matters: a person's parents must appear before they do.  "tested"
# marks the people who show up as DNA matches -- note that the father and both
# paternal grandparents are untested, which is the whole point of the demo.
_ROSTER: Tuple[Tuple, ...] = (
    # -- maternal line ------------------------------------------------------
    ("mggf", "Alfred", "Whitlock", "M", 1901, "Bristol, England", None, None, False),
    ("mggm", "Edith", "Ramsay", "F", 1904, "Bristol, England", None, None, False),
    ("mgf", "Harold", "Whitlock", "M", 1930, "Portland, Oregon, USA", "mggf", "mggm", False),
    ("mgm", "Marjorie", "Pike", "F", 1933, "Salem, Oregon, USA", None, None, False),
    ("mgreataunt", "Constance", "Whitlock", "F", 1935, "Portland, Oregon, USA", "mggf", "mggm", True),
    ("mgreataunt_spouse", "Frank", "Ellery", "M", 1932, "Boise, Idaho, USA", None, None, False),
    ("mother", "Diane", "Whitlock", "F", 1958, "Portland, Oregon, USA", "mgf", "mgm", True),
    ("maunt", "Susan", "Whitlock", "F", 1961, "Portland, Oregon, USA", "mgf", "mgm", True),
    ("maunt_spouse", "Gerald", "Doyle", "M", 1959, "Eugene, Oregon, USA", None, None, False),
    ("m1c1r", "Paul", "Ellery", "M", 1962, "Boise, Idaho, USA", "mgreataunt_spouse", "mgreataunt", False),
    ("m1c1r_spouse", "Nadine", "Krause", "F", 1964, "Boise, Idaho, USA", None, None, False),
    ("mcousin1", "Kevin", "Doyle", "M", 1989, "Eugene, Oregon, USA", "maunt_spouse", "maunt", True),
    ("mcousin2", "Rachel", "Ellery", "F", 1991, "Boise, Idaho, USA", "m1c1r", "m1c1r_spouse", True),
    # -- paternal line ------------------------------------------------------
    ("pggf", "Walter", "Ferris", "M", 1899, "Tulsa, Oklahoma, USA", None, None, False),
    ("pggm", "Ivy", "Marsh", "F", 1903, "Tulsa, Oklahoma, USA", None, None, False),
    ("pgf", "Raymond", "Ferris", "M", 1929, "Tulsa, Oklahoma, USA", "pggf", "pggm", False),
    ("pgm", "Eleanor", "Nash", "F", 1934, "Wichita, Kansas, USA", None, None, False),
    ("pgf_partner", "Lorraine", "Beckett", "F", 1927, "Tulsa, Oklahoma, USA", None, None, False),
    ("pgreatuncle", "Douglas", "Ferris", "M", 1932, "Tulsa, Oklahoma, USA", "pggf", "pggm", False),
    ("pgreatuncle_spouse", "Norma", "Sandoval", "F", 1935, "Wichita, Kansas, USA", None, None, False),
    ("phalfuncle", "Bernard", "Ferris", "M", 1950, "Tulsa, Oklahoma, USA", "pgf", "pgf_partner", False),
    ("phalfuncle_spouse", "Yvette", "Marchetti", "F", 1953, "Kansas City, Missouri, USA", None, None, False),
    ("father", "Martin", "Ferris", "M", 1956, "Tulsa, Oklahoma, USA", "pgf", "pgm", False),
    ("paunt", "Gloria", "Ferris", "F", 1959, "Tulsa, Oklahoma, USA", "pgf", "pgm", False),
    ("paunt_spouse", "Curtis", "Vaughn", "M", 1957, "Norman, Oklahoma, USA", None, None, False),
    ("p1c1r_a", "Janet", "Ferris", "F", 1960, "Tulsa, Oklahoma, USA", "pgreatuncle", "pgreatuncle_spouse", False),
    ("p1c1r_a_spouse", "Dennis", "Kane", "M", 1958, "Wichita, Kansas, USA", None, None, False),
    ("p1c1r_b", "Roy", "Ferris", "M", 1963, "Tulsa, Oklahoma, USA", "pgreatuncle", "pgreatuncle_spouse", False),
    ("p1c1r_b_spouse", "Sheila", "Otto", "F", 1965, "Wichita, Kansas, USA", None, None, False),
    ("pcousin1", "Terrence", "Vaughn", "M", 1987, "Norman, Oklahoma, USA", "paunt_spouse", "paunt", True),
    ("pcousin2a", "Melissa", "Kane", "F", 1990, "Wichita, Kansas, USA", "p1c1r_a_spouse", "p1c1r_a", True),
    ("pcousin2b", "Owen", "Ferris", "M", 1992, "Tulsa, Oklahoma, USA", "p1c1r_b", "p1c1r_b_spouse", True),
    ("phalfcousin1", "Dana", "Ferris", "F", 1979, "Kansas City, Missouri, USA", "phalfuncle", "phalfuncle_spouse", True),
    # -- the proband --------------------------------------------------------
    ("self", "Ada", "Whitlock", "F", 1986, "Portland, Oregon, USA", "father", "mother", True),
    ("sibling", "Nolan", "Whitlock", "M", 1983, "Portland, Oregon, USA", "father", "mother", True),
)

SELF = "self"
MOTHER = "mother"

#: The documented half of the tree.  Everything else is the mystery.
_GEDCOM_LABELS: Tuple[str, ...] = (
    "mggf", "mggm", "mgf", "mgm", "mother", "maunt", "mcousin1", "self",
)


# ---------------------------------------------------------------------------
# SNP panel


@dataclass
class _ChromPanel:
    """The master SNP layout for one chromosome."""

    chrom: str
    rsid_base: int
    positions: array  # 'l', ascending, unique
    ref: array  # 'B', index into ALLELES
    alt: array  # 'B', index into ALLELES
    freq: array  # 'f', alternate-allele frequency
    membership: array  # 'B', bit 1 = 23andMe panel, bit 2 = AncestryDNA panel
    idx_23andme: array = field(default_factory=lambda: array("i"))
    idx_ancestry: array = field(default_factory=lambda: array("i"))

    def __len__(self) -> int:
        return len(self.positions)

    def count_between(self, start_bp: int, end_bp: int) -> int:
        lo = bisect.bisect_left(self.positions, start_bp)
        hi = bisect.bisect_right(self.positions, end_bp)
        return hi - lo


Panel = Dict[str, _ChromPanel]


def _panel_sizes(cfg: DemoConfig) -> Dict[str, int]:
    """Split the SNP budget across chromosomes by physical length."""
    lengths = CHROM_BP[cfg.build]
    covered = AUTOSOMES + ("X",)
    total_bp = sum(lengths[c] for c in covered)
    sizes = {c: max(100, int(round(cfg.snp_count * lengths[c] / total_bp))) for c in covered}
    sizes["Y"] = min(cfg.y_snps, lengths["Y"] - 1)
    sizes["MT"] = min(cfg.mt_snps, lengths["MT"] - 1)
    return sizes


def _build_panel(cfg: DemoConfig, rng: random.Random) -> Panel:
    """Lay out SNPs, pick alleles and frequencies, assign vendor panels.

    Vendor membership is drawn so that the fraction of one kit's positions
    that also appear in the other kit is exactly ``panel_overlap``.  With
    ``c`` the chance a SNP is on both panels and ``e`` the chance it is
    exclusive to one, ``c / (c + e) == overlap``; forcing ``c + 2e == 1`` (so
    that no SNP is wasted) pins ``c = overlap / (2 - overlap)``.
    """
    overlap = min(1.0, max(1e-6, cfg.panel_overlap))
    p_both = overlap / (2.0 - overlap)
    p_one = p_both * (1.0 - overlap) / overlap

    lengths = CHROM_BP[cfg.build]
    sizes = _panel_sizes(cfg)
    panel: Panel = {}
    rsid_base = 1_000_000
    rnd = rng.random
    beta = rng.betavariate

    for chrom in CHROM_ORDER:
        n = sizes[chrom]
        # sample() guarantees distinct positions: a duplicated (chrom, pos)
        # is a real defect that raw.read_with_stats counts, and the demo
        # should not be shipping one.
        positions = array("l", sorted(rng.sample(range(1, lengths[chrom] + 1), n)))

        ref = array("B", bytes(n))
        alt = array("B", bytes(n))
        freq = array("f", [0.0]) * n
        membership = array("B", bytes(n))
        for i in range(n):
            a = int(rnd() * 4.0)
            b = (a + 1 + int(rnd() * 3.0)) & 3
            ref[i] = a
            alt[i] = b
            # Folded arcsine: minor-allele frequency piled up near zero, with
            # a long thin tail toward 0.5, which is roughly what real arrays
            # look like once rare sites have been dropped.
            x = beta(0.5, 0.5)
            maf = x if x <= 0.5 else 1.0 - x
            freq[i] = maf if maf >= 0.01 else 0.01
            u = rnd()
            if u < p_both:
                membership[i] = 3
            elif u < p_both + p_one:
                membership[i] = 1
            else:
                membership[i] = 2

        cp = _ChromPanel(
            chrom=chrom,
            rsid_base=rsid_base,
            positions=positions,
            ref=ref,
            alt=alt,
            freq=freq,
            membership=membership,
        )
        cp.idx_23andme = array("i", [i for i in range(n) if membership[i] & 1])
        cp.idx_ancestry = array("i", [i for i in range(n) if membership[i] & 2])
        panel[chrom] = cp
        rsid_base += n
    return panel


# ---------------------------------------------------------------------------
# inheritance


class _Meiosis:
    """Founder creation and gamete formation over the whole SNP panel.

    Founder-haplotype labels are handed out per chromosome, so two founders
    never collide and a person's own two copies never look identical by
    accident.
    """

    def __init__(self, panel: Panel, gmap: GeneticMap, rng: random.Random):
        self.panel = panel
        self.gmap = gmap
        self.rng = rng
        self._next_label: Dict[str, int] = {c: 0 for c in CHROM_ORDER}
        self.lengths: Dict[str, float] = {
            c: (gmap.chrom_cm(c) if c in CHROM_CM else _NONRECOMBINING_CM)
            for c in CHROM_ORDER
        }

    # -- founders -------------------------------------------------------

    def _label(self, chrom: str) -> int:
        lab = self._next_label[chrom]
        self._next_label[chrom] = lab + 1
        return lab

    def _fresh(self, chrom: str) -> Tuple[array, List[Interval]]:
        cp = self.panel[chrom]
        rnd = self.rng.random
        hap = array("B", [1 if rnd() < f else 0 for f in cp.freq])
        return hap, [(0.0, self.lengths[chrom], self._label(chrom))]

    def founder(self, sex: str) -> _Genome:
        g = _Genome()
        for chrom in AUTOSOMES:
            g.hap1[chrom], g.prov1[chrom] = self._fresh(chrom)
            g.hap2[chrom], g.prov2[chrom] = self._fresh(chrom)
        # A male carries one X (from his mother, hence slot 2) and one Y.
        g.hap2["X"], g.prov2["X"] = self._fresh("X")
        if sex == "F":
            g.hap1["X"], g.prov1["X"] = self._fresh("X")
            g.hap1["Y"], g.prov1["Y"] = None, []
        else:
            g.hap1["X"], g.prov1["X"] = None, []
            g.hap1["Y"], g.prov1["Y"] = self._fresh("Y")
        g.hap2["Y"], g.prov2["Y"] = None, []
        g.hap1["MT"], g.prov1["MT"] = None, []
        g.hap2["MT"], g.prov2["MT"] = self._fresh("MT")
        return g

    # -- gametes --------------------------------------------------------

    def gamete(self, chrom: str, g: _Genome) -> Tuple[array, List[Interval]]:
        """One recombined chromosome, alleles and provenance in step.

        Crossovers arrive as a Poisson process at one per Morgan, so the gaps
        between them are exponential with mean 100 cM.  Each crossover's cM
        coordinate is converted to a base-pair coordinate and then to a SNP
        index, which is what keeps the emitted genotypes consistent with the
        segment boundaries reported in the CSVs.
        """
        rng = self.rng
        length = self.lengths[chrom]
        positions = self.panel[chrom].positions
        bp_at_cm = self.gmap.bp_at_cm

        points: List[float] = []
        x = rng.expovariate(0.01)
        while x < length:
            points.append(x)
            x += rng.expovariate(0.01)

        bounds_cm = [0.0] + points + [length]
        bounds_idx = [0]
        for cm in points:
            bounds_idx.append(bisect.bisect_left(positions, bp_at_cm(chrom, cm)))
        bounds_idx.append(len(positions))

        hap_a, hap_b = g.hap1[chrom], g.hap2[chrom]
        prov_a, prov_b = g.prov1[chrom], g.prov2[chrom]
        cur = rng.getrandbits(1)
        hap = array("B")
        prov: List[Interval] = []
        for i in range(len(bounds_cm) - 1):
            src_h = hap_a if cur == 0 else hap_b
            src_p = prov_a if cur == 0 else prov_b
            hap.extend(src_h[bounds_idx[i]:bounds_idx[i + 1]])
            prov.extend(_slice(src_p, bounds_cm[i], bounds_cm[i + 1]))
            cur ^= 1
        return hap, _coalesce(prov)

    def child(self, father: _Genome, mother: _Genome, sex: str) -> _Genome:
        g = _Genome()
        for chrom in AUTOSOMES:
            g.hap1[chrom], g.prov1[chrom] = self.gamete(chrom, father)
            g.hap2[chrom], g.prov2[chrom] = self.gamete(chrom, mother)

        # The X is the whole reason sex has to be modelled.  A father has one
        # X and passes it whole -- no recombination partner -- to a daughter,
        # and passes his Y to a son instead.  A mother recombines her two X
        # copies exactly as she does an autosome.
        g.hap2["X"], g.prov2["X"] = self.gamete("X", mother)
        if sex == "F":
            g.hap1["X"] = array("B", father.hap2["X"])
            g.prov1["X"] = list(father.prov2["X"])
            g.hap1["Y"], g.prov1["Y"] = None, []
        else:
            g.hap1["X"], g.prov1["X"] = None, []
            g.hap1["Y"] = array("B", father.hap1["Y"])
            g.prov1["Y"] = list(father.prov1["Y"])
        g.hap2["Y"], g.prov2["Y"] = None, []
        g.hap1["MT"], g.prov1["MT"] = None, []
        g.hap2["MT"] = array("B", mother.hap2["MT"])
        g.prov2["MT"] = list(mother.prov2["MT"])
        return g


def _simulate_people(
    people: Dict[str, DemoPerson], panel: Panel, gmap: GeneticMap, rng: random.Random
) -> None:
    """Give everybody a genome, parents before children."""
    meiosis = _Meiosis(panel, gmap, rng)
    for label, person in people.items():
        if person.father is None and person.mother is None:
            person.genome = meiosis.founder(person.sex)
        else:
            father = people[person.father].genome
            mother = people[person.mother].genome
            assert father is not None and mother is not None, f"{label}: parents unsimulated"
            person.genome = meiosis.child(father, mother, person.sex)


# ---------------------------------------------------------------------------
# pedigree arithmetic


def _ancestors(people: Dict[str, DemoPerson], label: str) -> Dict[str, int]:
    """Every ancestor of ``label`` mapped to its shortest generation depth.

    Depth 0 is the person themself, which is what makes direct lines fall out
    of the same most-recent-common-ancestor calculation as everything else.
    """
    depths: Dict[str, int] = {}
    frontier = [(label, 0)]
    while frontier:
        cur, d = frontier.pop()
        if cur in depths and depths[cur] <= d:
            continue
        depths[cur] = d
        person = people[cur]
        for parent in (person.father, person.mother):
            if parent is not None:
                frontier.append((parent, d + 1))
    return depths


def _mrcas(
    people: Dict[str, DemoPerson], a: str, b: str
) -> Tuple[List[str], int, int]:
    """Most recent common ancestors of two people, with the meiosis counts.

    Returns ``(labels, up, down)`` where ``up`` is the number of generations
    from ``a`` to the common ancestors and ``down`` from ``b``.  Ties are kept
    together: a full relationship descends from two ancestors, a half
    relationship from one, and that count is what distinguishes them.
    """
    da, db = _ancestors(people, a), _ancestors(people, b)
    common = set(da) & set(db)
    if not common:
        return [], 0, 0
    best = min(da[c] + db[c] for c in common)
    labels = sorted(c for c in common if da[c] + db[c] == best)
    return labels, da[labels[0]], db[labels[0]]


def _relationship(people: Dict[str, DemoPerson], viewer: str, other: str) -> str:
    labels, up, down = _mrcas(people, viewer, other)
    if not labels:
        return "unrelated"
    return relationship_for(up, down, len(labels)).name


def _side(people: Dict[str, DemoPerson], viewer: str, other: str) -> str:
    """Which of the viewer's parents the relationship runs through."""
    labels, _up, _down = _mrcas(people, viewer, other)
    if not labels:
        return "Unknown"
    if viewer in labels:
        return "Direct line"
    v = people[viewer]
    pat = _ancestors(people, v.father) if v.father else {}
    mat = _ancestors(people, v.mother) if v.mother else {}
    on_pat = any(c in pat for c in labels)
    on_mat = any(c in mat for c in labels)
    if on_pat and on_mat:
        return "Both"
    if on_mat:
        return "Maternal"
    if on_pat:
        return "Paternal"
    return "Unknown"


# ---------------------------------------------------------------------------
# segments


@dataclass
class DemoSegment:
    """One stretch of chromosome two people inherited from the same founder."""

    chrom: str
    start_bp: int
    end_bp: int
    cm: float
    snps: int


@dataclass
class DemoMatch:
    """What a testing company would show for one relative."""

    label: str
    name: str
    total_cm: float
    segments: List[DemoSegment]
    true_ibd_cm: float
    relationship: str
    side: str

    @property
    def longest_cm(self) -> float:
        return max((s.cm for s in self.segments), default=0.0)

    @property
    def segment_count(self) -> int:
        return len(self.segments)


def _shared(
    a: _Genome, b: _Genome, panel: Panel, gmap: GeneticMap, cfg: DemoConfig
) -> Tuple[List[DemoSegment], float]:
    """Autosomal half-identical regions between two simulated genomes.

    Returns the segments that clear both thresholds, plus the *un*thresholded
    total, because the gap between the two is exactly the "small segments the
    companies cannot see" effect that makes distant relatives unreliable.
    """
    kept: List[DemoSegment] = []
    true_total = 0.0
    min_snps = cfg.segment_snp_threshold
    for chrom in AUTOSOMES:
        spans: List[Span] = _union(
            _matching(a.prov1[chrom], b.prov1[chrom])
            + _matching(a.prov1[chrom], b.prov2[chrom])
            + _matching(a.prov2[chrom], b.prov1[chrom])
            + _matching(a.prov2[chrom], b.prov2[chrom])
        )
        cp = panel[chrom]
        for start_cm, end_cm in spans:
            cm = end_cm - start_cm
            true_total += cm
            if cm < cfg.min_segment_cm:
                continue
            start_bp = gmap.bp_at_cm(chrom, start_cm)
            end_bp = gmap.bp_at_cm(chrom, end_cm)
            snps = cp.count_between(start_bp, end_bp)
            if snps < min_snps:
                continue
            kept.append(DemoSegment(chrom, start_bp, end_bp, cm, snps))
    return kept, true_total


def _match_list(
    people: Dict[str, DemoPerson],
    viewer: str,
    panel: Panel,
    gmap: GeneticMap,
    cfg: DemoConfig,
) -> List[DemoMatch]:
    """Everyone tested who shares detectable DNA with ``viewer``."""
    vg = people[viewer].genome
    assert vg is not None
    out: List[DemoMatch] = []
    for label, person in people.items():
        if label == viewer or not person.tested:
            continue
        assert person.genome is not None
        segs, true_total = _shared(vg, person.genome, panel, gmap, cfg)
        if not segs:
            continue
        out.append(
            DemoMatch(
                label=label,
                name=person.name,
                total_cm=sum(s.cm for s in segs),
                segments=segs,
                true_ibd_cm=true_total,
                relationship=_relationship(people, viewer, label),
                side=_side(people, viewer, label),
            )
        )
    out.sort(key=lambda m: -m.total_cm)
    return out


def _predicted_relationship(cm: float) -> str:
    """The coarse bucket a testing company would print.

    Deliberately not `roots.dna.predict`: the point of the demo is to give the
    real predictor something to improve on, so this is the crude cM-range
    lookup the vendors ship, ranges and all.
    """
    if cm >= 3300:
        return "Parent / Child"
    if cm >= 2100:
        return "Sibling"
    if cm >= 1300:
        return "Close family - 1st cousin"
    if cm >= 850:
        return "1st cousin"
    if cm >= 400:
        return "1st - 2nd cousin"
    if cm >= 200:
        return "2nd cousin"
    if cm >= 90:
        return "2nd - 3rd cousin"
    if cm >= 40:
        return "3rd - 4th cousin"
    return "Distant cousin"


# ---------------------------------------------------------------------------
# raw file emission


def _sparse_indices(n: int, p: float, rng: random.Random) -> set:
    """Indices in ``range(n)`` independently selected with probability ``p``.

    Drawn by jumping a geometric number of positions at a time rather than
    rolling a die per SNP.  At the default rates that is ~3000 draws per
    chromosome instead of ~30000, which is the difference between the writer
    being noticeable and not.
    """
    if p <= 0.0 or n <= 0:
        return set()
    if p >= 1.0:
        return set(range(n))
    step = math.log1p(-p)
    out = set()
    i = -1
    while True:
        # 1 - random() lands in (0, 1], so the log is always defined.
        i += 1 + int(math.log(1.0 - rng.random()) / step)
        if i >= n:
            return out
        out.add(i)


def _calls(
    person: DemoPerson, panel: Panel, vendor_bit: int, cfg: DemoConfig, rng: random.Random
) -> Iterator[Tuple[str, str, int, Optional[str], Optional[str]]]:
    """Yield ``(rsid, chrom, pos, allele1, allele2)`` for one kit.

    ``allele2`` is None where the site is hemizygous (a male's X and Y, and
    everybody's mitochondrion); both are None for a no-call.  Errors and
    no-calls are drawn per kit, so the same person's two vendor files disagree
    in different places -- which is the only reason cross-kit concordance
    checking is worth running.
    """
    g = person.genome
    assert g is not None
    for chrom in CHROM_ORDER:
        cp = panel[chrom]
        idx = cp.idx_23andme if vendor_bit == 1 else cp.idx_ancestry
        n = len(idx)
        errors = _sparse_indices(n, cfg.error_rate, rng)
        no_calls = _sparse_indices(n, cfg.no_call_rate, rng)

        h1, h2 = g.hap1[chrom], g.hap2[chrom]
        if chrom == "MT":
            hemi, diploid = h2, None
        elif chrom == "Y":
            hemi, diploid = h1, None
        elif chrom == "X" and h1 is None:
            hemi, diploid = h2, None
        else:
            hemi, diploid = None, (h1, h2)

        base = cp.rsid_base
        positions = cp.positions
        ref, alt = cp.ref, cp.alt
        for j in range(n):
            i = idx[j]
            rsid = "rs" + str(base + i)
            pos = positions[i]
            if j in no_calls or (hemi is None and diploid is None):
                yield rsid, chrom, pos, None, None
                continue
            flip = j in errors
            if hemi is not None:
                v = hemi[i]
                if flip:
                    v ^= 1
                yield rsid, chrom, pos, ALLELES[ref[i] if v == 0 else alt[i]], None
                continue
            v1, v2 = diploid[0][i], diploid[1][i]
            if flip:
                # A single flipped allele is what turns a het into a
                # homozygote and so manufactures the opposite-homozygote
                # "Mendelian violation" that downstream error handling exists
                # to absorb.
                if rng.getrandbits(1):
                    v1 ^= 1
                else:
                    v2 ^= 1
            yield (
                rsid,
                chrom,
                pos,
                ALLELES[ref[i] if v1 == 0 else alt[i]],
                ALLELES[ref[i] if v2 == 0 else alt[i]],
            )


_23ANDME_HEADER = """\
# This data file was generated by roots.demo, not by 23andMe.  The genotypes
# below are simulated and belong to no living person; the layout copies the
# 23andMe raw data export so that tools can be exercised against it.
#
# Below is a text version of your data.  Fields are TAB-separated.
# Each line corresponds to a single SNP.  For each SNP, we provide its
# identifier (an rsid), its location on the reference human genome, and the
# genotype call oriented with respect to the plus strand on the human
# reference sequence.
#
# reference human assembly build 37 (also known as Annotation Release 104)
#
# rsid\tchromosome\tposition\tgenotype
"""

_ANCESTRY_HEADER = """\
#AncestryDNA raw data download
#This file was generated by roots.demo.  It is synthetic data in the
#AncestryDNA export layout and contains no real human genotypes.
#
#Data was collected using AncestryDNA array version: V2.0
#Data is formatted using human reference build 37.
#
#Below is a text version of the data.  Fields are TAB-separated.
"""


def _write_23andme(path: str, person: DemoPerson, panel: Panel, cfg: DemoConfig,
                   rng: random.Random) -> int:
    """Write a 23andMe-format kit; returns the number of SNP rows."""
    written = 0
    chunk: List[str] = []
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(_23ANDME_HEADER)
        for rsid, chrom, pos, a1, a2 in _calls(person, panel, 1, cfg, rng):
            if a1 is None:
                gt = "--"
            elif a2 is None:
                gt = a1  # hemizygous: one letter, the way the vendor does it
            else:
                gt = a1 + a2
            chunk.append(f"{rsid}\t{chrom}\t{pos}\t{gt}")
            written += 1
            if len(chunk) >= 50_000:
                fh.write("\n".join(chunk))
                fh.write("\n")
                chunk.clear()
        if chunk:
            fh.write("\n".join(chunk))
            fh.write("\n")
    return written


def _write_ancestry(path: str, person: DemoPerson, panel: Panel, cfg: DemoConfig,
                    rng: random.Random) -> int:
    """Write an AncestryDNA-format kit; returns the number of SNP rows.

    Ancestry has two allele columns and no way to spell a hemizygous call, so
    a male's X and Y and everyone's mitochondrion get the single allele
    written twice -- which is what the real export does.
    """
    written = 0
    chunk: List[str] = []
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(_ANCESTRY_HEADER)
        fh.write("rsid\tchromosome\tposition\tallele1\tallele2\n")
        for rsid, chrom, pos, a1, a2 in _calls(person, panel, 2, cfg, rng):
            if a1 is None:
                c1 = c2 = "0"
            elif a2 is None:
                c1 = c2 = a1
            else:
                c1, c2 = a1, a2
            code = _ANCESTRY_CHROM.get(chrom, chrom)
            chunk.append(f"{rsid}\t{code}\t{pos}\t{c1}\t{c2}")
            written += 1
            if len(chunk) >= 50_000:
                fh.write("\n".join(chunk))
                fh.write("\n")
                chunk.clear()
        if chunk:
            fh.write("\n".join(chunk))
            fh.write("\n")
    return written


# ---------------------------------------------------------------------------
# report emission


MATCH_HEADER = [
    "Name", "Shared cM", "Shared Segments", "Longest cM",
    "Predicted Relationship", "Side", "True Relationship",
]
SEGMENT_HEADER = [
    "Match Name", "Chromosome", "Start Position", "End Position",
    "Centimorgans", "Matching SNPs",
]
ICW_HEADER = ["Match Name A", "Match Name B"]


def _write_matches(path: str, matches: Sequence[DemoMatch]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(MATCH_HEADER)
        for m in matches:
            w.writerow([
                m.name,
                f"{m.total_cm:.1f}",
                m.segment_count,
                f"{m.longest_cm:.1f}",
                _predicted_relationship(m.total_cm),
                m.side,
                m.relationship,
            ])


def _write_segments(path: str, matches: Sequence[DemoMatch]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(SEGMENT_HEADER)
        for m in matches:
            for s in sorted(m.segments, key=lambda x: (int(x.chrom), x.start_bp)):
                w.writerow([m.name, s.chrom, s.start_bp, s.end_bp, f"{s.cm:.2f}", s.snps])


def _write_icw(
    path: str,
    matches: Sequence[DemoMatch],
    people: Dict[str, DemoPerson],
    panel: Panel,
    gmap: GeneticMap,
    cfg: DemoConfig,
) -> int:
    """Which of SELF's matches also match each other.

    This is the "in common with" grid, and it is the raw material for every
    clustering method: the paternal matches will fall into one block and the
    maternal ones into another, without anybody having said which is which.
    """
    rows: List[Tuple[str, str]] = []
    for i in range(len(matches)):
        gi = people[matches[i].label].genome
        assert gi is not None
        for j in range(i + 1, len(matches)):
            gj = people[matches[j].label].genome
            assert gj is not None
            segs, _ = _shared(gi, gj, panel, gmap, cfg)
            if segs:
                rows.append((matches[i].name, matches[j].name))
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(ICW_HEADER)
        for a, b in rows:
            w.writerow([a, b])
    return len(rows)


def _write_gedcom(path: str, people: Dict[str, DemoPerson], labels: Sequence[str]) -> None:
    """Emit the documented maternal side only.

    SELF's family record has a WIFE and a CHIL and no HUSB.  That absence is
    not a bug in the file; it is the research question the rest of the
    toolkit is pointed at.  Line endings are LF rather than the CRLF the
    5.5.1 spec asks for, because every parser in practice accepts LF and a
    text fixture with CRLF in it is a nuisance to diff.
    """
    ids = {lab: f"@I{i + 1}@" for i, lab in enumerate(labels)}
    included = set(labels)

    fam_order: List[Tuple[Optional[str], Optional[str]]] = []
    fam_children: Dict[Tuple[Optional[str], Optional[str]], List[str]] = {}
    for lab in labels:
        p = people[lab]
        f = p.father if p.father in included else None
        m = p.mother if p.mother in included else None
        if f is None and m is None:
            continue
        key = (f, m)
        if key not in fam_children:
            fam_children[key] = []
            fam_order.append(key)
        fam_children[key].append(lab)
    fam_ids = {k: f"@F{i + 1}@" for i, k in enumerate(fam_order)}

    famc: Dict[str, str] = {}
    fams: Dict[str, List[str]] = {}
    for key, kids in fam_children.items():
        fid = fam_ids[key]
        for kid in kids:
            famc[kid] = fid
        for spouse in key:
            if spouse is not None:
                fams.setdefault(spouse, []).append(fid)

    today = datetime.now(timezone.utc)
    out: List[str] = [
        "0 HEAD",
        "1 SOUR ROOTS",
        "2 NAME roots genealogy toolkit",
        "2 VERS demo",
        "1 DEST ANY",
        f"1 DATE {today.day} {_MONTHS[today.month - 1]} {today.year}",
        "1 SUBM @SUB1@",
        "1 FILE maternal_tree.ged",
        "1 GEDC",
        "2 VERS 5.5.1",
        "2 FORM LINEAGE-LINKED",
        "1 CHAR UTF-8",
        "1 NOTE Documented maternal line only; the paternal side is unknown.",
        "0 @SUB1@ SUBM",
        "1 NAME roots demo generator",
    ]
    for lab in labels:
        p = people[lab]
        out.append(f"0 {ids[lab]} INDI")
        out.append(f"1 NAME {p.gedcom_name}")
        out.append(f"2 GIVN {p.given}")
        out.append(f"2 SURN {p.surname}")
        out.append(f"1 SEX {p.sex}")
        out.append("1 BIRT")
        out.append(f"2 DATE {p.birth_date}")
        out.append(f"2 PLAC {p.birth_place}")
        if lab in famc:
            out.append(f"1 FAMC {famc[lab]}")
        for fid in fams.get(lab, []):
            out.append(f"1 FAMS {fid}")
    for key in fam_order:
        fid = fam_ids[key]
        husb, wife = key
        out.append(f"0 {fid} FAM")
        if husb is not None:
            out.append(f"1 HUSB {ids[husb]}")
        if wife is not None:
            out.append(f"1 WIFE {ids[wife]}")
        for kid in fam_children[key]:
            out.append(f"1 CHIL {ids[kid]}")
        if husb is not None and wife is not None:
            eldest = min(people[k].birth_year for k in fam_children[key])
            out.append("1 MARR")
            out.append(f"2 DATE {_MONTHS[(eldest % 12)]} {eldest - 2}")
    out.append("0 TRLR")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out))
        fh.write("\n")


# ---------------------------------------------------------------------------
# entry point


def _build_people(cfg: DemoConfig) -> Dict[str, DemoPerson]:
    rng = random.Random(cfg.seed ^ 0x5EED)
    people: Dict[str, DemoPerson] = {}
    for label, given, surname, sex, year, place, father, mother, tested in _ROSTER:
        people[label] = DemoPerson(
            label=label,
            given=given,
            surname=surname,
            sex=sex,
            birth_year=year,
            birth_place=place,
            father=father,
            mother=mother,
            tested=tested,
            birth_month=rng.randint(1, 12),
            birth_day=rng.randint(1, 28),
        )
    return people


def generate(outdir: str, config: Optional[DemoConfig] = None) -> Dict[str, str]:
    """Build the demo family and write every file.

    Returns a mapping of logical name (``'self_23andme'``, ``'truth'``, ...)
    to the path written.  The whole thing is a pure function of
    ``config.seed``: same seed, byte-identical output.
    """
    cfg = config or DemoConfig()
    os.makedirs(outdir, exist_ok=True)

    rng = random.Random(cfg.seed)
    gmap = GeneticMap.linear(cfg.build)

    panel = _build_panel(cfg, rng)
    people = _build_people(cfg)
    _simulate_people(people, panel, gmap, rng)

    for label in people:
        people[label].relationship_to_self = _relationship(people, SELF, label)
        people[label].relationship_to_mother = _relationship(people, MOTHER, label)
        people[label].side = _side(people, SELF, label)

    self_matches = _match_list(people, SELF, panel, gmap, cfg)
    mother_matches = _match_list(people, MOTHER, panel, gmap, cfg)

    paths: Dict[str, str] = {}

    def out(key: str) -> str:
        return os.path.join(outdir, cfg.filenames[key])

    # Raw kits.  Errors and no-calls are drawn from a writer-local generator
    # per kit so that adding or removing a kit does not perturb the others.
    kit_specs = [
        ("self_23andme", SELF, _write_23andme),
        ("self_ancestry", SELF, _write_ancestry),
        ("mother_23andme", MOTHER, _write_23andme),
        ("mother_ancestry", MOTHER, _write_ancestry),
    ]
    if cfg.include_sibling_kit:
        kit_specs.append(("sibling_23andme", "sibling", _write_23andme))

    kit_rows: Dict[str, int] = {}
    for n, (key, label, writer) in enumerate(kit_specs):
        path = out(key)
        kit_rows[key] = writer(path, people[label], panel, cfg, random.Random(cfg.seed + 1013 * (n + 1)))
        paths[key] = path

    _write_matches(out("matches_self"), self_matches)
    paths["matches_self"] = out("matches_self")
    _write_matches(out("matches_mother"), mother_matches)
    paths["matches_mother"] = out("matches_mother")
    _write_segments(out("segments_self"), self_matches)
    paths["segments_self"] = out("segments_self")
    icw_pairs = _write_icw(out("icw_self"), self_matches, people, panel, gmap, cfg)
    paths["icw_self"] = out("icw_self")
    _write_gedcom(out("maternal_tree"), people, _GEDCOM_LABELS)
    paths["maternal_tree"] = out("maternal_tree")

    truth = {
        "generator": "roots.demo",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config": asdict(cfg),
        "self": SELF,
        "mother": MOTHER,
        "detection": {
            "min_segment_cm": cfg.min_segment_cm,
            "min_segment_snps": cfg.segment_snp_threshold,
            "chromosomes": list(AUTOSOMES),
        },
        "panel": {
            "master_snps": sum(len(panel[c]) for c in CHROM_ORDER),
            "per_chromosome": {c: len(panel[c]) for c in CHROM_ORDER},
            "kit_rows": kit_rows,
        },
        "people": [people[label].public() for label in people],
        "edges": [
            {"child": label, "father": p.father, "mother": p.mother}
            for label, p in people.items()
            if p.father or p.mother
        ],
        "matches_self": [_match_truth(m, people, SELF) for m in self_matches],
        "matches_mother": [_match_truth(m, people, MOTHER) for m in mother_matches],
        "icw_pairs": icw_pairs,
        "files": {k: os.path.basename(v) for k, v in paths.items()},
    }
    truth_path = out("truth")
    with open(truth_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(truth, fh, indent=2, sort_keys=False)
        fh.write("\n")
    paths["truth"] = truth_path
    return paths


def _match_truth(m: DemoMatch, people: Dict[str, DemoPerson], viewer: str) -> Dict[str, object]:
    labels, up, down = _mrcas(people, viewer, m.label)
    return {
        "label": m.label,
        "name": m.name,
        "relationship": m.relationship,
        "side": m.side,
        "common_ancestors": labels,
        "up": up,
        "down": down,
        "detected_cm": round(m.total_cm, 2),
        "true_ibd_cm": round(m.true_ibd_cm, 2),
        "segments": m.segment_count,
        "longest_cm": round(m.longest_cm, 2),
        "predicted_relationship": _predicted_relationship(m.total_cm),
    }

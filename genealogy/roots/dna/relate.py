"""Detect IBD segments directly from two kits' genotypes.

This is the part that does not depend on any testing company: given raw data
for two people, find the stretches of chromosome where they are half- or
fully-identical, and add them up.  It is how the project verifies that the
mother kit really is the mother, and how any two kits in the project (yours
plus a cousin who sent you their download) get compared without uploading
anything anywhere.

Method, which is the standard one used by the third-party comparison tools:

1.  Walk the SNPs typed on both kits, in position order, and score each as
    IBS2 (identical genotypes), IBS1 (one allele in common), or IBS0
    (opposite homozygotes).
2.  IBS0 is impossible inside a shared segment, so IBS0 sites are the
    segment boundaries -- except that real data has genotyping errors.  An
    isolated IBS0 surrounded by long clean stretches is treated as an error;
    IBS0 sites that cluster are treated as real boundaries.  This is the
    "mismatch bunching" rule.
3.  Runs between surviving boundaries become candidate half-identical
    regions, kept if they clear both a genetic-length and a SNP-count
    threshold.  Both thresholds matter: a long stretch of low SNP density
    can clear 7 cM by accident, which is the classic false-positive segment.
4.  Inside each half-identical region, runs with no IBS1 at all are
    fully-identical regions -- both chromosome copies matching, which
    happens for full siblings and closer.

Thresholds default to 7 cM / 500 SNPs, comparable to what the testing
companies apply.  Below about 7 cM the majority of "segments" in consumer
data are not genuine recent inheritance, so lowering these is usually a way
to generate noise rather than discoveries.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from ..genome import AUTOSOMES, GeneticMap
from ..store import Store
from .qc import ibs


@dataclass
class Segment:
    chrom: str
    start_bp: int
    end_bp: int
    cm: float
    snps: int
    kind: str = "HIR"

    def as_row(self, kit_a: int, kit_b: int) -> Tuple:
        return (kit_a, kit_b, self.chrom, self.start_bp, self.end_bp, self.cm, self.snps, self.kind)

    def overlaps(self, other: "Segment", min_overlap_bp: int = 1) -> bool:
        if self.chrom != other.chrom:
            return False
        lo = max(self.start_bp, other.start_bp)
        hi = min(self.end_bp, other.end_bp)
        return (hi - lo) >= min_overlap_bp


@dataclass
class PairComparison:
    kit_a: int
    kit_b: int
    segments: List[Segment] = field(default_factory=list)
    shared_snps: int = 0
    ibs0: int = 0
    ibs1: int = 0
    ibs2: int = 0
    masked_errors: int = 0

    @property
    def hir(self) -> List[Segment]:
        return [s for s in self.segments if s.kind == "HIR"]

    @property
    def fir(self) -> List[Segment]:
        return [s for s in self.segments if s.kind == "FIR"]

    @property
    def total_cm(self) -> float:
        return sum(s.cm for s in self.hir)

    @property
    def total_fir_cm(self) -> float:
        return sum(s.cm for s in self.fir)

    @property
    def largest_cm(self) -> float:
        return max((s.cm for s in self.hir), default=0.0)

    @property
    def ibs0_rate(self) -> float:
        return self.ibs0 / self.shared_snps if self.shared_snps else 0.0


def compare(
    store: Store,
    kit_a: int,
    kit_b: int,
    gmap: GeneticMap,
    min_cm: float = 7.0,
    min_snps: int = 500,
    fir_min_cm: float = 5.0,
    fir_min_snps: int = 400,
    error_gap: int = 200,
    chroms: Optional[Sequence[str]] = None,
) -> PairComparison:
    """Compare two kits and return their shared segments.

    ``error_gap`` is the mismatch-bunching window: an IBS0 site with no
    other IBS0 within this many SNPs on either side is written off as a
    genotyping error instead of splitting a segment.
    """
    result = PairComparison(kit_a=kit_a, kit_b=kit_b)
    for chrom in (chroms or AUTOSOMES):
        positions: List[int] = []
        states: List[int] = []
        for pos, ga, gb in store.paired_genotypes(kit_a, kit_b, chrom):
            state = ibs(ga, gb)
            if state < 0:
                continue
            positions.append(pos)
            states.append(state)
        if len(positions) < min_snps:
            continue
        result.shared_snps += len(positions)
        result.ibs0 += states.count(0)
        result.ibs1 += states.count(1)
        result.ibs2 += states.count(2)

        boundaries, masked = _real_boundaries(states, 0, error_gap)
        result.masked_errors += masked
        for lo, hi in _runs_between(len(states), boundaries):
            seg = _make_segment(chrom, positions, lo, hi, gmap, "HIR")
            if seg.cm < min_cm or seg.snps < min_snps:
                continue
            result.segments.append(seg)
            # Fully-identical stretches live inside half-identical ones.
            sub_states = states[lo : hi + 1]
            sub_bounds, _ = _real_boundaries(sub_states, 1, error_gap)
            for slo, shi in _runs_between(len(sub_states), sub_bounds):
                fir = _make_segment(chrom, positions, lo + slo, lo + shi, gmap, "FIR")
                if fir.cm >= fir_min_cm and fir.snps >= fir_min_snps:
                    result.segments.append(fir)
    return result


def _real_boundaries(states: Sequence[int], breaker: int, error_gap: int) -> Tuple[List[int], int]:
    """Indices of breaking sites, after discarding isolated (error) ones."""
    hits = [i for i, s in enumerate(states) if s == breaker]
    if not hits:
        return [], 0
    kept: List[int] = []
    masked = 0
    for idx, i in enumerate(hits):
        prev_gap = i - hits[idx - 1] if idx > 0 else error_gap + 1
        next_gap = hits[idx + 1] - i if idx + 1 < len(hits) else error_gap + 1
        if prev_gap > error_gap and next_gap > error_gap:
            masked += 1
        else:
            kept.append(i)
    return kept, masked


def _runs_between(n: int, boundaries: Sequence[int]):
    """Yield inclusive index ranges of the gaps between boundary indices."""
    start = 0
    for b in boundaries:
        if b - 1 >= start:
            yield start, b - 1
        start = b + 1
    if start <= n - 1:
        yield start, n - 1


def _make_segment(
    chrom: str, positions: Sequence[int], lo: int, hi: int, gmap: GeneticMap, kind: str
) -> Segment:
    start_bp, end_bp = positions[lo], positions[hi]
    return Segment(
        chrom=chrom,
        start_bp=start_bp,
        end_bp=end_bp,
        cm=gmap.length_cm(chrom, start_bp, end_bp),
        snps=hi - lo + 1,
        kind=kind,
    )


# ---------------------------------------------------------------------------
# interpreting a comparison


@dataclass
class CloseVerdict:
    label: str
    confidence: str
    reasoning: str
    hir_fraction: float
    fir_fraction: float


def classify_close(cmp: PairComparison, gmap: GeneticMap) -> CloseVerdict:
    """Name the relationship for close pairs, where the pattern is decisive.

    The genome-wide *pattern* of sharing, not just the total, identifies
    everyone out to full siblings:

    * identical twins  -- everything shared on both copies
    * parent/child     -- everything shared on exactly one copy, and
                          crucially zero opposite-homozygote sites genome-wide
    * full siblings    -- about three quarters shared, a quarter of it on
                          both copies
    * half sib / grandparent / aunt / uncle -- about half, none on both copies

    Beyond that the pattern stops discriminating and you need the segment
    total plus genealogical context, which is what `roots.dna.predict` does.
    """
    genome = sum(gmap.chrom_cm(c) for c in AUTOSOMES) or 1.0
    hir_f = cmp.total_cm / genome
    fir_f = cmp.total_fir_cm / genome
    ibs0 = cmp.ibs0_rate

    if hir_f > 0.92 and fir_f > 0.85:
        return CloseVerdict(
            "identical twin or duplicate kit", "high",
            "the entire genome matches on both copies", hir_f, fir_f,
        )
    # Thresholds allow for real data: no-calls and imperfect segment
    # boundaries cost a few percent of genome-wide coverage, and genotyping
    # error puts opposite homozygotes at a few tenths of a percent even
    # between a true parent and child. What no other relationship can do is
    # combine near-total half-identity with an IBS0 rate this low; a full
    # sibling sits near 5%.
    if hir_f > 0.90 and ibs0 < 0.005:
        return CloseVerdict(
            "parent/child", "high",
            f"whole genome half-identical with essentially no opposite "
            f"homozygotes ({ibs0:.4%}), which only a parent-child pair produces",
            hir_f, fir_f,
        )
    if 0.60 <= hir_f <= 0.90 and fir_f > 0.10:
        return CloseVerdict(
            "full siblings", "high",
            f"{hir_f:.0%} half-identical with {fir_f:.0%} fully identical -- the "
            "signature of two children of the same two parents",
            hir_f, fir_f,
        )
    if 0.35 <= hir_f <= 0.60 and fir_f < 0.05:
        return CloseVerdict(
            "half sibling, grandparent/grandchild, or aunt/uncle-niece/nephew",
            "medium",
            "about half the genome shared on one copy only; these three "
            "relationships are indistinguishable by total sharing alone and "
            "need ages or known genealogy to separate",
            hir_f, fir_f,
        )
    return CloseVerdict(
        "more distant than second degree", "low",
        f"{cmp.total_cm:.0f} cM shared ({hir_f:.1%} of the autosomes); "
        "use relationship prediction rather than pattern matching",
        hir_f, fir_f,
    )


def verify_parent(cmp: PairComparison, gmap: GeneticMap) -> Tuple[bool, str]:
    """Confirm a claimed parent-child pair before phasing against it.

    Phasing assumes the parent really is a parent; if that assumption is
    wrong every downstream side assignment is wrong too, so this gate is
    worth failing loudly.
    """
    verdict = classify_close(cmp, gmap)
    if verdict.label == "parent/child":
        return True, verdict.reasoning
    return False, (
        f"expected a parent-child pattern but found: {verdict.label} "
        f"({verdict.reasoning})"
    )

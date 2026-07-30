"""Phase a child's genotypes against a tested parent.

This is the single highest-value thing you can do with a mother-child pair,
and it is why having your mother tested is worth more than testing yourself
twice.

At every site where you are heterozygous (say A/G) and your mother is
homozygous (A/A), she can only have given you A -- so the G came from your
father.  Repeat genome-wide and you have separated your two chromosome
copies: a maternal haplotype, and a paternal haplotype that is a direct
readout of half of a man who was never tested.  The only sites that stay
ambiguous are those where you are both heterozygous, since either of her
alleles could have been the transmitted one.

Two concrete payoffs:

* The paternal haplotype can be exported as a standalone pseudo-kit.  A
  match against it is a paternal-side match, full stop -- no cluster
  analysis, no guessing.  It also kills false segments, which arise when a
  comparison zig-zags between your maternal and paternal copies and reports
  a "shared" stretch that no single ancestor ever transmitted.  Those false
  segments dominate the small-segment range, so phased comparison lets you
  trust shorter matches than raw comparison ever can.

* The Mendelian error rate falls out for free and is a hard data-quality
  check: a genuine parent-child pair produces errors at a few tenths of a
  percent, all genotyping noise.  A percent or more means the pair is not
  parent and child.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

from ..genome import AUTOSOMES
from ..store import Store
from .raw import NO_CALL


@dataclass
class PhaseStats:
    sites: int = 0
    homozygous: int = 0
    het_resolved: int = 0
    het_ambiguous: int = 0
    parent_missing: int = 0
    mendelian_errors: int = 0
    x_maternal: int = 0
    per_chrom_resolved: Dict[str, int] = field(default_factory=dict)

    @property
    def het_total(self) -> int:
        return self.het_resolved + self.het_ambiguous

    @property
    def het_resolution_rate(self) -> float:
        return self.het_resolved / self.het_total if self.het_total else 0.0

    @property
    def mendelian_rate(self) -> float:
        return self.mendelian_errors / self.sites if self.sites else 0.0

    @property
    def phased_fraction(self) -> float:
        """Share of called sites where both parental alleles are known."""
        usable = self.homozygous + self.het_total
        return (self.homozygous + self.het_resolved) / usable if usable else 0.0

    def notes(self) -> List[str]:
        out: List[str] = []
        if self.mendelian_rate > 0.01:
            out.append(
                f"Mendelian error rate {self.mendelian_rate:.2%} is far above the "
                "~0.1-0.5% expected from genotyping noise -- these two kits are "
                "probably not parent and child. Do not trust anything phased "
                "from this pair."
            )
        elif self.mendelian_rate > 0.005:
            out.append(
                f"Mendelian error rate {self.mendelian_rate:.2%} is on the high "
                "side; check call rates on both kits."
            )
        if self.het_total and self.het_resolution_rate < 0.4:
            out.append(
                f"only {self.het_resolution_rate:.0%} of heterozygous sites "
                "resolved, which is low -- expect around 50-60%"
            )
        return out


def phase_against_parent(
    store: Store,
    child_kit: int,
    parent_kit: int,
    parent_role: str = "mother",
    child_sex: Optional[str] = None,
    persist: bool = True,
) -> PhaseStats:
    """Split a child's genotypes into the two parental haplotypes.

    ``parent_role`` decides which slot the tested parent's contribution
    fills; with a tested mother, the inferred haplotype is the father's.
    Returns counters; the per-site result is written to the ``phased`` table
    unless ``persist`` is false.
    """
    if parent_role not in ("mother", "father"):
        raise ValueError("parent_role must be 'mother' or 'father'")
    stats = PhaseStats()
    if persist:
        store.clear_phase(child_kit)

    def rows() -> Iterator[Tuple[str, int, Optional[str], Optional[str], str]]:
        chroms = list(AUTOSOMES) + ["X"]
        for chrom in chroms:
            for pos, child_gt, parent_gt in store.paired_genotypes(child_kit, parent_kit, chrom):
                if child_gt == NO_CALL:
                    continue
                known, inferred, method = _phase_site(
                    chrom, child_gt, parent_gt, parent_role, child_sex, stats
                )
                if method == "skip":
                    continue
                stats.sites += 1
                if known is not None or inferred is not None:
                    stats.per_chrom_resolved[chrom] = stats.per_chrom_resolved.get(chrom, 0) + 1
                if parent_role == "mother":
                    yield chrom, pos, known, inferred, method
                else:
                    yield chrom, pos, inferred, known, method

    if persist:
        store.insert_phase(child_kit, rows())
    else:
        for _ in rows():
            pass
    return stats


def _phase_site(
    chrom: str,
    child_gt: str,
    parent_gt: str,
    parent_role: str,
    child_sex: Optional[str],
    stats: PhaseStats,
) -> Tuple[Optional[str], Optional[str], str]:
    """Resolve one site.  Returns (from_tested_parent, from_other_parent, method)."""
    # A son's single X copy is entirely his mother's; there is no paternal X
    # to infer, since he got his father's Y instead.
    if chrom == "X" and child_sex == "M":
        allele = child_gt[0]
        if parent_role == "mother":
            stats.x_maternal += 1
            return allele, None, "x-hemizygous"
        return None, allele, "x-hemizygous"

    if len(child_gt) != 2:
        return None, None, "skip"

    if parent_gt == NO_CALL or len(parent_gt) != 2:
        stats.parent_missing += 1
        if child_gt[0] == child_gt[1]:
            stats.homozygous += 1
            return child_gt[0], child_gt[1], "hom-no-parent"
        stats.het_ambiguous += 1
        return None, None, "no-parent"

    child_alleles = set(child_gt)
    parent_alleles = set(parent_gt)

    if not (child_alleles & parent_alleles):
        # Opposite homozygotes, or no allele in common at all: the parent
        # could not have transmitted either of the child's alleles.
        stats.mendelian_errors += 1
        return None, None, "mendelian-error"

    if child_gt[0] == child_gt[1]:
        stats.homozygous += 1
        return child_gt[0], child_gt[1], "hom"

    # Child heterozygous.  A homozygous parent pins the transmitted allele.
    if parent_gt[0] == parent_gt[1]:
        from_parent = parent_gt[0]
        other = [a for a in child_gt if a != from_parent]
        if not other:
            stats.mendelian_errors += 1
            return None, None, "mendelian-error"
        stats.het_resolved += 1
        return from_parent, other[0], "het-parent-hom"

    # Both heterozygous: either allele could have come from either side.
    stats.het_ambiguous += 1
    return None, None, "het-het-ambiguous"


# ---------------------------------------------------------------------------
# exporting an inferred parent as a pseudo-kit


def export_phased_kit(
    store: Store,
    child_kit: int,
    side: str,
    out_path: str,
    label: str = "",
) -> Dict[str, int]:
    """Write one phased haplotype as a 23andMe-format pseudo-kit.

    Each resolved allele is written doubled (an inferred ``G`` becomes
    ``GG``) and unresolved sites become no-calls.  That is the convention
    every comparison tool understands, and because only one haplotype is
    present, a comparison against this file cannot produce a
    false compound segment.

    Note what this file is and is not.  The paternal export is half of your
    father's genome -- specifically the half he transmitted to you.  Total
    sharing against it runs about half what the same relative shows against
    your unphased kit, so interpret the numbers on that scale.
    """
    if side not in ("maternal", "paternal"):
        raise ValueError("side must be 'maternal' or 'paternal'")
    column = "mat" if side == "maternal" else "pat"
    counts = {"written": 0, "resolved": 0}
    kit = store.kit(child_kit)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(
            "# Phased pseudo-kit generated by roots (genealogy toolkit)\n"
            f"# source kit: {kit.label if kit else child_kit}\n"
            f"# haplotype: {side}\n"
            "# This is ONE inherited haplotype written as homozygous calls.\n"
            "# It represents half of the corresponding parent's genome.\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
        )
        cur = store.db.execute(
            f"SELECT p.chrom, p.pos, g.rsid, p.{column} AS allele "
            "FROM phased p LEFT JOIN genotype g "
            "  ON g.kit_id=p.kit_id AND g.chrom=p.chrom AND g.pos=p.pos "
            "WHERE p.kit_id=? ORDER BY p.chrom, p.pos",
            (child_kit,),
        )
        for row in cur:
            allele = row["allele"]
            gt = (allele * 2) if allele else NO_CALL
            if allele:
                counts["resolved"] += 1
            rsid = row["rsid"] or f"{row['chrom']}:{row['pos']}"
            fh.write(f"{rsid}\t{row['chrom']}\t{row['pos']}\t{gt}\n")
            counts["written"] += 1
    return counts


def phased_coverage(store: Store, child_kit: int) -> Dict[str, Tuple[int, int]]:
    """Per-chromosome (resolved, total) phased-site counts."""
    out: Dict[str, Tuple[int, int]] = {}
    for row in store.db.execute(
        "SELECT chrom, COUNT(*) AS total, SUM(CASE WHEN pat IS NOT NULL THEN 1 ELSE 0 END) AS res "
        "FROM phased WHERE kit_id=? GROUP BY chrom",
        (child_kit,),
    ):
        out[row["chrom"]] = (row["res"] or 0, row["total"])
    return out

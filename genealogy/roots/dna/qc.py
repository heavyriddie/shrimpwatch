"""Quality control and cross-vendor kit merging."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

from ..genome import ALL_CHROMS, AUTOSOMES
from ..store import Kit, Store
from .raw import NO_CALL


def ibs(gt_a: str, gt_b: str) -> int:
    """Identity-by-state between two diploid genotypes: 2, 1, or 0.

    Returns -1 when either side is a no-call or hemizygous, which callers
    must skip rather than treat as a mismatch.
    """
    if len(gt_a) != 2 or len(gt_b) != 2 or gt_a == NO_CALL or gt_b == NO_CALL:
        return -1
    if gt_a == gt_b:
        return 2
    if set(gt_a) & set(gt_b):
        return 1
    return 0


def is_het(gt: str) -> bool:
    return len(gt) == 2 and gt[0] != gt[1]


def same_call(a: str, b: str) -> bool:
    """Whether two spellings describe the same call.

    Vendors disagree about hemizygous sites: a male X or a mitochondrial
    genotype is written ``T`` by one and ``TT`` by another. Those are the
    same call, and treating them as a conflict would blank real data during
    a merge.

    The parser cannot resolve this itself on the X: a lone ``G`` and a
    ``GG`` are the same call for a man and different spellings of the same
    call for a woman, and sex is not known until the whole file has been
    read. So the ambiguity is settled here, at the point of comparison,
    rather than guessed at parse time.
    """
    if a == b:
        return True
    if NO_CALL in (a, b):
        return False
    return set(a) == set(b) and len(set(a)) == 1


@dataclass
class KitQC:
    kit: Kit
    total: int = 0
    called: int = 0
    het: int = 0
    per_chrom: Dict[str, int] = field(default_factory=dict)
    per_chrom_called: Dict[str, int] = field(default_factory=dict)
    x_het: int = 0
    x_called: int = 0
    y_called: int = 0
    y_total: int = 0
    indel: int = 0

    @property
    def call_rate(self) -> float:
        return self.called / self.total if self.total else 0.0

    @property
    def het_rate(self) -> float:
        return self.het / self.called if self.called else 0.0

    @property
    def inferred_sex(self) -> Optional[str]:
        x_het_rate = self.x_het / self.x_called if self.x_called else 0.0
        y_rate = self.y_called / self.y_total if self.y_total else 0.0
        if self.x_called < 500 and self.y_total == 0:
            return None
        if x_het_rate > 0.10:
            return "F"
        if y_rate > 0.30 or x_het_rate < 0.02:
            return "M"
        return None

    def warnings(self) -> List[str]:
        out: List[str] = []
        if self.call_rate < 0.97:
            out.append(
                f"call rate {self.call_rate:.1%} is low; below ~97% the segment "
                "detector will lose sensitivity and may fragment real segments"
            )
        if self.total < 300_000:
            out.append(
                f"only {self.total:,} SNPs -- fewer than a typical consumer chip "
                "(600k-1.8M); check that the whole file imported"
            )
        # Autosomal heterozygosity outside 0.20-0.40 usually means a sample
        # problem or a non-standard file rather than unusual ancestry.
        if self.called and not (0.15 <= self.het_rate <= 0.45):
            out.append(
                f"heterozygosity {self.het_rate:.1%} is outside the usual range; "
                "possible file corruption, sample mixture, or a non-diploid source"
            )
        missing = [c for c in AUTOSOMES if self.per_chrom.get(c, 0) == 0]
        if missing:
            out.append(f"no data on autosome(s) {', '.join(missing)}")
        return out


def qc_kit(store: Store, kit: Kit) -> KitQC:
    """Recompute per-kit QC from what actually landed in the database."""
    q = KitQC(kit=kit)
    for row in store.db.execute(
        "SELECT chrom, gt, COUNT(*) AS n FROM genotype WHERE kit_id=? GROUP BY chrom, gt",
        (kit.id,),
    ):
        chrom, gt, n = row["chrom"], row["gt"], row["n"]
        q.total += n
        q.per_chrom[chrom] = q.per_chrom.get(chrom, 0) + n
        if chrom == "Y":
            q.y_total += n
        if gt == NO_CALL:
            continue
        q.called += n
        q.per_chrom_called[chrom] = q.per_chrom_called.get(chrom, 0) + n
        if "I" in gt or "D" in gt:
            q.indel += n
        if chrom == "Y":
            q.y_called += n
        if chrom == "X":
            q.x_called += n
            if is_het(gt):
                q.x_het += n
        if is_het(gt):
            q.het += n
    return q


@dataclass
class Concordance:
    shared: int = 0
    identical: int = 0
    ibs1: int = 0
    ibs0: int = 0
    only_a: int = 0
    only_b: int = 0

    @property
    def rate(self) -> float:
        return self.identical / self.shared if self.shared else 0.0

    def verdict(self) -> str:
        """Interpret the concordance rate.

        Two kits from the same person should agree at essentially every
        shared SNP; the residual disagreement is genotyping error, a few
        parts per thousand.  Anything materially worse means the kits are
        not the same person, or the files disagree about strand.
        """
        if self.shared < 5_000:
            return "too few shared SNPs to judge"
        if self.rate >= 0.995 and self.ibs0 <= self.shared * 0.001:
            return "same person"
        if self.rate >= 0.90:
            return (
                "high but imperfect agreement -- same person with poor data, or "
                "strand inconsistency between vendors"
            )
        return "different people"


def concordance(store: Store, kit_a: int, kit_b: int) -> Concordance:
    """Compare two kits SNP by SNP.  Used to validate same-person merges."""
    c = Concordance()
    for chrom in AUTOSOMES:
        for _pos, ga, gb in store.paired_genotypes(kit_a, kit_b, chrom):
            state = ibs(ga, gb)
            if state < 0:
                continue
            c.shared += 1
            if state == 2:
                c.identical += 1
            elif state == 1:
                c.ibs1 += 1
            else:
                c.ibs0 += 1
    return c


def merge_kits(
    store: Store,
    kit_ids: List[int],
    new_label: str,
    person: str,
    build: str = "37",
) -> Tuple[int, Dict[str, int]]:
    """Build a consensus kit from several kits belonging to one person.

    A consensus kit is worth having because vendors type overlapping but
    different SNP sets: combining a 23andMe and an Ancestry download
    typically yields 20-40% more usable positions than either alone, which
    directly improves segment-boundary resolution.

    Conflict policy is deliberately conservative: where two source kits are
    both called and disagree, the merged site becomes a no-call.  Silently
    picking a winner would manufacture false heterozygous sites, which is
    exactly what breaks phasing.
    """
    stats = {"positions": 0, "agree": 0, "conflict": 0, "single_source": 0}
    consensus: Dict[Tuple[str, int], Tuple[Optional[str], str]] = {}
    for kid in kit_ids:
        for row in store.genotypes(kid):
            key = (row["chrom"], row["pos"])
            gt = row["gt"]
            prior = consensus.get(key)
            if prior is None:
                consensus[key] = (row["rsid"], gt)
                continue
            prsid, pgt = prior
            if pgt == gt or same_call(pgt, gt):
                continue
            if pgt == NO_CALL:
                consensus[key] = (prsid or row["rsid"], gt)
            elif gt == NO_CALL:
                continue
            else:
                consensus[key] = (prsid or row["rsid"], NO_CALL)
                stats["conflict"] += 1

    new_id = store.create_kit(new_label, person, "merged", build, "consensus")
    def rows():
        for (chrom, pos), (rsid, gt) in consensus.items():
            yield chrom, pos, rsid, gt

    n = store.insert_genotypes(new_id, sorted(rows(), key=lambda r: (r[0], r[1])))
    stats["positions"] = n
    q = qc_kit(store, store.kit(new_id))
    store.finalize_kit(new_id, q.total, q.called, q.inferred_sex)
    return new_id, stats

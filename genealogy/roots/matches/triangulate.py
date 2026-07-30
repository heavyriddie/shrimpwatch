"""Group matches by the chromosome region they share with you.

If several matches all share the same stretch of your chromosome 7, that
stretch probably came from one ancestral couple, and identifying the common
ancestor of any two of them identifies the source of the segment for all of
them.  This is the most direct route from anonymous matches to a named
ancestor, and unlike clustering it tells you *where* the DNA came from, not
just that it travelled together.

Two cautions are built into the output rather than left to the reader.

First, overlapping is not triangulating.  You have two copies of chromosome
7, one from each parent.  Two matches can share the same coordinates with
you while sitting on opposite copies, descending from entirely unrelated
families.  Genuine triangulation requires that the two matches also match
*each other* across that region.  Groups are labelled by which of those two
standards they actually meet, and a group whose members split across
maternal and paternal is called out as the two-copies artefact it is.

Second, a shared region is only as trustworthy as its size.  Below about
7 cM most reported segments are not recent inheritance, so a large group of
matches over a short region is usually an artefact of the region itself --
some parts of the genome are simply prone to spurious matching -- rather
than a big undiscovered branch.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

from ..genome import GeneticMap
from ..store import Store


@dataclass
class OverlapGroup:
    chrom: str
    start_bp: int
    end_bp: int
    cm: float
    match_ids: List[int] = field(default_factory=list)
    names: List[str] = field(default_factory=list)
    sides: List[Optional[str]] = field(default_factory=list)
    support: str = "overlap only"
    icw_fraction: float = 0.0
    conflict: bool = False

    @property
    def size(self) -> int:
        return len(self.match_ids)

    def side(self) -> Optional[str]:
        known = [s for s in self.sides if s in ("maternal", "paternal")]
        if not known:
            return None
        if all(s == known[0] for s in known):
            return known[0]
        return "conflicting"

    def describe(self) -> str:
        loc = f"chr{self.chrom}:{self.start_bp / 1e6:.1f}-{self.end_bp / 1e6:.1f} Mb"
        return (
            f"{loc} ({self.cm:.1f} cM), {self.size} matches, {self.support}"
            + (f", side {self.side()}" if self.side() else "")
        )


def overlap_groups(
    store: Store,
    kit_id: int,
    gmap: GeneticMap,
    min_members: int = 2,
    min_overlap_cm: float = 7.0,
    source: Optional[str] = None,
    max_match_cm: float = 2000.0,
) -> List[OverlapGroup]:
    """Find regions where several matches share DNA with the same kit.

    Uses a sweep over segment boundaries, so a group is a maximal region
    covered by exactly the same set of matches -- not a chain of segments
    that merely touch end to end.

    ``max_match_cm`` drops relatives who are too close to be informative.
    A parent matches you across every chromosome, and a full sibling across
    three quarters of one, so they join every group that exists and tell you
    nothing about where a segment came from. Excluding them is standard
    practice, not a shortcut: triangulation asks which ancestor a segment
    came from, and a relative who shares everything cannot narrow that down.
    Second-degree relatives and further are kept, since an aunt sharing a
    specific region genuinely locates it.
    """
    match_rows = {
        m["id"]: m for m in store.matches(kit_id, source=source)
        if (m["total_cm"] or 0.0) <= max_match_cm
    }
    by_chrom: Dict[str, List[Tuple[int, int, int]]] = {}
    for seg in store.match_segments():
        mid = seg["match_id"]
        if mid not in match_rows:
            continue
        by_chrom.setdefault(seg["chrom"], []).append((seg["start_bp"], seg["end_bp"], mid))

    icw = _icw_pairs(store, kit_id)
    triangulated = _icw_pairs(store, kit_id, kind="triangulated")

    groups: List[OverlapGroup] = []
    for chrom, segs in by_chrom.items():
        points: Set[int] = set()
        for start, end, _mid in segs:
            points.add(start)
            points.add(end)
        edges = sorted(points)
        current: Optional[Tuple[int, frozenset]] = None
        runs: List[Tuple[int, int, frozenset]] = []
        for i in range(len(edges) - 1):
            lo, hi = edges[i], edges[i + 1]
            if hi <= lo:
                continue
            covering = frozenset(
                mid for start, end, mid in segs if start <= lo and end >= hi
            )
            if current and current[1] == covering:
                continue
            if current:
                runs.append((current[0], lo, current[1]))
            current = (lo, covering)
        if current:
            runs.append((current[0], edges[-1], current[1]))

        for start, end, members in runs:
            if len(members) < min_members:
                continue
            cm = gmap.length_cm(chrom, start, end)
            if cm < min_overlap_cm:
                continue
            member_list = sorted(members)
            g = OverlapGroup(
                chrom=chrom, start_bp=start, end_bp=end, cm=cm,
                match_ids=member_list,
                names=[match_rows[m]["name"] or str(m) for m in member_list],
                sides=[match_rows[m]["side"] for m in member_list],
            )
            g.icw_fraction, g.support = _support(member_list, icw, triangulated)
            g.conflict = g.side() == "conflicting"
            groups.append(g)

    groups.sort(key=lambda g: (-g.size, -g.cm))
    return groups


def _icw_pairs(store: Store, kit_id: int, kind: str = "icw") -> Set[Tuple[int, int]]:
    out: Set[Tuple[int, int]] = set()
    for row in store.shared_matches(kit_id, kind=kind):
        out.add((min(row["a_id"], row["b_id"]), max(row["a_id"], row["b_id"])))
    return out


def _support(
    members: Sequence[int],
    icw: Set[Tuple[int, int]],
    triangulated: Set[Tuple[int, int]],
) -> Tuple[float, str]:
    pairs = [
        (min(a, b), max(a, b))
        for i, a in enumerate(members)
        for b in members[i + 1:]
    ]
    if not pairs:
        return 0.0, "overlap only"
    tri = sum(1 for p in pairs if p in triangulated)
    common = sum(1 for p in pairs if p in icw or p in triangulated)
    frac = common / len(pairs)
    if tri == len(pairs):
        return frac, "triangulated"
    if tri:
        return frac, "partly triangulated"
    if frac > 0.8:
        return frac, "shared-match supported"
    if frac > 0:
        return frac, "partly shared-match supported"
    return 0.0, "overlap only"


# ---------------------------------------------------------------------------
# chromosome painting


@dataclass
class PaintedSegment:
    chrom: str
    start_bp: int
    end_bp: int
    cm: float
    side: Optional[str]
    match_id: int
    name: str


def paint(
    store: Store,
    kit_id: int,
    gmap: GeneticMap,
    min_cm: float = 7.0,
    source: Optional[str] = None,
) -> Dict[str, List[PaintedSegment]]:
    """Lay every match's segments out along the chromosomes.

    The result is the data behind a chromosome map: which stretches of your
    genome you can currently attribute to a side, and -- more usefully --
    which stretches you cannot attribute to anyone yet.  Those gaps are where
    your unexplored ancestry lives.
    """
    match_rows = {m["id"]: m for m in store.matches(kit_id, source=source)}
    out: Dict[str, List[PaintedSegment]] = {}
    for seg in store.match_segments():
        m = match_rows.get(seg["match_id"])
        if not m:
            continue
        cm = seg["cm"] if seg["cm"] is not None else gmap.length_cm(
            seg["chrom"], seg["start_bp"], seg["end_bp"]
        )
        if cm < min_cm:
            continue
        out.setdefault(seg["chrom"], []).append(
            PaintedSegment(
                chrom=seg["chrom"],
                start_bp=seg["start_bp"],
                end_bp=seg["end_bp"],
                cm=cm,
                side=m["side"],
                match_id=m["id"],
                name=m["name"] or str(m["remote_id"]),
            )
        )
    for segs in out.values():
        segs.sort(key=lambda s: (s.start_bp, s.end_bp))
    return out


def coverage(
    painted: Dict[str, List[PaintedSegment]], gmap: GeneticMap
) -> Dict[str, float]:
    """How much of the genome is attributed to each side, in centimorgans.

    Coverage counts each stretch once however many matches cover it.  The
    ``unattributed`` figure is the honest headline: it is the portion of your
    DNA that no current match explains.
    """
    from ..genome import AUTOSOMES

    totals = {"maternal": 0.0, "paternal": 0.0, "any": 0.0}
    for chrom, segs in painted.items():
        for side in ("maternal", "paternal", "any"):
            spans = [
                (s.start_bp, s.end_bp)
                for s in segs
                if side == "any" or s.side == side
            ]
            for start, end in _merge_spans(spans):
                totals[side] += gmap.length_cm(chrom, start, end)
    genome = sum(gmap.chrom_cm(c) for c in AUTOSOMES)
    totals["genome"] = genome
    totals["unattributed"] = max(genome - totals["any"], 0.0)
    return totals


def _merge_spans(spans: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []
    for start, end in sorted(spans):
        if out and start <= out[-1][1]:
            if end > out[-1][1]:
                out[-1] = (out[-1][0], end)
        else:
            out.append((start, end))
    return out

"""Sort DNA matches onto the maternal or paternal side.

This is the biggest single payoff from having a parent tested, and it is why
a mother's kit is worth more than any amount of extra analysis on your own.

The rule is almost embarrassingly simple.  Every one of your matches is
related to you through your mother or through your father (or, rarely, both).
If a match also appears on your mother's match list, they are maternal.  If
they do not appear on her list, and they share enough DNA with you that they
would certainly have appeared had they been maternal, they are paternal.
That second half is the powerful one: it identifies your father's relatives
using only your mother's data, which is exactly the situation when the
paternal line is the unknown one.

Three conditions have to hold for the inference to be sound, and this module
enforces all three rather than assuming them:

1.  Both lists must come from the same testing platform.  A person who
    tested at 23andMe will not appear on an Ancestry list no matter how
    closely related they are, so a cross-platform absence means nothing.

2.  The match must share more DNA with you than the platform's reporting
    floor, with margin.  Sites truncate their match lists; a 9 cM match
    missing from your mother's list may simply have fallen off the end of
    it.  We take the smallest total on the parent's list as the observed
    floor and refuse to infer "paternal" below a configurable multiple of it.

3.  The pair should not be shared through both parents.  If your parents are
    themselves related -- which is the norm in endogamous populations and
    common in any isolated community -- a match can be genuinely maternal
    and paternal at once.  We flag the pattern rather than forcing a side.

The quantitative cross-check is the ratio.  A relative on your mother's side
sits one generation further from you than from her, so they should share
roughly twice as much DNA with her as with you.  A match who appears on her
list but shares *less* with her than with you is not behaving like a
maternal relative, and gets flagged.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from ..store import Store
from .ingest import normalize_name

#: Sharing at or above this is whole-genome: a parent, a child, or an
#: identical twin. Full siblings top out well below it, around 2900 cM.
PARENT_LEVEL_CM = 3000.0

#: Above this, two children of the same mother also share a father. Full
#: siblings average about 2613 cM and half siblings about 1759 cM, so the
#: midpoint separates them with room to spare on both sides.
FULL_SIBLING_FLOOR_CM = 2200.0


@dataclass
class SideCall:
    match_id: int
    name: str
    side: Optional[str]
    confidence: str
    reason: str
    self_cm: float
    parent_cm: Optional[float] = None


@dataclass
class SideReport:
    parent_role: str
    source: str
    calls: List[SideCall] = field(default_factory=list)
    parent_floor_cm: float = 0.0
    warnings: List[str] = field(default_factory=list)

    def counts(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for c in self.calls:
            key = c.side or "unknown"
            out[key] = out.get(key, 0) + 1
        return out

    def flagged(self) -> List[SideCall]:
        return [c for c in self.calls if c.side == "both" or "flag" in c.confidence]

    def summary(self) -> str:
        counts = self.counts()
        parts = [f"{v} {k}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1])]
        return ", ".join(parts) if parts else "no matches processed"


def assign_from_parent_list(
    store: Store,
    child_kit: int,
    parent_kit: int,
    parent_role: str = "mother",
    source: Optional[str] = None,
    floor_multiple: float = 2.0,
    min_paternal_cm: float = 20.0,
    persist: bool = True,
) -> SideReport:
    """Assign sides by comparing a child's match list with a parent's.

    ``floor_multiple`` sets how far above the parent list's observed
    truncation point a match must sit before its *absence* from that list is
    treated as evidence.  ``min_paternal_cm`` is a hard floor on the same
    logic.
    """
    other_side = "paternal" if parent_role == "mother" else "maternal"
    same_side = "maternal" if parent_role == "mother" else "paternal"

    child_matches = store.matches(child_kit)
    parent_matches = store.matches(parent_kit)
    sources_child = {m["source"] for m in child_matches}
    sources_parent = {m["source"] for m in parent_matches}
    usable_sources = sources_child & sources_parent
    if source:
        usable_sources = usable_sources & {source}

    report = SideReport(parent_role=parent_role, source=", ".join(sorted(usable_sources)))
    if not usable_sources:
        report.warnings.append(
            f"no platform appears in both match lists (child: "
            f"{sorted(sources_child) or 'none'}; {parent_role}: "
            f"{sorted(sources_parent) or 'none'}). Absence from the "
            f"{parent_role}'s list cannot be interpreted across platforms, so "
            "no sides were inferred."
        )
        return report

    parent_index: Dict[Tuple[str, str], float] = {}
    for m in parent_matches:
        if m["source"] not in usable_sources:
            continue
        cm = m["total_cm"] or 0.0
        for key in _keys(m):
            prior = parent_index.get((m["source"], key))
            if prior is None or cm > prior:
                parent_index[(m["source"], key)] = cm

    floors = [m["total_cm"] for m in parent_matches
              if m["source"] in usable_sources and m["total_cm"]]
    report.parent_floor_cm = min(floors) if floors else 0.0
    inference_floor = max(report.parent_floor_cm * floor_multiple, min_paternal_cm)
    if report.parent_floor_cm > 50:
        report.warnings.append(
            f"the {parent_role}'s list only goes down to "
            f"{report.parent_floor_cm:.0f} cM, which is unusually high -- it "
            "looks truncated. Absence below "
            f"{inference_floor:.0f} cM will not be used as evidence."
        )

    for m in child_matches:
        if m["source"] not in usable_sources:
            continue
        self_cm = m["total_cm"] or 0.0
        parent_cm = None
        for key in _keys(m):
            hit = parent_index.get((m["source"], key))
            if hit is not None:
                parent_cm = hit if parent_cm is None else max(parent_cm, hit)

        call = _decide(
            m, self_cm, parent_cm, parent_role, same_side, other_side,
            inference_floor,
        )
        report.calls.append(call)
        if persist and call.side:
            store.set_match_side(call.match_id, call.side, call.reason)
    if persist:
        store.commit()
    return report


def _keys(match_row) -> List[str]:
    out = []
    if match_row["name"]:
        out.append(normalize_name(match_row["name"]))
    if match_row["remote_id"]:
        out.append(normalize_name(str(match_row["remote_id"])))
    return [k for k in out if k]


def _decide(
    m,
    self_cm: float,
    parent_cm: Optional[float],
    parent_role: str,
    same_side: str,
    other_side: str,
    inference_floor: float,
) -> SideCall:
    name = m["name"] or str(m["remote_id"])

    if parent_cm is None and self_cm >= PARENT_LEVEL_CM:
        # Nobody appears on their own match list, so the tested parent shows
        # up here as "shares your whole genome, absent from the list" -- which
        # the absence rule would otherwise read as paternal. Anyone at this
        # level is a parent, a child, or an identical twin, and none of them
        # belong to one side.
        return SideCall(
            m["id"], name, None, "not applicable",
            f"shares {self_cm:.0f} cM, essentially your whole genome, so this is "
            f"a parent, a child, or an identical twin -- most likely the tested "
            f"{parent_role} herself, who cannot appear on her own match list. "
            "Sides do not apply to them",
            self_cm, None,
        )

    if parent_cm is not None and parent_cm >= PARENT_LEVEL_CM:
        # The parent shares their whole genome with this person, so the person
        # is that parent's child: your sibling. How much you share with them
        # then says whether you share a father as well.
        if self_cm >= FULL_SIBLING_FLOOR_CM:
            return SideCall(
                m["id"], name, "both", "high",
                f"is your {parent_role}'s child ({parent_cm:.0f} cM with her) and "
                f"shares {self_cm:.0f} cM with you, which is full-sibling range -- "
                "so you share a father too, and they are related to you on both "
                "sides",
                self_cm, parent_cm,
            )
        return SideCall(
            m["id"], name, same_side, "high",
            f"is your {parent_role}'s child ({parent_cm:.0f} cM with her) but "
            f"shares only {self_cm:.0f} cM with you -- half-sibling range, so you "
            f"share a {parent_role} and not a father",
            self_cm, parent_cm,
        )

    if parent_cm is not None:
        # Present on the parent's list.  Expected ratio is about 2:1 in the
        # parent's favour, since they sit one generation closer.
        ratio = parent_cm / self_cm if self_cm else 0.0
        if ratio >= 0.8:
            return SideCall(
                m["id"], name, same_side, "high",
                f"appears on the {parent_role}'s match list at {parent_cm:.0f} cM "
                f"versus {self_cm:.0f} cM with you (ratio {ratio:.1f}, as expected "
                f"for a {same_side} relative)",
                self_cm, parent_cm,
            )
        return SideCall(
            m["id"], name, "both", "flagged",
            f"shares {self_cm:.0f} cM with you but only {parent_cm:.0f} cM with "
            f"your {parent_role} (ratio {ratio:.1f}). A {same_side} relative "
            "should share more with her than with you, so this pair is likely "
            "related through both parents, or one of the two segments is a "
            "false positive",
            self_cm, parent_cm,
        )

    if self_cm >= inference_floor:
        return SideCall(
            m["id"], name, other_side, "high",
            f"shares {self_cm:.0f} cM with you and is absent from the "
            f"{parent_role}'s list, which reaches well below that level -- so "
            f"the relationship runs through your {other_side} line",
            self_cm, None,
        )

    return SideCall(
        m["id"], name, None, "low",
        f"only {self_cm:.0f} cM shared, below the {inference_floor:.0f} cM point "
        f"where absence from the {parent_role}'s truncated list is informative",
        self_cm, None,
    )


# ---------------------------------------------------------------------------
# segment-level confirmation


def confirm_with_segments(
    store: Store,
    child_kit: int,
    parent_kit: int,
    parent_role: str = "mother",
    min_overlap_bp: int = 1_000_000,
) -> List[SideCall]:
    """Confirm side calls using segment overlap, where segment data exists.

    Stronger evidence than list membership alone: if the same person shares a
    segment with you at chr7:20-45 Mb and shares an overlapping segment with
    your mother, that specific segment is maternal.  It also catches the case
    where someone matches you on two segments, one from each side.
    """
    same_side = "maternal" if parent_role == "mother" else "paternal"
    child_segs = _segments_by_name(store, child_kit)
    parent_segs = _segments_by_name(store, parent_kit)
    calls: List[SideCall] = []
    for key, (mid, name, segs) in child_segs.items():
        parent_entry = parent_segs.get(key)
        if not parent_entry:
            continue
        _pmid, _pname, psegs = parent_entry
        confirmed = 0.0
        total = 0.0
        for chrom, start, end, cm in segs:
            total += cm or 0.0
            for pchrom, pstart, pend, _pcm in psegs:
                if pchrom != chrom:
                    continue
                overlap = min(end, pend) - max(start, pstart)
                if overlap >= min_overlap_bp:
                    confirmed += cm or 0.0
                    break
        if total <= 0:
            continue
        frac = confirmed / total
        if frac > 0.9:
            reason = (
                f"every shared segment overlaps a segment this person also "
                f"shares with your {parent_role}"
            )
            side = same_side
        elif frac > 0.05:
            reason = (
                f"{frac:.0%} of shared segments overlap segments shared with "
                f"your {parent_role}; the remainder come from your other side, "
                "so this person is related on both"
            )
            side = "both"
        else:
            continue
        calls.append(SideCall(mid, name, side, "segment-confirmed", reason, total, confirmed))
    return calls


def _segments_by_name(store: Store, kit_id: int):
    out: Dict[str, Tuple[int, str, List[Tuple[str, int, int, Optional[float]]]]] = {}
    for m in store.matches(kit_id):
        segs = [
            (s["chrom"], s["start_bp"], s["end_bp"], s["cm"])
            for s in store.match_segments(m["id"])
        ]
        if segs:
            out[normalize_name(m["name"] or "")] = (m["id"], m["name"] or "", segs)
    return out


def side_summary(store: Store, kit_id: int) -> Dict[str, Dict[str, float]]:
    """Per-side counts and total centimorgans, for a report header."""
    out: Dict[str, Dict[str, float]] = {}
    for m in store.matches(kit_id):
        key = m["side"] or "unassigned"
        bucket = out.setdefault(key, {"matches": 0, "total_cm": 0.0, "largest_cm": 0.0})
        bucket["matches"] += 1
        bucket["total_cm"] += m["total_cm"] or 0.0
        bucket["largest_cm"] = max(bucket["largest_cm"], m["total_cm"] or 0.0)
    return out

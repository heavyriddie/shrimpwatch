"""Documentary evidence: sources, citations, and what to look up next.

DNA tells you that a connection exists and roughly how close it is. It never
tells you a name, a date, or a place. Those come from records -- civil
registration, censuses, parish registers, wills -- and the whole point of the
DNA work is to tell you *which* records are worth paying for.

This module is the bridge. It records what documents you have actually seen
and what each of them asserts, and it generates a prioritised list of records
to look up next, derived from where the tree runs out and where the DNA says
the answers must be.

Three ideas shape it.

**A source and a claim are different things.** The individual table holds
your current best understanding of a person. A source holds a document. An
evidence row holds what that document says about that person -- which may
contradict what another document says. Genealogy without that separation
degenerates into a tree of unattributed assertions, and the entire value of
having DNA is that it can adjudicate between conflicting documents.

**The research frontier is where the parents are unknown.** Every ancestor
in your tree with no recorded parents is a wall, and the record that most
often breaks it is not a birth record but a *marriage* record, because in
England and Wales a marriage certificate names the fathers of both parties.
One certificate can therefore advance two lines at once. That single fact
determines most of the ordering this module produces.

**Where you are matters more than when.** Scottish statutory records are far
richer than English ones -- a Scottish death certificate names both parents
of the deceased, so a death record alone can extend a line, which is never
true in England. Irish civil registration starts later and most pre-1901
censuses were destroyed. The suggestions below are jurisdiction-aware because
following English assumptions in Scotland wastes money.

Costs are described in bands rather than figures, because prices change and a
stale number in code is worse than no number. See docs/RECORD_SOURCES.md for
current pricing.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

from .store import Store
from .tree import kinship

#: Cost bands, deliberately coarse. "low" means the price of a coffee or two
#: per record; "subscription" means you need a paid site but the marginal
#: record is then free.
COST_BANDS = ("free", "low", "subscription", "unknown")

# Jurisdiction milestones that decide which record set can answer a question.
CIVIL_REGISTRATION_START = {
    "england": 1837, "wales": 1837, "england and wales": 1837,
    "scotland": 1855,
    "ireland": 1864, "northern ireland": 1864,
}

CENSUS_YEARS = {
    "england": (1841, 1851, 1861, 1871, 1881, 1891, 1901, 1911, 1921),
    "wales": (1841, 1851, 1861, 1871, 1881, 1891, 1901, 1911, 1921),
    "scotland": (1841, 1851, 1861, 1871, 1881, 1891, 1901, 1911, 1921),
    "ireland": (1901, 1911),
}

_COUNTRY_HINTS = (
    ("scotland", ("scotland", "scottish", "aberdeen", "glasgow", "edinburgh",
                  "lanark", "ayrshire", "fife", "dundee", "inverness")),
    ("ireland", ("ireland", "irish", "dublin", "cork", "galway", "mayo",
                 "kerry", "tipperary", "limerick")),
    ("northern ireland", ("antrim", "armagh", "belfast", "derry", "londonderry",
                          "down", "fermanagh", "tyrone")),
    ("wales", ("wales", "welsh", "glamorgan", "cardiff", "swansea", "gwynedd",
               "denbigh", "carmarthen", "pembroke")),
    ("england", ("england", "london", "yorkshire", "lancashire", "kent",
                 "devon", "essex", "norfolk", "durham", "surrey", "sussex",
                 "leeds", "manchester", "liverpool", "birmingham")),
)


def infer_country(place: Optional[str]) -> Optional[str]:
    """Guess a jurisdiction from a place string.

    Crude by design: it only has to be right often enough to pick the correct
    record set, and it returns None rather than guessing when nothing matches.
    """
    if not place:
        return None
    p = place.lower()
    for country, needles in _COUNTRY_HINTS:
        for needle in needles:
            if needle in p:
                return country
    return None


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ---------------------------------------------------------------------------
# recording what you have


@dataclass
class Source:
    id: Optional[int]
    title: str
    repository: Optional[str] = None
    record_type: Optional[str] = None
    reference: Optional[str] = None
    url: Optional[str] = None
    accessed: Optional[str] = None
    cost: Optional[float] = None
    notes: Optional[str] = None


def add_source(store: Store, title: str, **fields: Any) -> int:
    cols = ["repository", "record_type", "reference", "url", "accessed", "cost", "notes"]
    payload = {c: fields.get(c) for c in cols}
    payload["accessed"] = payload["accessed"] or _now()[:10]
    keys = ["title"] + cols
    vals = [title] + [payload[c] for c in cols]
    cur = store.db.execute(
        f"INSERT INTO source({', '.join(keys)}) VALUES({', '.join('?' * len(keys))})",
        vals,
    )
    store.commit()
    return int(cur.lastrowid)


def add_evidence(
    store: Store,
    source_id: int,
    subject_xref: str,
    claim_type: str,
    claim: str,
    supports: bool = True,
    confidence: str = "direct",
    notes: Optional[str] = None,
) -> int:
    """Record what a document asserts about a person.

    ``supports=False`` is the interesting case: it records a document that
    contradicts the tree as it stands. Those are kept rather than resolved,
    because a contradiction between two sources is exactly the situation DNA
    is qualified to settle.
    """
    cur = store.db.execute(
        "INSERT INTO evidence(source_id, subject_xref, claim_type, claim, "
        "supports, confidence, notes) VALUES(?,?,?,?,?,?,?)",
        (source_id, subject_xref, claim_type, claim, 1 if supports else 0,
         confidence, notes),
    )
    store.commit()
    return int(cur.lastrowid)


def sources(store: Store) -> List[Any]:
    return list(store.db.execute("SELECT * FROM source ORDER BY id"))


def evidence_for(store: Store, xref: str) -> List[Any]:
    return list(
        store.db.execute(
            "SELECT e.*, s.title, s.repository, s.record_type FROM evidence e "
            "LEFT JOIN source s ON s.id=e.source_id WHERE e.subject_xref=? "
            "ORDER BY e.id",
            (xref,),
        )
    )


def contradictions(store: Store) -> List[Any]:
    return list(
        store.db.execute(
            "SELECT e.*, s.title, s.repository FROM evidence e "
            "LEFT JOIN source s ON s.id=e.source_id WHERE e.supports=0 "
            "ORDER BY e.subject_xref"
        )
    )


@dataclass
class EvidenceCoverage:
    individuals: int = 0
    with_evidence: int = 0
    sources: int = 0
    total_cost: float = 0.0
    unsourced_ancestors: List[Tuple[str, str]] = field(default_factory=list)

    @property
    def fraction(self) -> float:
        return self.with_evidence / self.individuals if self.individuals else 0.0

    def summary(self) -> str:
        return (
            f"{self.with_evidence} of {self.individuals} people in the tree "
            f"({self.fraction:.0%}) have at least one document behind them, "
            f"from {self.sources} sources"
        )


def coverage(store: Store, idx: Optional[kinship.TreeIndex] = None) -> EvidenceCoverage:
    """How much of the tree rests on documents rather than assertion.

    Usually a sobering number on a tree assembled from other people's online
    trees, and worth knowing before DNA results are used to "confirm" any of
    it.
    """
    idx = idx or kinship.load_tree(store)
    cov = EvidenceCoverage(individuals=len(idx.individuals))
    rows = store.db.execute("SELECT DISTINCT subject_xref FROM evidence").fetchall()
    sourced = {r["subject_xref"] for r in rows if r["subject_xref"]}
    cov.with_evidence = len(sourced & set(idx.individuals))
    row = store.db.execute(
        "SELECT COUNT(*) AS n, COALESCE(SUM(cost), 0) AS c FROM source"
    ).fetchone()
    cov.sources = row["n"]
    cov.total_cost = row["c"] or 0.0
    return cov


# ---------------------------------------------------------------------------
# working out what to look up next


@dataclass
class ResearchTask:
    subject: str
    subject_xref: Optional[str]
    question: str
    record_set: str
    repository: str
    cost_band: str = "unknown"
    priority: float = 0.5
    rationale: str = ""

    def line(self) -> str:
        return (
            f"[{self.priority:.2f}] {self.subject}: {self.question}\n"
            f"        look in: {self.record_set} ({self.repository}, {self.cost_band})"
        )


def frontier(
    idx: kinship.TreeIndex, root: str, max_generations: int = 12
) -> List[Tuple[str, int]]:
    """Ancestors with no recorded parents: the edge of the documented tree.

    Returned with their generation depth, closest first, because a wall two
    generations up blocks far more of the tree than one eight generations up
    and is nearly always cheaper to break.
    """
    depths = kinship.ancestors(idx, root, max_generations)
    out: List[Tuple[str, int]] = []
    for xref, depth in depths.items():
        if not idx.parents.get(xref):
            out.append((xref, depth))
    out.sort(key=lambda t: t[1])
    return out


GENERATION_YEARS = 30


def estimate_birth_year(
    idx: kinship.TreeIndex, xref: str, max_depth: int = 6
) -> Tuple[Optional[int], bool]:
    """A person's birth year, estimated from descendants when not recorded.

    Returns ``(year, is_estimated)``.  Wall ancestors are precisely the people
    with no dates -- that is usually why they are walls -- so refusing to
    suggest anything without a recorded year would stay silent exactly where
    help is most needed.  Counting back a generation at a time from the
    nearest dated descendant is the same arithmetic a researcher does in
    their head, and it is good enough to pick the right record set and the
    right decade of an index.
    """
    ind = idx.individuals.get(xref, {})
    recorded = ind.get("birth_year")
    if recorded:
        return recorded, False

    frontier_nodes: List[Tuple[str, int]] = [(xref, 0)]
    seen: Set[str] = {xref}
    while frontier_nodes:
        current, depth = frontier_nodes.pop(0)
        if depth >= max_depth:
            continue
        for fam in idx.families:
            fam_row = idx.families[fam]
            if fam_row.get("husb") != current and fam_row.get("wife") != current:
                continue
            for child, _rel in idx.children.get(fam, []):
                if child in seen:
                    continue
                seen.add(child)
                year = (idx.individuals.get(child, {}) or {}).get("birth_year")
                if year:
                    return year - GENERATION_YEARS * (depth + 1), True
                frontier_nodes.append((child, depth + 1))
    return None, False


def _person_facts(
    idx: kinship.TreeIndex, xref: str
) -> Tuple[Optional[int], Optional[str], Optional[str], bool]:
    ind = idx.individuals.get(xref, {})
    year, estimated = estimate_birth_year(idx, xref)
    place = ind.get("birth_place") or ind.get("death_place")
    if not place:
        place = _inherited_place(idx, xref)
    return year, place, ind.get("surname"), estimated


def _inherited_place(idx: kinship.TreeIndex, xref: str) -> Optional[str]:
    """Borrow a place from a child, since families rarely move far in one step."""
    for fam, fam_row in idx.families.items():
        if fam_row.get("husb") != xref and fam_row.get("wife") != xref:
            continue
        for child, _rel in idx.children.get(fam, []):
            ind = idx.individuals.get(child, {}) or {}
            place = ind.get("birth_place") or ind.get("death_place")
            if place:
                return place
    return None


def suggest_for_person(
    idx: kinship.TreeIndex,
    xref: str,
    depth: int,
    side: Optional[str] = None,
) -> List[ResearchTask]:
    """Records that could identify one wall ancestor's parents.

    Ordered by yield per pound, which in England and Wales means the
    marriage certificate first: it names both fathers, so it advances two
    lines for one purchase, and it is usually easier to locate than a birth
    because the couple's names narrow the index sharply.
    """
    birth_year, place, surname, estimated = _person_facts(idx, xref)
    country = infer_country(place) or "england"
    name = idx.name(xref)
    about = "estimated " if estimated else ""
    tasks: List[ResearchTask] = []

    if birth_year is None:
        # Nothing at all is known about when this person lived, so the only
        # move is to work from a relative who is dated.
        return [ResearchTask(
            name, xref,
            "when and where did they live? Nothing is dated for them yet",
            "the birth or marriage record of their known child, which names "
            "the parents and fixes the generation",
            "GRO / FreeBMD / FamilySearch", "low", max(0.35, 0.90 - depth * 0.08),
            "an undated ancestor cannot be searched for directly; the way in "
            "is always through a descendant whose record names them",
        )]

    # Closer walls block more of the tree, and a known side that matters to
    # the investigation raises the value further.
    base = max(0.30, 0.95 - depth * 0.08)
    if side == "paternal":
        base = min(0.99, base + 0.05)

    civil_start = CIVIL_REGISTRATION_START.get(country, 1837)
    marriage_era = (birth_year + 22) if birth_year else None

    if country == "scotland":
        tasks.append(ResearchTask(
            name, xref,
            "who were their parents?",
            "statutory death register (Scotland) -- Scottish death entries name "
            "both parents of the deceased, including the mother's maiden name",
            "ScotlandsPeople", "low", base,
            "in Scotland a death record alone can extend a line, which is never "
            "true in England; it is usually the cheapest way up a generation",
        ))
        if marriage_era and marriage_era >= 1855:
            tasks.append(ResearchTask(
                name, xref, "who were their parents?",
                f"statutory marriage register (Scotland), around {about}{marriage_era}",
                "ScotlandsPeople", "low", base - 0.03,
                "Scottish marriage entries name both sets of parents",
            ))
    elif country in ("ireland", "northern ireland"):
        if marriage_era and marriage_era >= civil_start:
            tasks.append(ResearchTask(
                name, xref, "who was their father?",
                f"civil marriage register, around {about}{marriage_era}",
                "irishgenealogy.ie / GRONI", "free", base,
                "Irish civil registration begins in 1864 (1845 for non-Catholic "
                "marriages); the free state index often carries register images",
            ))
        tasks.append(ResearchTask(
            name, xref, "where was the household?",
            "1901 and 1911 censuses of Ireland",
            "National Archives of Ireland", "free", base - 0.10,
            "the surviving Irish censuses are free and name every household "
            "member, which is how most Irish lines are reconstructed given the "
            "destruction of the earlier returns",
        ))
    else:
        if marriage_era and marriage_era >= civil_start:
            tasks.append(ResearchTask(
                name, xref, "who were their fathers?",
                f"GRO marriage index and certificate, around {about}{marriage_era}",
                "GRO / FreeBMD", "low", base,
                "an England and Wales marriage certificate names the fathers of "
                "both bride and groom, so one purchase advances two lines",
            ))
        if birth_year and birth_year >= civil_start:
            tasks.append(ResearchTask(
                name, xref, "who was their mother?",
                f"GRO birth index and PDF, around {about}{birth_year}",
                "GRO / FreeBMD", "low", base - 0.05,
                "the GRO birth index carries the mother's maiden surname from "
                "1837 onward, which often identifies the parents without buying "
                "anything at all",
            ))
        if birth_year and birth_year < civil_start:
            tasks.append(ResearchTask(
                name, xref, "who were their parents?",
                f"parish baptism registers, around {about}{birth_year}, "
                f"{place or 'place unknown'}",
                "FamilySearch / FindMyPast / county record office", "subscription",
                base - 0.05,
                "before civil registration the parish register is the only "
                "systematic record of a birth; FamilySearch holds much of it free",
            ))

    census_years = CENSUS_YEARS.get(country, CENSUS_YEARS["england"])
    if birth_year:
        candidates = [y for y in census_years if birth_year <= y <= birth_year + 70]
        if candidates:
            first = candidates[0]
            tasks.append(ResearchTask(
                name, xref,
                "who else was in the household, and where?",
                f"{first} census" + (f" (and {candidates[-1]})" if len(candidates) > 1 else ""),
                "FindMyPast / Ancestry / FamilySearch", "subscription",
                base - 0.12,
                "a census entry places the family, gives ages and birthplaces "
                "for everyone under the roof, and frequently names a widowed "
                "parent living with them",
            ))

    if birth_year and birth_year < 1910:
        tasks.append(ResearchTask(
            name, xref, "when and where did they die, and what did they leave?",
            "probate calendar (England and Wales, 1858 onward)",
            "gov.uk Find a Will", "low", base - 0.20,
            "the probate calendar is indexed by name and year and names "
            "executors, who are very often the children",
        ))
    return tasks


def generate_tasks(
    store: Store,
    root: Optional[str] = None,
    idx: Optional[kinship.TreeIndex] = None,
    kit_id: Optional[int] = None,
    max_generations: int = 12,
    limit: int = 40,
    persist: bool = True,
) -> List[ResearchTask]:
    """Build a prioritised research plan from the tree and the DNA.

    Two sources of work. Walls in the documented tree produce record
    lookups. Clusters of DNA matches with no tree connection produce a
    different kind of task -- there is no record to order, the work is to
    build the matches' own trees until they meet yours.
    """
    idx = idx or kinship.load_tree(store)
    tasks: List[ResearchTask] = []

    if root and root in idx.individuals:
        sourced = {
            r["subject_xref"]
            for r in store.db.execute("SELECT DISTINCT subject_xref FROM evidence")
            if r["subject_xref"]
        }
        for xref, depth in frontier(idx, root, max_generations):
            if xref == root:
                continue
            side = _side_of_ancestor(idx, root, xref)
            for task in suggest_for_person(idx, xref, depth, side):
                if xref in sourced:
                    task.priority -= 0.08
                tasks.append(task)

    if kit_id is not None:
        tasks.extend(_cluster_tasks(store, kit_id, idx))

    tasks.sort(key=lambda t: -t.priority)
    tasks = tasks[:limit]

    if persist:
        store.db.execute("DELETE FROM research_task WHERE status='open'")
        store.db.executemany(
            "INSERT INTO research_task(subject, subject_xref, question, record_set, "
            "repository, cost_band, priority, created_at) VALUES(?,?,?,?,?,?,?,?)",
            [
                (t.subject, t.subject_xref, t.question, t.record_set,
                 t.repository, t.cost_band, t.priority, _now())
                for t in tasks
            ],
        )
        store.commit()
    return tasks


def _side_of_ancestor(
    idx: kinship.TreeIndex, root: str, xref: str
) -> Optional[str]:
    """Whether an ancestor sits on the maternal or paternal line."""
    parents = idx.parents.get(root, [])
    for parent, _rel in parents:
        sex = (idx.individuals.get(parent, {}) or {}).get("sex")
        depths = kinship.ancestors(idx, parent)
        if xref in depths:
            if sex == "F":
                return "maternal"
            if sex == "M":
                return "paternal"
    return None


def _cluster_tasks(
    store: Store, kit_id: int, idx: kinship.TreeIndex
) -> List[ResearchTask]:
    out: List[ResearchTask] = []
    clusters = list(
        store.db.execute("SELECT * FROM cluster WHERE kit_id=?", (kit_id,))
    )
    for c in clusters:
        members = list(
            store.db.execute(
                "SELECT m.* FROM cluster_member cm JOIN match m ON m.id=cm.match_id "
                "WHERE cm.cluster_id=? ORDER BY m.total_cm DESC",
                (c["id"],),
            )
        )
        if not members:
            continue
        linked = [m for m in members if m["tree_xref"] and m["tree_xref"] in idx.individuals]
        if linked:
            continue
        top = members[0]
        side = c["side"] or "unknown"
        out.append(ResearchTask(
            subject=c["name"],
            subject_xref=None,
            question=(
                f"who is the common ancestor of this {side} cluster? "
                f"{len(members)} matches, none connected to your tree"
            ),
            record_set=(
                f"build {top['name']}'s ancestry back four or five generations "
                "and look for a surname or parish shared with the rest of the "
                "cluster"
            ),
            repository="Ancestry / FindMyPast / FamilySearch trees",
            cost_band="subscription",
            priority=0.88 if side == "paternal" else 0.72,
            rationale=(
                "a cluster with no tree links is a branch you have no "
                "documentation for; its members collectively do"
            ),
        ))
    return out


# ---------------------------------------------------------------------------
# importing citations from a GEDCOM


def import_gedcom_sources(store: Store, path: str) -> Dict[str, int]:
    """Pull SOUR records and their citations out of a GEDCOM.

    Exports from the subscription sites carry citations for every fact they
    attached, and those citations are the record of what has already been
    looked up -- which is exactly what stops you paying twice for the same
    certificate. The tree importer deliberately keeps only scalar values, so
    the citation structure is read separately here.
    """
    from .tree.gedcom import _parse_nodes, _read_text, _bare

    warnings: List[str] = []
    text = _read_text(path, warnings)
    nodes = _parse_nodes(text, warnings)
    stats = {"sources": 0, "citations": 0}

    source_ids: Dict[str, int] = {}
    for node in nodes:
        if node.tag != "SOUR" or not node.xref:
            continue
        title = node.value_of("TITL") or node.value_of("ABBR") or node.xref
        repo = node.value_of("REPO") or node.value_of("AUTH")
        sid = add_source(
            store, title,
            repository=repo,
            reference=node.value_of("PUBL"),
            notes=node.value_of("NOTE"),
            record_type=_guess_record_type(title),
        )
        source_ids[node.xref] = sid
        stats["sources"] += 1

    for node in nodes:
        if node.tag != "INDI" or not node.xref:
            continue
        subject = _bare(node.xref)
        for citation, claim_type in _walk_citations(node):
            sid = source_ids.get(citation)
            if sid is None:
                continue
            add_evidence(
                store, sid, subject, claim_type,
                f"cited for {claim_type}", supports=True, confidence="direct",
            )
            stats["citations"] += 1
    return stats


def _walk_citations(node, claim_type: str = "identity"):
    """Yield ``(source_xref, claim_type)`` for every SOUR pointer beneath a node."""
    tag_to_claim = {
        "BIRT": "birth", "CHR": "birth", "DEAT": "death", "BURI": "death",
        "MARR": "marriage", "RESI": "residence", "CENS": "residence",
        "OCCU": "occupation", "NAME": "identity",
    }
    for child in node.children:
        if child.tag == "SOUR" and child.value.startswith("@"):
            yield child.value.strip(), claim_type
        else:
            yield from _walk_citations(child, tag_to_claim.get(child.tag, claim_type))


def _guess_record_type(title: str) -> Optional[str]:
    t = (title or "").lower()
    for needle, kind in (
        ("census", "census"), ("birth", "birth"), ("baptis", "birth"),
        ("marriage", "marriage"), ("death", "death"), ("burial", "death"),
        ("probate", "will"), ("will", "will"), ("parish", "parish"),
        ("newspaper", "newspaper"), ("military", "military"),
    ):
        if needle in t:
            return kind
    return None

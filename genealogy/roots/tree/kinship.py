"""Navigate the family tree: ancestors, descendants, and how two people relate.

The tree side of the project answers a different question from the DNA side.
DNA tells you *how much* two people share and therefore how closely they are
related; the tree tells you *through whom*.  The interesting work happens
where the two disagree.

Two design decisions are worth stating.

Adoptive, foster, and step relationships are tracked but excluded from
genetic reckoning by default.  A GEDCOM records them with a PEDI tag under
the child link, and treating an adoptive parent as a genetic one would
silently invalidate every prediction downstream.  Pass ``genetic_only=False``
when you want the documented family as recorded rather than the genetic one.

Relationship finding returns *all* minimal paths, not one.  In any tree that
goes back far enough -- and in any endogamous community, immediately -- two
people are related through several lines at once, and their DNA reflects the
sum of those lines.  Reporting only the closest connection systematically
underestimates expected sharing, which is the usual reason a documented
relationship looks like it "shares too much DNA".
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from ..dna.simulate import Relationship, relationship_for
from ..store import Store

NON_GENETIC_LINKS = {"adopted", "foster", "step", "sealing"}


@dataclass
class TreeIndex:
    """In-memory view of the tree, built once and queried many times."""

    individuals: Dict[str, Dict] = field(default_factory=dict)
    families: Dict[str, Dict] = field(default_factory=dict)
    parents: Dict[str, List[Tuple[str, str]]] = field(default_factory=dict)
    children: Dict[str, List[Tuple[str, str]]] = field(default_factory=dict)
    family_of_child: Dict[str, List[str]] = field(default_factory=dict)
    couple_of: Dict[str, Set[str]] = field(default_factory=dict)

    def name(self, xref: str) -> str:
        ind = self.individuals.get(xref)
        if not ind:
            return xref
        given = (ind.get("given") or "").strip()
        surname = (ind.get("surname") or "").strip()
        full = f"{given} {surname}".strip()
        years = ""
        b, d = ind.get("birth_year"), ind.get("death_year")
        if b or d:
            years = f" ({b or '?'}-{d or '?'})"
        return (full or xref) + years

    def size(self) -> Tuple[int, int]:
        return len(self.individuals), len(self.families)


def load_tree(store: Store) -> TreeIndex:
    idx = TreeIndex()
    for row in store.individuals():
        idx.individuals[row["xref"]] = dict(row)
    for row in store.families():
        fam = dict(row)
        idx.families[row["xref"]] = fam
        husb, wife = fam.get("husb"), fam.get("wife")
        if husb and wife:
            idx.couple_of.setdefault(husb, set()).add(wife)
            idx.couple_of.setdefault(wife, set()).add(husb)
    for row in store.children():
        fam_xref, indi, rel = row["fam_xref"], row["indi_xref"], (row["rel"] or "birth")
        fam = idx.families.get(fam_xref)
        idx.family_of_child.setdefault(indi, []).append(fam_xref)
        idx.children.setdefault(fam_xref, []).append((indi, rel))
        if not fam:
            continue
        for parent in (fam.get("husb"), fam.get("wife")):
            if parent:
                idx.parents.setdefault(indi, []).append((parent, rel))
    return idx


def ancestors(
    idx: TreeIndex,
    xref: str,
    max_generations: int = 20,
    genetic_only: bool = True,
) -> Dict[str, int]:
    """Map every ancestor to the fewest generations separating them.

    Includes the person themselves at depth 0, which makes direct-line
    relationships fall out of the same code path as collateral ones.
    """
    depths: Dict[str, int] = {xref: 0}
    queue = deque([(xref, 0)])
    while queue:
        current, depth = queue.popleft()
        if depth >= max_generations:
            continue
        for parent, rel in idx.parents.get(current, []):
            if genetic_only and rel in NON_GENETIC_LINKS:
                continue
            known = depths.get(parent)
            if known is None or depth + 1 < known:
                depths[parent] = depth + 1
                queue.append((parent, depth + 1))
    return depths


def descendants(
    idx: TreeIndex, xref: str, max_generations: int = 20, genetic_only: bool = True
) -> Dict[str, int]:
    depths: Dict[str, int] = {xref: 0}
    queue = deque([(xref, 0)])
    while queue:
        current, depth = queue.popleft()
        if depth >= max_generations:
            continue
        for fam in _families_as_spouse(idx, current):
            for child, rel in idx.children.get(fam, []):
                if genetic_only and rel in NON_GENETIC_LINKS:
                    continue
                known = depths.get(child)
                if known is None or depth + 1 < known:
                    depths[child] = depth + 1
                    queue.append((child, depth + 1))
    return depths


def _families_as_spouse(idx: TreeIndex, xref: str) -> List[str]:
    return [
        fam_xref
        for fam_xref, fam in idx.families.items()
        if fam.get("husb") == xref or fam.get("wife") == xref
    ]


@dataclass
class KinPath:
    """One way in which two people are related."""

    mrcas: List[str]
    up: int
    down: int
    ancestors: int
    relationship: Relationship

    def describe(self, idx: TreeIndex) -> str:
        who = " and ".join(idx.name(x) for x in self.mrcas)
        return f"{self.relationship.name} through {who}"


def relationship_between(
    idx: TreeIndex,
    a: str,
    b: str,
    max_generations: int = 20,
    genetic_only: bool = True,
    max_paths: int = 6,
) -> List[KinPath]:
    """Every minimal-distance path connecting two people in the tree.

    Returns an empty list when the tree records no connection, which is the
    normal case for a new match and is not evidence of anything.
    """
    if a == b:
        return []
    anc_a = ancestors(idx, a, max_generations, genetic_only)
    anc_b = ancestors(idx, b, max_generations, genetic_only)
    common = set(anc_a) & set(anc_b)
    if not common:
        return []

    # Discard ancestors that are only reachable *through* another common
    # ancestor: a great-grandparent is not an independent connection when
    # their child is already a shared ancestor.
    minimal: List[str] = []
    for x in common:
        redundant = False
        for y in common:
            if x == y:
                continue
            if anc_a[y] < anc_a[x] and anc_b[y] < anc_b[x]:
                anc_of_y = ancestors(idx, y, max_generations, genetic_only)
                if x in anc_of_y:
                    redundant = True
                    break
        if not redundant:
            minimal.append(x)

    grouped: Dict[Tuple[int, int], List[str]] = {}
    for x in minimal:
        grouped.setdefault((anc_a[x], anc_b[x]), []).append(x)

    paths: List[KinPath] = []
    for (up, down), people in sorted(grouped.items(), key=lambda kv: sum(kv[0])):
        for couple in _split_into_couples(idx, people):
            n = len(couple)
            paths.append(
                KinPath(
                    mrcas=couple,
                    up=up,
                    down=down,
                    ancestors=min(n, 2),
                    relationship=relationship_for(up, down, min(n, 2)),
                )
            )
    paths.sort(key=lambda p: (p.up + p.down, -p.ancestors))
    return paths[:max_paths]


def _split_into_couples(idx: TreeIndex, people: Sequence[str]) -> List[List[str]]:
    """Group equally-distant common ancestors into couples where they married.

    Two shared grandparents who were married to each other are one
    connection (a full relationship); two who were not are two separate
    half relationships.
    """
    remaining = list(people)
    out: List[List[str]] = []
    while remaining:
        person = remaining.pop(0)
        partner = None
        for candidate in remaining:
            if candidate in idx.couple_of.get(person, set()):
                partner = candidate
                break
        if partner:
            remaining.remove(partner)
            out.append([person, partner])
        else:
            out.append([person])
    return out


def combined_relationship(paths: Sequence[KinPath]) -> Optional[Relationship]:
    """The single closest path, for labelling purposes."""
    if not paths:
        return None
    return paths[0].relationship


def find_people(idx: TreeIndex, query: str, limit: int = 20) -> List[Tuple[str, str]]:
    """Loose name search across the tree."""
    q = query.strip().lower()
    hits: List[Tuple[str, str]] = []
    for xref, ind in idx.individuals.items():
        name = f"{ind.get('given') or ''} {ind.get('surname') or ''}".strip().lower()
        if q in name or q in xref.lower():
            hits.append((xref, idx.name(xref)))
        if len(hits) >= limit:
            break
    return hits


def ancestral_couples(
    idx: TreeIndex, root: str, generation: int, genetic_only: bool = True
) -> List[List[str]]:
    """The ancestral couples sitting exactly ``generation`` steps above someone.

    Used to enumerate candidate common ancestors when the DNA says "your
    shared ancestor is about four generations up" and you want the list of
    couples that could be it.
    """
    depths = ancestors(idx, root, generation, genetic_only)
    at_level = [x for x, d in depths.items() if d == generation]
    return _split_into_couples(idx, at_level)


def line_of_descent(
    idx: TreeIndex, ancestor: str, descendant: str, genetic_only: bool = True
) -> List[str]:
    """The chain of people from an ancestor down to a descendant."""
    queue = deque([[ancestor]])
    seen: Set[str] = {ancestor}
    while queue:
        path = queue.popleft()
        current = path[-1]
        if current == descendant:
            return path
        for fam in _families_as_spouse(idx, current):
            for child, rel in idx.children.get(fam, []):
                if genetic_only and rel in NON_GENETIC_LINKS:
                    continue
                if child in seen:
                    continue
                seen.add(child)
                queue.append(path + [child])
    return []


def tree_gaps(idx: TreeIndex, root: str, max_generations: int = 8) -> Dict[int, Tuple[int, int]]:
    """How complete the tree is at each generation above a person.

    Returns ``{generation: (known, possible)}``.  Completeness falls off a
    cliff somewhere, and where it falls off is exactly where DNA evidence
    has to take over from documents.
    """
    depths = ancestors(idx, root, max_generations)
    out: Dict[int, Tuple[int, int]] = {}
    for gen in range(1, max_generations + 1):
        known = sum(1 for d in depths.values() if d == gen)
        out[gen] = (known, 2 ** gen)
    return out

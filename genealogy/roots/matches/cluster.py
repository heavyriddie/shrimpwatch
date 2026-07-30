"""Cluster matches into ancestral lines from shared-match data.

This is the Leeds method, automated.  The observation behind it: two of your
matches who descend from the same ancestral couple will usually match each
other as well as matching you, while two matches from unrelated branches
will not.  Group your matches by who-matches-whom and the groups fall out
along ancestral lines -- typically four large ones corresponding to your four
grandparents, plus smaller ones for more distant branches.

Its great virtue is that it needs no segment data, which makes it the main
analytical tool for Ancestry results, where segment data is not available at
all.  All it needs is the shared-match list for each of your matches.

Choice of cM band matters.  Too high and everyone is in everyone's shared
list because they are all close relatives; too low and the shared-match
lists are dominated by false or unresolvably distant connections.  The
classic Leeds band is 90-400 cM; automated clustering tools usually open it
up somewhat.  We default to 40-400 cM and let you move it.

Clustering method: average-linkage agglomerative clustering on Jaccard
distance between shared-match sets, computed with the nearest-neighbour
chain algorithm so it stays O(n^2) rather than O(n^3) and can handle a few
thousand matches.  Average linkage suits this problem because clusters
correspond to descent groups of varying size and density; single linkage
would chain unrelated branches together through a single doubly-related
person, which is exactly the failure mode you hit in any endogamous line.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from ..store import Store


@dataclass
class Cluster:
    index: int
    match_ids: List[int] = field(default_factory=list)
    names: List[str] = field(default_factory=list)
    side: Optional[str] = None
    side_support: Tuple[int, int] = (0, 0)
    total_cm: float = 0.0
    top_cm: float = 0.0
    cohesion: float = 0.0

    @property
    def size(self) -> int:
        return len(self.match_ids)

    def describe(self) -> str:
        side = self.side or "unassigned"
        return (
            f"cluster {self.index}: {self.size} matches, {side}, "
            f"strongest {self.top_cm:.0f} cM, cohesion {self.cohesion:.0%}"
        )


@dataclass
class ClusterResult:
    clusters: List[Cluster] = field(default_factory=list)
    unclustered: List[int] = field(default_factory=list)
    order: List[int] = field(default_factory=list)
    edges: Set[Tuple[int, int]] = field(default_factory=set)
    labels: Dict[int, str] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    band: Tuple[float, float] = (0.0, 0.0)

    def summary(self) -> str:
        sided = sum(1 for c in self.clusters if c.side)
        return (
            f"{len(self.clusters)} clusters over "
            f"{sum(c.size for c in self.clusters)} matches "
            f"({len(self.unclustered)} unclustered, {sided} side-assigned)"
        )


def cluster_matches(
    store: Store,
    kit_id: int,
    min_cm: float = 40.0,
    max_cm: float = 400.0,
    distance_cutoff: float = 0.80,
    min_cluster_size: int = 2,
    max_matches: int = 2000,
    source: Optional[str] = None,
) -> ClusterResult:
    """Group a kit's matches into probable ancestral lines."""
    result = ClusterResult(band=(min_cm, max_cm))
    rows = [
        m for m in store.matches(kit_id, source=source)
        if m["total_cm"] is not None and min_cm <= m["total_cm"] <= max_cm
    ]
    if len(rows) > max_matches:
        rows = sorted(rows, key=lambda m: -(m["total_cm"] or 0))[:max_matches]
        result.warnings.append(
            f"limited to the {max_matches} strongest matches in the band; "
            "raise --max-matches if you need the full list"
        )
    if len(rows) < 3:
        result.warnings.append(
            f"only {len(rows)} matches fall between {min_cm:.0f} and "
            f"{max_cm:.0f} cM; widen the band"
        )
        return result

    ids = [m["id"] for m in rows]
    index_of = {mid: i for i, mid in enumerate(ids)}
    result.labels = {m["id"]: (m["name"] or str(m["remote_id"])) for m in rows}
    sides = {m["id"]: m["side"] for m in rows}
    cms = {m["id"]: (m["total_cm"] or 0.0) for m in rows}

    neighbours: List[Set[int]] = [set() for _ in ids]
    for edge in store.shared_matches(kit_id):
        a, b = edge["a_id"], edge["b_id"]
        ia, ib = index_of.get(a), index_of.get(b)
        if ia is None or ib is None or ia == ib:
            continue
        neighbours[ia].add(ib)
        neighbours[ib].add(ia)
        result.edges.add((min(ia, ib), max(ia, ib)))

    if not result.edges:
        result.warnings.append(
            "no shared-match data for these matches -- clustering needs an "
            "in-common-with export, not just the match list"
        )
        return result

    # A match is its own shared match; including self makes two matches that
    # match each other but nobody else come out as similar rather than
    # maximally distant.
    for i, nb in enumerate(neighbours):
        nb.add(i)

    n = len(ids)
    dist = _distance_matrix(neighbours)
    merges = _nn_chain(n, dist)
    groups = _cut(n, merges, distance_cutoff)

    ordered_groups = sorted(
        (g for g in groups if len(g) >= min_cluster_size),
        key=lambda g: (-len(g), -max(cms[ids[i]] for i in g)),
    )
    for pos, g in enumerate(ordered_groups, start=1):
        members = sorted(g, key=lambda i: -cms[ids[i]])
        c = Cluster(index=pos)
        for i in members:
            c.match_ids.append(ids[i])
            c.names.append(result.labels[ids[i]])
        c.total_cm = sum(cms[mid] for mid in c.match_ids)
        c.top_cm = max(cms[mid] for mid in c.match_ids)
        c.cohesion = _cohesion(members, neighbours)
        c.side, c.side_support = _majority_side([sides[mid] for mid in c.match_ids])
        result.clusters.append(c)
        result.order.extend(c.match_ids)

    clustered = {mid for c in result.clusters for mid in c.match_ids}
    result.unclustered = [mid for mid in ids if mid not in clustered]
    result.order.extend(result.unclustered)
    return result


def _distance_matrix(neighbours: Sequence[Set[int]]) -> List[List[float]]:
    """Jaccard distance between every pair of shared-match sets."""
    n = len(neighbours)
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        ni = neighbours[i]
        for j in range(i + 1, n):
            nj = neighbours[j]
            inter = len(ni & nj)
            union = len(ni) + len(nj) - inter
            d = 1.0 - (inter / union if union else 0.0)
            dist[i][j] = dist[j][i] = d
    return dist


def _nn_chain(n: int, dist: List[List[float]]) -> List[Tuple[int, int, float]]:
    """Average-linkage agglomeration via the nearest-neighbour chain algorithm.

    Average linkage is "reducible", which is what makes the chain algorithm
    valid: merging two mutual nearest neighbours never brings some third
    cluster closer than either was before, so the merges the chain finds are
    exactly the merges a naive closest-pair search would find, at O(n^2)
    instead of O(n^3).
    """
    active = set(range(n))
    size = [1] * n
    chain: List[int] = []
    merges: List[Tuple[int, int, float]] = []

    while len(active) > 1:
        if not chain:
            chain = [next(iter(active))]
        a = chain[-1]
        best, b = math.inf, -1
        row = dist[a]
        for c in active:
            if c == a:
                continue
            d = row[c]
            if d < best or (d == best and (b < 0 or c < b)):
                best, b = d, c
        if len(chain) >= 2 and b == chain[-2]:
            chain.pop()
            chain.pop()
            merges.append((a, b, best))
            na, nb = size[a], size[b]
            row_b = dist[b]
            for c in active:
                if c == a or c == b:
                    continue
                merged = (na * row[c] + nb * row_b[c]) / (na + nb)
                row[c] = dist[c][a] = merged
            size[a] = na + nb
            active.discard(b)
        else:
            chain.append(b)
    return merges


def _cut(n: int, merges: Sequence[Tuple[int, int, float]], cutoff: float) -> List[Set[int]]:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, d in sorted(merges, key=lambda m: m[2]):
        if d > cutoff:
            break
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    groups: Dict[int, Set[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), set()).add(i)
    return list(groups.values())


def _cohesion(members: Sequence[int], neighbours: Sequence[Set[int]]) -> float:
    """Share of possible within-cluster shared-match links that are present.

    A tight cluster descends from one ancestral couple and approaches 1.0.
    A loose one has probably swallowed two branches and should be split by
    lowering the distance cutoff.
    """
    k = len(members)
    if k < 2:
        return 1.0
    present = 0
    member_set = set(members)
    for i in members:
        present += len((neighbours[i] & member_set) - {i})
    possible = k * (k - 1)
    return present / possible if possible else 0.0


def _majority_side(sides: Iterable[Optional[str]]) -> Tuple[Optional[str], Tuple[int, int]]:
    maternal = sum(1 for s in sides if s == "maternal")
    paternal = sum(1 for s in sides if s == "paternal")
    if maternal == 0 and paternal == 0:
        return None, (0, 0)
    if maternal > paternal * 2:
        return "maternal", (maternal, paternal)
    if paternal > maternal * 2:
        return "paternal", (maternal, paternal)
    return "mixed", (maternal, paternal)


def persist(store: Store, kit_id: int, result: ClusterResult, run: str = "default") -> None:
    """Save clusters so later commands and reports can refer to them."""
    store.db.execute(
        "DELETE FROM cluster_member WHERE cluster_id IN "
        "(SELECT id FROM cluster WHERE kit_id=? AND run=?)",
        (kit_id, run),
    )
    store.db.execute("DELETE FROM cluster WHERE kit_id=? AND run=?", (kit_id, run))
    for c in result.clusters:
        cur = store.db.execute(
            "INSERT INTO cluster(kit_id, run, name, side, ancestral_hint) VALUES(?,?,?,?,?)",
            (kit_id, run, f"cluster {c.index}", c.side, None),
        )
        cid = int(cur.lastrowid)
        store.db.executemany(
            "INSERT OR REPLACE INTO cluster_member(cluster_id, match_id) VALUES(?,?)",
            [(cid, mid) for mid in c.match_ids],
        )
    store.commit()


def grandparent_quadrants(result: ClusterResult) -> List[Cluster]:
    """The four largest clusters, which usually map to the four grandparents.

    Worth stating plainly: this is a heuristic about typical results, not a
    property of the method.  Branches that tested heavily produce big
    clusters and branches that did not produce none, so a person can easily
    see six clusters from one grandparent and zero from another.
    """
    return sorted(result.clusters, key=lambda c: -c.size)[:4]

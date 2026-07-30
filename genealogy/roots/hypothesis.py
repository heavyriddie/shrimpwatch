"""Turn data into questions worth investigating.

The other modules produce facts: this person shares 212 cM, that cluster has
eleven members, this region of chromosome 4 is unattributed.  This module
does the part a researcher actually wants -- reading those facts against each
other and surfacing the places where they disagree, or where a small amount
of extra work would resolve a lot of uncertainty.

The checks it runs, and what each is really looking for:

*Documented relationships against measured sharing.*  A paper-trail second
cousin who shares 900 cM is not a second cousin.  Somewhere in the four
links between you there is either an undocumented additional relationship or
a wrong parent.  Conversely a documented first cousin sharing nothing at all
is a near-certain misattributed parentage, because first cousins always
share DNA.  Distance matters here: at third cousin and beyond, sharing
nothing is unremarkable and means nothing.

*Projection onto an untested parent.*  Every paternal-side match of yours is
also a match of your father's, one generation closer.  If a match is your
paternal first cousin once removed, they are your father's first cousin.
Translating your whole paternal match list into your father's match list is
the single most useful move when the father is the unknown, because it puts
the closest relatives into relationships that name actual people:
"your father's aunt" is a searchable claim in a way that "your first cousin
once removed" is not.

*Clusters against the tree.*  A cluster whose members trace to a known
ancestral couple identifies the whole cluster.  A large cluster with no tree
links at all is the highest-value target in the whole dataset: it is a
branch of your family you have no documentation for, and its members
collectively do.

*Unattributed DNA.*  The stretches of your genome that no match explains are
a direct map of where your research has not reached.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .dna.predict import fit_hypothesis, predict
from .dna.simulate import CATALOGUE, Relationship, relationship_for
from .genome import GeneticMap
from .matches import cluster as cluster_mod
from .matches import triangulate
from .store import Store
from .tree import kinship


@dataclass
class Finding:
    kind: str
    subject: str
    summary: str
    score: float
    detail: Dict[str, Any] = field(default_factory=dict)
    actions: List[str] = field(default_factory=list)

    def line(self) -> str:
        return f"[{self.score:.2f}] {self.subject}: {self.summary}"


def investigate(
    store: Store,
    kit_id: int,
    gmap: GeneticMap,
    *,
    iterations: int = 1500,
    min_cm: float = 7.0,
    persist: bool = True,
) -> List[Finding]:
    """Run every check and return findings ordered by how much they matter."""
    findings: List[Finding] = []
    idx = kinship.load_tree(store)
    kit = store.kit(kit_id)
    root = _tree_xref_for_kit(store, kit_id)

    findings.extend(
        check_documented_relationships(store, kit_id, gmap, idx, root, iterations, min_cm)
    )
    findings.extend(check_claimed_relationships(store, kit_id, gmap, iterations, min_cm))
    findings.extend(project_onto_untested_parent(store, kit_id, gmap, iterations, min_cm))
    findings.extend(label_clusters(store, kit_id, idx))
    findings.extend(unattributed_regions(store, kit_id, gmap))
    findings.extend(tree_completeness(store, idx, root))
    findings.extend(conflicting_overlap_groups(store, kit_id, gmap))

    findings.sort(key=lambda f: -f.score)
    if persist:
        store.clear_hypotheses()
        for f in findings:
            store.add_hypothesis(f.kind, f.subject, f.summary, f.score, f.detail)
    return findings


def _tree_xref_for_kit(store: Store, kit_id: int) -> Optional[str]:
    kit = store.kit(kit_id)
    if not kit or not kit.person_id:
        return None
    row = store.db.execute(
        "SELECT tree_xref FROM person WHERE id=?", (kit.person_id,)
    ).fetchone()
    return row["tree_xref"] if row else None


# ---------------------------------------------------------------------------
# documented relationships vs measured sharing


def check_documented_relationships(
    store: Store,
    kit_id: int,
    gmap: GeneticMap,
    idx: kinship.TreeIndex,
    root: Optional[str],
    iterations: int = 1500,
    min_cm: float = 7.0,
) -> List[Finding]:
    out: List[Finding] = []
    if not root:
        return out
    for m in store.matches(kit_id):
        xref = m["tree_xref"]
        if not xref or xref not in idx.individuals:
            continue
        paths = kinship.relationship_between(idx, root, xref)
        if not paths:
            continue
        observed = m["total_cm"] or 0.0
        name = m["name"] or xref

        if len(paths) > 1:
            out.append(
                Finding(
                    kind="multiple-lines",
                    subject=name,
                    summary=(
                        f"related through {len(paths)} separate lines in the tree "
                        f"({'; '.join(p.relationship.name for p in paths)}), so "
                        "expected sharing is the sum of all of them, not just the "
                        "closest"
                    ),
                    score=0.45,
                    detail={
                        "paths": [
                            {"relationship": p.relationship.name,
                             "mrcas": [idx.name(x) for x in p.mrcas]}
                            for p in paths
                        ],
                        "observed_cm": observed,
                    },
                )
            )

        best = paths[0]
        fit = fit_hypothesis(
            store, best.relationship, observed, gmap=gmap, min_cm=min_cm,
            iterations=iterations,
        )
        if fit.plausible:
            continue

        degree = best.relationship.meioses
        if observed > fit.expected_cm:
            score = 0.9 if degree >= 4 else 0.75
            out.append(
                Finding(
                    kind="excess-sharing",
                    subject=name,
                    summary=(
                        f"documented as {best.relationship.name} but shares "
                        f"{observed:.0f} cM, against a simulated average of "
                        f"{fit.expected_cm:.0f} cM. {fit.comment}"
                    ),
                    score=score,
                    detail={
                        "documented": best.relationship.name,
                        "observed_cm": observed,
                        "expected_cm": fit.expected_cm,
                        "percentile": fit.percentile,
                        "mrcas": [idx.name(x) for x in best.mrcas],
                    },
                    actions=[
                        "look for a second connection between the lines",
                        "verify each parent-child link on the documented path",
                    ],
                )
            )
        elif degree <= 4 and observed < fit.expected_cm:
            # Sharing far too little only means something for close kin;
            # third cousins and beyond genuinely miss each other.
            out.append(
                Finding(
                    kind="deficient-sharing",
                    subject=name,
                    summary=(
                        f"documented as {best.relationship.name} but shares only "
                        f"{observed:.0f} cM, below almost every simulated pair at "
                        f"that relationship (average {fit.expected_cm:.0f} cM). "
                        "A wrong parent somewhere on the path between you is the "
                        "usual explanation"
                    ),
                    score=0.85,
                    detail={
                        "documented": best.relationship.name,
                        "observed_cm": observed,
                        "expected_cm": fit.expected_cm,
                        "percentile": fit.percentile,
                        "path": [idx.name(x) for x in best.mrcas],
                    },
                    actions=[
                        "test another descendant of each intervening couple to "
                        "find which link fails",
                    ],
                )
            )
    return out


def check_claimed_relationships(
    store: Store,
    kit_id: int,
    gmap: GeneticMap,
    iterations: int = 1500,
    min_cm: float = 7.0,
) -> List[Finding]:
    """Test relationships the user asserted but has not placed in the tree.

    A claimed relationship is evidence to be checked, not a premise.  People
    routinely inherit a family story about how they are related to someone,
    and a story that survives three generations of retelling has had ample
    opportunity to drift.  Checking it costs nothing and occasionally
    overturns everything built on top of it.
    """
    out: List[Finding] = []
    for m in store.matches(kit_id):
        if not m["extra"]:
            continue
        try:
            extra = json.loads(m["extra"])
        except (ValueError, TypeError):
            continue
        key = extra.get("claimed_relationship") if isinstance(extra, dict) else None
        rel = CATALOGUE.get(key) if key else None
        if not rel:
            continue
        observed = m["total_cm"] or 0.0
        fit = fit_hypothesis(
            store, rel, observed, gmap=gmap, min_cm=min_cm, iterations=iterations
        )
        if fit.plausible:
            continue
        out.append(
            Finding(
                kind="claim-conflict",
                subject=m["name"] or str(m["remote_id"]),
                summary=(
                    f"recorded as {rel.name}, but {observed:.0f} cM does not fit: "
                    f"{fit.comment}"
                ),
                score=0.8,
                detail={
                    "claimed": rel.name,
                    "observed_cm": observed,
                    "expected_cm": fit.expected_cm,
                    "percentile": fit.percentile,
                },
                actions=["re-examine how this relationship was established"],
            )
        )
    return out


# ---------------------------------------------------------------------------
# projecting matches onto an untested parent


def project_to_parent(rel: Relationship) -> Optional[Relationship]:
    """Restate a relationship from the point of view of one's parent.

    A relationship reaches the parent by dropping one meiosis from your side
    of the path.  Your paternal first cousin once removed (two steps up to
    your shared great-grandparents, three steps down) is your father's first
    cousin (one step up on his side).
    """
    if rel.special or rel.up < 1:
        return None
    return relationship_for(rel.up - 1, rel.down, rel.ancestors)


def project_onto_untested_parent(
    store: Store,
    kit_id: int,
    gmap: GeneticMap,
    iterations: int = 1500,
    min_cm: float = 7.0,
    side: str = "paternal",
    top_n: int = 12,
) -> List[Finding]:
    """Restate the strongest matches on one side as that parent's relatives.

    This is the workhorse for an unknown parent.  Each finding names what the
    match most likely is *to the untested parent*, which is the form you can
    act on: a great-aunt has a documented family, a first cousin has a
    surname.
    """
    out: List[Finding] = []
    rows = [
        m for m in store.matches(kit_id)
        if m["side"] == side and (m["total_cm"] or 0) > 0
    ]
    if not rows:
        return out
    rows.sort(key=lambda m: -(m["total_cm"] or 0))

    for m in rows[:top_n]:
        observed = m["total_cm"] or 0.0
        pred = predict(
            store, observed, gmap=gmap, min_cm=min_cm, iterations=iterations,
            segments=m["seg_count"],
        )
        projected: List[Tuple[str, float]] = []
        for cand in pred.credible_set(0.85):
            parent_rel = project_to_parent(cand.rel)
            if parent_rel:
                projected.append((parent_rel.name, cand.posterior))
        if not projected:
            continue
        merged: Dict[str, float] = {}
        for name, p in projected:
            merged[name] = merged.get(name, 0.0) + p
        ranked = sorted(merged.items(), key=lambda kv: -kv[1])
        best_name, best_p = ranked[0]
        parent_word = "father" if side == "paternal" else "mother"

        out.append(
            Finding(
                kind="parent-projection",
                subject=m["name"] or str(m["remote_id"]),
                summary=(
                    f"shares {observed:.0f} cM with you on the {side} side, which "
                    f"makes them most likely your {parent_word}'s {best_name} "
                    f"({best_p:.0%} of the probability mass)"
                ),
                score=min(0.95, 0.4 + observed / 2000.0),
                detail={
                    "observed_cm": observed,
                    "side": side,
                    "candidates_for_parent": ranked[:5],
                    "candidates_for_self": [
                        (c.rel.name, round(c.posterior, 4)) for c in pred.top(5)
                    ],
                },
                actions=[
                    f"build this person's tree back far enough to reach your "
                    f"{parent_word}'s generation",
                    "check who else clusters with them",
                ],
            )
        )
    return out


# ---------------------------------------------------------------------------
# clusters against the tree


def label_clusters(
    store: Store, kit_id: int, idx: kinship.TreeIndex, run: str = "default"
) -> List[Finding]:
    """Attach an ancestral couple to each saved cluster, where possible."""
    out: List[Finding] = []
    clusters = list(
        store.db.execute("SELECT * FROM cluster WHERE kit_id=? AND run=?", (kit_id, run))
    )
    for c in clusters:
        members = list(
            store.db.execute(
                "SELECT m.* FROM cluster_member cm JOIN match m ON m.id=cm.match_id "
                "WHERE cm.cluster_id=?",
                (c["id"],),
            )
        )
        if not members:
            continue
        linked = [m for m in members if m["tree_xref"] and m["tree_xref"] in idx.individuals]
        total_cm = sum(m["total_cm"] or 0 for m in members)
        top = max(members, key=lambda m: m["total_cm"] or 0)

        if not linked:
            out.append(
                Finding(
                    kind="unidentified-cluster",
                    subject=c["name"],
                    summary=(
                        f"{len(members)} matches totalling {total_cm:.0f} cM, "
                        f"strongest {top['total_cm']:.0f} cM "
                        f"({top['name']}), and not one of them connects to your "
                        f"tree. This is a branch of your family you have no "
                        f"documentation for"
                    ),
                    score=min(0.95, 0.5 + len(members) / 40.0),
                    detail={
                        "side": c["side"],
                        "members": [m["name"] for m in members][:25],
                        "total_cm": total_cm,
                    },
                    actions=[
                        f"contact {top['name']} first -- they share the most DNA",
                        "look for a surname or birthplace common to the cluster",
                    ],
                )
            )
            continue

        common = _common_ancestors_of(idx, [m["tree_xref"] for m in linked])
        if common:
            names = ", ".join(idx.name(x) for x in common[:2])
            store.db.execute(
                "UPDATE cluster SET ancestral_hint=? WHERE id=?", (names, c["id"])
            )
            out.append(
                Finding(
                    kind="identified-cluster",
                    subject=c["name"],
                    summary=(
                        f"{len(members)} matches trace to {names}; the "
                        f"{len(members) - len(linked)} unlinked members of this "
                        "cluster almost certainly descend from the same couple"
                    ),
                    score=0.6,
                    detail={
                        "ancestors": [idx.name(x) for x in common],
                        "linked": len(linked),
                        "unlinked": len(members) - len(linked),
                    },
                    actions=["work the unlinked members' trees toward that couple"],
                )
            )
    store.commit()
    return out


def _common_ancestors_of(idx: kinship.TreeIndex, xrefs: Sequence[str]) -> List[str]:
    sets = []
    for x in xrefs:
        if x in idx.individuals:
            sets.append(set(kinship.ancestors(idx, x)))
    if not sets:
        return []
    shared = set.intersection(*sets) if len(sets) > 1 else sets[0]
    if not shared:
        return []
    depths = kinship.ancestors(idx, xrefs[0])
    return sorted(shared, key=lambda x: depths.get(x, 99))[:4]


# ---------------------------------------------------------------------------
# coverage and completeness


def unattributed_regions(store: Store, kit_id: int, gmap: GeneticMap) -> List[Finding]:
    painted = triangulate.paint(store, kit_id, gmap)
    if not painted:
        return []
    cov = triangulate.coverage(painted, gmap)
    genome = cov.get("genome", 0.0) or 1.0
    unattributed = cov.get("unattributed", 0.0)
    return [
        Finding(
            kind="coverage",
            subject="chromosome coverage",
            summary=(
                f"{cov['any']:.0f} cM of {genome:.0f} cM ({cov['any'] / genome:.0%}) "
                f"is covered by at least one match; {unattributed:.0f} cM "
                f"({unattributed / genome:.0%}) is explained by nobody. "
                f"Maternal side accounts for {cov['maternal']:.0f} cM, paternal "
                f"{cov['paternal']:.0f} cM"
            ),
            score=0.35,
            detail=cov,
            actions=["the unexplained portion is where the undiscovered branches are"],
        )
    ]


def tree_completeness(
    store: Store, idx: kinship.TreeIndex, root: Optional[str]
) -> List[Finding]:
    if not root or root not in idx.individuals:
        return []
    gaps = kinship.tree_gaps(idx, root)
    first_gap = None
    for gen in sorted(gaps):
        known, possible = gaps[gen]
        if known < possible:
            first_gap = (gen, known, possible)
            break
    if not first_gap:
        return []
    gen, known, possible = first_gap
    return [
        Finding(
            kind="tree-gap",
            subject="documented ancestry",
            summary=(
                f"generation {gen} above you is {known}/{possible} complete; this "
                "is the first level where the paper trail runs out, and matches "
                "whose common ancestor sits at or beyond it cannot be placed "
                "until it is filled"
            ),
            score=0.5,
            detail={"generations": {str(k): v for k, v in gaps.items()}},
        )
    ]


def conflicting_overlap_groups(
    store: Store, kit_id: int, gmap: GeneticMap, limit: int = 5
) -> List[Finding]:
    groups = triangulate.overlap_groups(store, kit_id, gmap)
    out: List[Finding] = []
    for g in groups[:limit]:
        if g.conflict:
            out.append(
                Finding(
                    kind="overlap-conflict",
                    subject=g.describe(),
                    summary=(
                        "this region has matches from both sides overlapping, so "
                        "the segments sit on opposite copies of the chromosome. "
                        "Treating the group as one ancestral line would merge two "
                        "unrelated families"
                    ),
                    score=0.55,
                    detail={"names": g.names, "sides": g.sides},
                )
            )
        elif g.size >= 3 and g.support in ("triangulated", "shared-match supported"):
            out.append(
                Finding(
                    kind="triangulation-group",
                    subject=g.describe(),
                    summary=(
                        f"{g.size} matches share this region and also match each "
                        f"other ({g.support}), so the segment descends from one "
                        "ancestral couple. Identifying the common ancestor of any "
                        "two of them identifies it for all"
                    ),
                    score=0.65,
                    detail={"names": g.names, "cm": g.cm, "side": g.side()},
                )
            )
    return out


# ---------------------------------------------------------------------------
# research planning


def research_plan(findings: Sequence[Finding], limit: int = 8) -> List[str]:
    """Condense findings into an ordered list of things worth doing next."""
    plan: List[str] = []
    seen: set = set()
    for f in findings:
        for action in f.actions:
            key = (f.kind, action)
            if key in seen:
                continue
            seen.add(key)
            plan.append(f"{action} ({f.subject})")
            if len(plan) >= limit:
                return plan
    return plan

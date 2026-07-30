"""Command line interface.

Every command operates on a single project file (a SQLite database, by
default ``roots.db`` in the working directory) which holds all genotypes,
matches, tree data and results.  Nothing in this toolkit makes a network
request: the project file is the whole system, and it stays on your machine.

Run ``roots`` with no arguments for the command list, or ``roots guide`` for
the recommended order of operations starting from two raw DNA downloads.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import textwrap
from typing import List, Optional, Sequence

from . import hypothesis as hypo
from .dna import phase as phase_mod
from .dna import predict as predict_mod
from .dna import qc as qc_mod
from .dna import raw as raw_mod
from .dna import relate as relate_mod
from .dna import simulate as sim_mod
from .genome import AUTOSOMES, GeneticMap
from .matches import cluster as cluster_mod
from .matches import ingest as ingest_mod
from .matches import sides as sides_mod
from .matches import triangulate as tri_mod
from .store import Store

DEFAULT_PROJECT = "roots.db"


# ---------------------------------------------------------------------------
# helpers


def _out(text: str = "") -> None:
    print(text)


def _rule(title: str) -> None:
    _out()
    _out(title)
    _out("-" * len(title))


def _open_store(args) -> Store:
    return Store(args.project)


def _load_map(store: Store, args) -> GeneticMap:
    path = getattr(args, "map", None) or store.get_meta("genetic_map")
    build = getattr(args, "build", None) or store.get_meta("build") or "37"
    if path and os.path.exists(path):
        return GeneticMap.load(path, build=build)
    if path:
        _out(f"warning: genetic map {path!r} not found; using the built-in approximation")
    return GeneticMap.linear(build=build)

def _warn_if_approximate(gmap: GeneticMap) -> None:
    if gmap.is_approximate:
        _out(
            "note: using a uniform-rate genetic map. Centimorgan figures computed "
            "here are approximate. Load a real recombination map with "
            "--map for publication-grade numbers (see docs/DATA_SOURCES.md)."
        )


def _resolve_kit(store: Store, ref: str):
    kit = store.kit(ref)
    if not kit:
        raise SystemExit(f"no kit named {ref!r}. Run 'roots kits' to list them.")
    return kit


def _table(rows: Sequence[Sequence[object]], headers: Sequence[str]) -> None:
    if not rows:
        _out("(none)")
        return
    cells = [[str(c) for c in row] for row in rows]
    widths = [len(h) for h in headers]
    for row in cells:
        for i, c in enumerate(row):
            if i < len(widths):
                widths[i] = max(widths[i], len(c))
    _out("  ".join(h.ljust(widths[i]) for i, h in enumerate(headers)))
    _out("  ".join("-" * w for w in widths))
    for row in cells:
        _out("  ".join(c.ljust(widths[i]) for i, c in enumerate(row)))


# ---------------------------------------------------------------------------
# commands


def cmd_init(args) -> int:
    store = _open_store(args)
    if args.map:
        store.set_meta("genetic_map", os.path.abspath(args.map))
    store.set_meta("build", args.build)
    _out(f"project ready at {os.path.abspath(args.project)}")
    _out("next: roots import-dna <file> --person self")
    store.close()
    return 0


def cmd_import_dna(args) -> int:
    store = _open_store(args)
    info = raw_mod.inspect(args.path)
    if args.vendor:
        info.vendor = args.vendor
    label = args.label or os.path.splitext(os.path.basename(args.path))[0]

    # Read a prefix to settle the assembly before committing anything.
    _info, records, _stats = raw_mod.read_with_stats(args.path, info)
    sample = []
    for i, rec in enumerate(records):
        sample.append((rec[0], rec[1]))
        if i > 20000:
            break
    build = args.build_override or raw_mod.resolve_build(info, sample)

    info2, records2, stats = raw_mod.read_with_stats(args.path, info)
    kit_id = store.create_kit(label, args.person, info.vendor, build, os.path.abspath(args.path))
    written = store.insert_genotypes(kit_id, records2)
    kit = store.kit(kit_id)
    q = qc_mod.qc_kit(store, kit)
    store.finalize_kit(kit_id, q.total, q.called, q.inferred_sex)

    _out(f"imported {info.summary()}")
    _out(f"  kit '{label}' for person '{args.person}': {written:,} positions, "
         f"call rate {q.call_rate:.2%}, inferred sex {q.inferred_sex or 'unknown'}")
    if stats.duplicates:
        _out(f"  {stats.duplicates:,} duplicate positions collapsed")
    for w in q.warnings():
        _out(f"  warning: {w}")
    store.close()
    return 0


def cmd_kits(args) -> int:
    store = _open_store(args)
    rows = []
    for k in store.kits():
        q = qc_mod.qc_kit(store, k)
        rows.append([
            k.id, k.label, k.person or "?", k.vendor or "?", f"b{k.build}",
            f"{q.total:,}", f"{q.call_rate:.1%}", f"{q.het_rate:.1%}",
            k.inferred_sex or "?",
        ])
    _table(rows, ["id", "label", "person", "vendor", "build", "snps", "call", "het", "sex"])
    store.close()
    return 0


def cmd_qc(args) -> int:
    store = _open_store(args)
    kits = [_resolve_kit(store, args.kit)] if args.kit else store.kits()
    for k in kits:
        q = qc_mod.qc_kit(store, k)
        _rule(f"{k.display}")
        _out(f"positions   {q.total:,} ({q.called:,} called, {q.call_rate:.2%})")
        _out(f"heterozygosity {q.het_rate:.2%}")
        _out(f"inferred sex   {q.inferred_sex or 'unknown'}"
             f"  (X het {q.x_het}/{q.x_called}, Y called {q.y_called}/{q.y_total})")
        if q.indel:
            _out(f"indel calls    {q.indel:,} (not comparable across vendors)")
        for w in q.warnings():
            _out(f"warning: {w}")
        if not q.warnings():
            _out("no problems found")
    store.close()
    return 0


def cmd_merge_kits(args) -> int:
    store = _open_store(args)
    kits = store.kits_for_person(args.person)
    if len(kits) < 2:
        raise SystemExit(
            f"person {args.person!r} has {len(kits)} kit(s); merging needs at least two"
        )
    _out(f"checking that the {len(kits)} kits really are the same person")
    base = kits[0]
    for other in kits[1:]:
        c = qc_mod.concordance(store, base.id, other.id)
        _out(f"  {base.label} vs {other.label}: {c.shared:,} shared SNPs, "
             f"{c.rate:.3%} identical -> {c.verdict()}")
        if c.verdict() == "different people":
            raise SystemExit("refusing to merge kits that are not the same person")
    label = args.label or f"{args.person}-merged"
    new_id, stats = qc_mod.merge_kits(
        store, [k.id for k in kits], label, args.person, build=base.build or "37"
    )
    _out(f"merged into '{label}': {stats['positions']:,} positions "
         f"({stats['conflict']:,} conflicting sites set to no-call)")
    gained = stats["positions"] - max(k.snp_count or 0 for k in kits)
    _out(f"that is {gained:,} more usable positions than the best single kit")
    carried = ingest_mod.carry_matches(store, [k.id for k in kits], new_id)
    if carried["matches"]:
        _out(f"carried {carried['matches']} matches, {carried['segments']} segments "
             f"and {carried['shared']} shared-match links onto the merged kit")
    store.close()
    return 0


def cmd_compare(args) -> int:
    store = _open_store(args)
    gmap = _load_map(store, args)
    _warn_if_approximate(gmap)
    a = _resolve_kit(store, args.kit_a)
    b = _resolve_kit(store, args.kit_b)
    _out(f"comparing {a.display} with {b.display} ...")
    cmp = relate_mod.compare(
        store, a.id, b.id, gmap, min_cm=args.min_cm, min_snps=args.min_snps
    )
    verdict = relate_mod.classify_close(cmp, gmap)
    _rule("result")
    _out(f"shared SNPs      {cmp.shared_snps:,}")
    _out(f"half-identical   {cmp.total_cm:,.0f} cM over {len(cmp.hir)} segments "
         f"(largest {cmp.largest_cm:.0f} cM)")
    _out(f"fully identical  {cmp.total_fir_cm:,.0f} cM over {len(cmp.fir)} segments")
    _out(f"opposite homozygotes {cmp.ibs0:,} ({cmp.ibs0_rate:.4%})")
    if cmp.masked_errors:
        _out(f"({cmp.masked_errors:,} isolated mismatches treated as genotyping error)")
    _rule("interpretation")
    _out(f"{verdict.label}  [confidence: {verdict.confidence}]")
    _out(textwrap.fill(verdict.reasoning, 78))
    if args.save:
        store.clear_ibd(a.id, b.id)
        store.insert_ibd([s.as_row(a.id, b.id) for s in cmp.segments])
        _out(f"saved {len(cmp.segments)} segments to the project")
    store.close()
    return 0


def cmd_phase(args) -> int:
    store = _open_store(args)
    gmap = _load_map(store, args)
    child = _resolve_kit(store, args.child)
    parent = _resolve_kit(store, args.parent)

    if not args.skip_verify:
        _out("verifying the parent-child relationship before phasing ...")
        cmp = relate_mod.compare(store, child.id, parent.id, gmap, min_cm=args.min_cm)
        ok, why = relate_mod.verify_parent(cmp, gmap)
        _out(f"  {why}")
        if not ok:
            raise SystemExit(
                "refusing to phase: everything downstream assumes a true "
                "parent-child pair. Re-run with --skip-verify to override."
            )

    _out(f"phasing {child.display} against {parent.display} ({args.role}) ...")
    stats = phase_mod.phase_against_parent(
        store, child.id, parent.id, parent_role=args.role, child_sex=child.inferred_sex
    )
    _rule("phasing result")
    _out(f"sites processed        {stats.sites:,}")
    _out(f"homozygous (trivial)   {stats.homozygous:,}")
    _out(f"heterozygous resolved  {stats.het_resolved:,} of {stats.het_total:,} "
         f"({stats.het_resolution_rate:.1%})")
    _out(f"ambiguous (both het)   {stats.het_ambiguous:,}")
    _out(f"Mendelian errors       {stats.mendelian_errors:,} ({stats.mendelian_rate:.3%})")
    if stats.x_maternal:
        _out(f"X hemizygous (son)     {stats.x_maternal:,} sites, wholly maternal")
    _out(f"fully phased           {stats.phased_fraction:.1%} of usable sites")
    for note in stats.notes():
        _out(f"warning: {textwrap.fill(note, 74)}")

    other = "paternal" if args.role == "mother" else "maternal"
    if args.export_other:
        counts = phase_mod.export_phased_kit(store, child.id, other, args.export_other)
        _out(f"wrote {other} pseudo-kit to {args.export_other} "
             f"({counts['resolved']:,} resolved of {counts['written']:,} rows)")
        _out(f"this file is half of your {'father' if other == 'paternal' else 'mother'}'s "
             "genome. Upload it to a site that accepts them and every match it "
             f"returns is a {other}-side relative.")
    if args.export_same:
        same = "maternal" if args.role == "mother" else "paternal"
        counts = phase_mod.export_phased_kit(store, child.id, same, args.export_same)
        _out(f"wrote {same} pseudo-kit to {args.export_same} "
             f"({counts['resolved']:,} resolved of {counts['written']:,} rows)")
    store.close()
    return 0


def cmd_import_matches(args) -> int:
    store = _open_store(args)
    kit = _resolve_kit(store, args.kit)
    report = ingest_mod.import_matches(store, kit.id, args.path, source=args.source)
    _out(f"{os.path.basename(args.path)} -> {kit.label} [{report.source}]: {report.summary()}")
    if report.columns:
        _out("  columns recognised: " + ", ".join(
            f"{k}={v!r}" for k, v in sorted(report.columns.items())
        ))
    for w in report.warnings:
        _out(f"  warning: {textwrap.fill(w, 74)}")
    store.close()
    return 0


def cmd_import_segments(args) -> int:
    store = _open_store(args)
    kit = _resolve_kit(store, args.kit)
    report = ingest_mod.import_segments(store, kit.id, args.path, source=args.source)
    _out(f"{os.path.basename(args.path)} -> {kit.label}: {report.summary()}")
    if report.unmatched_names:
        _out("  names with no match-list entry: "
             + ", ".join(report.unmatched_names[:8])
             + (" ..." if len(report.unmatched_names) > 8 else ""))
    for w in report.warnings:
        _out(f"  warning: {w}")
    store.close()
    return 0


def cmd_import_icw(args) -> int:
    store = _open_store(args)
    kit = _resolve_kit(store, args.kit)
    report = ingest_mod.import_shared_matches(store, kit.id, args.path, kind=args.kind)
    _out(f"{os.path.basename(args.path)} -> {kit.label}: {report.summary()}")
    for w in report.warnings:
        _out(f"  warning: {w}")
    store.close()
    return 0


def cmd_sides(args) -> int:
    store = _open_store(args)
    child = _resolve_kit(store, args.child)
    parent = _resolve_kit(store, args.parent)
    report = sides_mod.assign_from_parent_list(
        store, child.id, parent.id, parent_role=args.role,
        source=args.source, min_paternal_cm=args.floor,
        floor_multiple=args.floor_multiple,
    )
    _rule(f"side assignment for {child.display}")
    for w in report.warnings:
        _out(f"warning: {textwrap.fill(w, 76)}")
    _out(f"platforms compared: {report.source or 'none'}")
    _out(f"{args.role}'s list reaches down to {report.parent_floor_cm:.0f} cM")
    _out(f"result: {report.summary()}")

    flagged = report.flagged()
    if flagged:
        _rule("flagged")
        for c in flagged[:10]:
            _out(f"{c.name}: {textwrap.fill(c.reason, 72)}")

    if args.segments:
        extra = sides_mod.confirm_with_segments(store, child.id, parent.id, args.role)
        if extra:
            _rule("segment-level confirmation")
            for c in extra[:15]:
                _out(f"{c.name}: {c.side} -- {c.reason}")

    if args.verbose:
        _rule("all calls")
        rows = [[c.name[:32], f"{c.self_cm:.0f}",
                 f"{c.parent_cm:.0f}" if c.parent_cm is not None else "-",
                 c.side or "unknown", c.confidence] for c in report.calls[:60]]
        _table(rows, ["match", "cM(you)", f"cM({args.role})", "side", "confidence"])
    store.close()
    return 0


def cmd_cluster(args) -> int:
    store = _open_store(args)
    kit = _resolve_kit(store, args.kit)
    result = cluster_mod.cluster_matches(
        store, kit.id, min_cm=args.min_cm, max_cm=args.max_cm,
        distance_cutoff=args.cutoff, min_cluster_size=args.min_size,
        source=args.source,
    )
    for w in result.warnings:
        _out(f"warning: {textwrap.fill(w, 76)}")
    _rule(f"clusters for {kit.display} ({args.min_cm:.0f}-{args.max_cm:.0f} cM)")
    _out(result.summary())
    rows = []
    for c in result.clusters:
        rows.append([
            c.index, c.size, c.side or "-", f"{c.top_cm:.0f}",
            f"{c.total_cm:.0f}", f"{c.cohesion:.0%}",
            ", ".join(c.names[:3])[:46],
        ])
    _table(rows, ["#", "size", "side", "top cM", "total cM", "cohesion", "members"])
    if result.clusters:
        cluster_mod.persist(store, kit.id, result)
        _out()
        _out("saved. Clusters usually correspond to ancestral lines; the four "
             "largest often map to your four grandparents, though only if all "
             "four branches have tested relatives.")
    store.close()
    return 0


def cmd_triangulate(args) -> int:
    store = _open_store(args)
    gmap = _load_map(store, args)
    _warn_if_approximate(gmap)
    kit = _resolve_kit(store, args.kit)
    groups = tri_mod.overlap_groups(
        store, kit.id, gmap, min_members=args.min_members,
        min_overlap_cm=args.min_cm, source=args.source,
        max_match_cm=args.max_match_cm,
    )
    _rule(f"shared regions for {kit.display}")
    if not groups:
        _out("no overlapping segments found. This command needs segment data "
             "(roots import-segments); Ancestry does not provide any.")
        store.close()
        return 0
    rows = []
    for g in groups[: args.limit]:
        rows.append([
            f"chr{g.chrom}", f"{g.start_bp / 1e6:.1f}-{g.end_bp / 1e6:.1f}Mb",
            f"{g.cm:.1f}", g.size, g.side() or "-", g.support,
            ", ".join(g.names[:3])[:40],
        ])
    _table(rows, ["chr", "region", "cM", "n", "side", "support", "members"])
    conflicts = [g for g in groups if g.conflict]
    if conflicts:
        _out()
        _out(f"{len(conflicts)} group(s) mix maternal and paternal matches. Those "
             "segments sit on opposite copies of the chromosome and are not one "
             "ancestral line.")
    store.close()
    return 0


def cmd_paint(args) -> int:
    store = _open_store(args)
    gmap = _load_map(store, args)
    kit = _resolve_kit(store, args.kit)
    painted = tri_mod.paint(store, kit.id, gmap, min_cm=args.min_cm)
    cov = tri_mod.coverage(painted, gmap)
    _rule(f"chromosome coverage for {kit.display}")
    rows = []
    for chrom in AUTOSOMES:
        segs = painted.get(chrom, [])
        mat = sum(1 for s in segs if s.side == "maternal")
        pat = sum(1 for s in segs if s.side == "paternal")
        rows.append([f"chr{chrom}", len(segs), mat, pat, len(segs) - mat - pat])
    _table(rows, ["chr", "segments", "maternal", "paternal", "unassigned"])
    _out()
    _out(f"attributed to at least one match: {cov['any']:.0f} cM of "
         f"{cov['genome']:.0f} cM ({cov['any'] / cov['genome']:.0%})")
    _out(f"  maternal {cov['maternal']:.0f} cM, paternal {cov['paternal']:.0f} cM")
    _out(f"explained by nobody: {cov['unattributed']:.0f} cM "
         f"({cov['unattributed'] / cov['genome']:.0%}) -- this is where your "
         "undiscovered branches are")
    store.close()
    return 0


def cmd_predict(args) -> int:
    store = _open_store(args)
    gmap = _load_map(store, args)
    pred = predict_mod.predict(
        store, args.cm, gmap=gmap, min_cm=args.min_cm, iterations=args.iterations,
        segments=args.segments, ibd2_cm=args.ibd2, age_gap_years=args.age_gap,
        use_prior=not args.flat_prior,
    )
    _rule(f"relationships that could explain {args.cm:.0f} cM")
    rows = []
    for c in pred.top(args.limit):
        rows.append([
            c.rel.name[:40], f"{c.posterior:.1%}",
            f"{c.dist.mean:.0f}" if c.dist else "-",
            c.range_text(),
        ])
    _table(rows, ["relationship", "probability", "avg cM", "typical range"])
    _out()
    for n in pred.notes:
        _out(textwrap.fill("note: " + n, 78))
    if not args.flat_prior:
        _out()
        _out(textwrap.fill(
            "Probabilities are weighted by how many relatives of each type a "
            "person typically has, which is why distant relationships can "
            "outrank closer ones at low cM. Use --flat-prior to see the raw "
            "likelihoods instead.", 78))
    store.close()
    return 0


def cmd_check(args) -> int:
    store = _open_store(args)
    gmap = _load_map(store, args)
    rel = sim_mod.CATALOGUE.get(args.relationship)
    if not rel:
        raise SystemExit(
            f"unknown relationship {args.relationship!r}. "
            f"Known keys: {', '.join(sorted(sim_mod.CATALOGUE))}"
        )
    fit = predict_mod.fit_hypothesis(
        store, rel, args.cm, gmap=gmap, min_cm=args.min_cm, iterations=args.iterations
    )
    _rule(f"does {args.cm:.0f} cM fit '{rel.name}'?")
    _out(f"simulated average {fit.expected_cm:.0f} cM")
    _out(f"observed sits at the {fit.percentile:.0%} percentile")
    _out(f"verdict: {'plausible' if fit.plausible else 'implausible'}")
    _out(textwrap.fill(fit.comment, 78))
    store.close()
    return 0


def cmd_simulate(args) -> int:
    store = _open_store(args)
    gmap = _load_map(store, args)
    keys = [args.relationship] if args.relationship else list(sim_mod.CATALOGUE)
    rows = []
    for key in keys:
        rel = sim_mod.CATALOGUE.get(key)
        if not rel:
            continue
        d = sim_mod.cached_distribution(
            store, rel, iterations=args.iterations, min_cm=args.min_cm, gmap=gmap
        )
        lo, hi = d.range()
        rows.append([
            rel.name[:40], f"{d.mean:.0f}", f"{d.stdev:.0f}",
            f"{lo:.0f}-{hi:.0f}", f"{d.p_undetected:.0%}",
        ])
    _table(rows, ["relationship", "mean cM", "sd", "5-95%", "no match"])
    _out()
    _out(textwrap.fill(
        "'no match' is the share of simulated pairs sharing nothing detectable. "
        "It is the reason a missing relative is weak evidence at cousin level "
        "and no evidence at all further out.", 78))
    store.close()
    return 0


def cmd_import_tree(args) -> int:
    from .tree import gedcom as gedcom_mod

    store = _open_store(args)
    stats = gedcom_mod.import_gedcom(store, args.path, replace=args.replace)
    _out(f"{os.path.basename(args.path)}: "
         + ", ".join(f"{v} {k}" for k, v in stats.items()))
    store.close()
    return 0


def cmd_tree(args) -> int:
    from .tree import kinship

    store = _open_store(args)
    idx = kinship.load_tree(store)
    if args.action == "stats":
        n_ind, n_fam = idx.size()
        _out(f"{n_ind} individuals, {n_fam} families")
        surnames = {}
        for ind in idx.individuals.values():
            s = (ind.get("surname") or "?").strip()
            surnames[s] = surnames.get(s, 0) + 1
        top = sorted(surnames.items(), key=lambda kv: -kv[1])[:15]
        _table([[s, n] for s, n in top], ["surname", "count"])
    elif args.action == "find":
        hits = kinship.find_people(idx, args.query or "")
        _table([[x, n] for x, n in hits], ["xref", "name"])
    elif args.action == "relate":
        paths = kinship.relationship_between(idx, args.query, args.other)
        if not paths:
            _out("no connection found in the tree")
        for p in paths:
            _out(f"- {p.describe(idx)}")
    elif args.action == "gaps":
        root = args.query
        gaps = kinship.tree_gaps(idx, root)
        rows = [[g, k, t, f"{k / t:.0%}"] for g, (k, t) in sorted(gaps.items())]
        _table(rows, ["generation", "known", "possible", "complete"])
    store.close()
    return 0


def cmd_person(args) -> int:
    store = _open_store(args)
    store.upsert_person(
        args.label, sex=args.sex, birth_year=args.birth_year, tree_xref=args.tree_xref
    )
    rows = [[p["id"], p["label"], p["sex"] or "-", p["birth_year"] or "-",
             p["tree_xref"] or "-"] for p in store.people()]
    _table(rows, ["id", "label", "sex", "born", "tree xref"])
    store.close()
    return 0


def cmd_link(args) -> int:
    store = _open_store(args)
    m = store.match_by_id(args.match)
    if not m:
        raise SystemExit(f"no match with id {args.match}")
    store.db.execute("UPDATE match SET tree_xref=? WHERE id=?", (args.xref, args.match))
    store.commit()
    _out(f"linked match {args.match} ({m['name']}) to tree individual {args.xref}")
    _out("run 'roots investigate' to test whether the DNA agrees with the tree")
    store.close()
    return 0


def cmd_claim(args) -> int:
    store = _open_store(args)
    gmap = _load_map(store, args)
    m = store.match_by_id(args.match)
    if not m:
        raise SystemExit(f"no match with id {args.match}")
    rel = sim_mod.CATALOGUE.get(args.relationship)
    if not rel:
        raise SystemExit(
            f"unknown relationship {args.relationship!r}. "
            f"Known keys: {', '.join(sorted(sim_mod.CATALOGUE))}"
        )
    extra = json.loads(m["extra"]) if m["extra"] else {}
    extra["claimed_relationship"] = rel.key
    store.db.execute("UPDATE match SET extra=? WHERE id=?", (json.dumps(extra), args.match))
    store.commit()

    fit = predict_mod.fit_hypothesis(
        store, rel, m["total_cm"] or 0.0, gmap=gmap, min_cm=args.min_cm,
        iterations=args.iterations,
    )
    _out(f"recorded: {m['name']} claimed as {rel.name}")
    _out(f"testing that claim against {m['total_cm']:.0f} cM of shared DNA:")
    _out(f"  simulated average for {rel.name}: {fit.expected_cm:.0f} cM")
    _out(f"  {textwrap.fill(fit.comment, 74)}")
    store.close()
    return 0


def cmd_matches(args) -> int:
    store = _open_store(args)
    kit = _resolve_kit(store, args.kit)
    rows = []
    for m in store.matches(kit.id, source=args.source):
        if args.side and (m["side"] or "unassigned") != args.side:
            continue
        rows.append([
            m["id"], (m["name"] or "")[:30], f"{m['total_cm'] or 0:.0f}",
            m["seg_count"] or "-", m["side"] or "-",
            (m["predicted"] or "")[:22], m["tree_xref"] or "-",
        ])
        if len(rows) >= args.limit:
            break
    _table(rows, ["id", "name", "cM", "segs", "side", "site prediction", "tree"])
    store.close()
    return 0


def cmd_investigate(args) -> int:
    store = _open_store(args)
    gmap = _load_map(store, args)
    _warn_if_approximate(gmap)
    kit = _resolve_kit(store, args.kit)
    _out(f"investigating {kit.display} ...")
    findings = hypo.investigate(
        store, kit.id, gmap, iterations=args.iterations, min_cm=args.min_cm
    )
    if not findings:
        _out("nothing to report yet. Import matches, sides and a tree first.")
        store.close()
        return 0
    by_kind = {}
    for f in findings:
        by_kind.setdefault(f.kind, []).append(f)
    for kind, group in sorted(by_kind.items(), key=lambda kv: -max(f.score for f in kv[1])):
        _rule(kind.replace("-", " "))
        for f in group[: args.limit]:
            _out(f"* {f.subject}")
            _out(textwrap.fill(f.summary, 76, initial_indent="    ",
                               subsequent_indent="    "))
    plan = hypo.research_plan(findings)
    if plan:
        _rule("what to do next")
        for i, step in enumerate(plan, 1):
            _out(textwrap.fill(f"{i}. {step}", 76, subsequent_indent="   "))
    store.close()
    return 0


def cmd_source(args) -> int:
    from . import evidence as ev

    store = _open_store(args)
    if args.action == "add":
        if not args.title:
            raise SystemExit("--title is required when adding a source")
        sid = ev.add_source(
            store, args.title, repository=args.repository,
            record_type=args.type, reference=args.reference, url=args.url,
            cost=args.cost, notes=args.notes,
        )
        _out(f"source {sid} recorded: {args.title}")
        _out("link it to a person with: roots evidence add --source "
             f"{sid} --xref <xref> --type birth --claim '...'")
    else:
        rows = [[r["id"], (r["title"] or "")[:38], r["repository"] or "-",
                 r["record_type"] or "-", r["accessed"] or "-",
                 f"{r['cost']:.2f}" if r["cost"] else "-"]
                for r in ev.sources(store)]
        _table(rows, ["id", "title", "repository", "type", "accessed", "cost"])
        total = sum(r[5] != "-" and float(r[5]) or 0 for r in rows)
        if total:
            _out(f"\ntotal spent on records: {total:.2f}")
    store.close()
    return 0


def cmd_evidence(args) -> int:
    from . import evidence as ev
    from .tree import kinship

    store = _open_store(args)
    if args.action == "add":
        for required in ("source", "xref", "claim"):
            if getattr(args, required) is None:
                raise SystemExit(f"--{required} is required when adding evidence")
        eid = ev.add_evidence(
            store, args.source, args.xref, args.type, args.claim,
            supports=not args.contradicts, confidence=args.confidence,
            notes=args.notes,
        )
        verb = "contradicts" if args.contradicts else "supports"
        _out(f"evidence {eid} recorded: {verb} the tree for {args.xref}")
        if args.contradicts:
            _out("a contradiction between documents is exactly what DNA is "
                 "qualified to settle -- check 'roots investigate' for whether "
                 "the sharing agrees with the tree here")
    elif args.action == "conflicts":
        rows = [[r["subject_xref"], r["claim_type"] or "-",
                 (r["claim"] or "")[:40], (r["title"] or "")[:28]]
                for r in ev.contradictions(store)]
        _table(rows, ["person", "claim type", "claim", "source"])
    else:
        idx = kinship.load_tree(store)
        cov = ev.coverage(store, idx)
        _rule("documentary coverage")
        _out(cov.summary())
        if cov.total_cost:
            _out(f"spent so far: {cov.total_cost:.2f}")
        if cov.individuals and cov.fraction < 0.5:
            _out(textwrap.fill(
                "Most of this tree rests on assertion rather than documents. "
                "That matters before using DNA to 'confirm' any of it: DNA can "
                "only confirm a relationship you have stated correctly.", 78))
        if args.xref:
            _rule(f"evidence for {idx.name(args.xref)}")
            rows = [[r["claim_type"] or "-", (r["claim"] or "")[:40],
                     "supports" if r["supports"] else "CONTRADICTS",
                     (r["title"] or "")[:28]]
                    for r in ev.evidence_for(store, args.xref)]
            _table(rows, ["claim type", "claim", "stance", "source"])
    store.close()
    return 0


def cmd_research(args) -> int:
    from . import evidence as ev
    from .tree import kinship

    store = _open_store(args)
    idx = kinship.load_tree(store)
    root = args.root
    if not root:
        row = store.db.execute(
            "SELECT tree_xref FROM person WHERE tree_xref IS NOT NULL LIMIT 1"
        ).fetchone()
        root = row["tree_xref"] if row else None
    kit_id = None
    if args.kit:
        kit_id = _resolve_kit(store, args.kit).id

    if not root and not kit_id:
        raise SystemExit(
            "nothing to work from. Link yourself into the tree with "
            "'roots person --label self --tree-xref <xref>', or pass --kit to "
            "get tasks from your DNA clusters."
        )

    tasks = ev.generate_tasks(store, root=root, idx=idx, kit_id=kit_id,
                              limit=args.limit)
    if not tasks:
        _out("no research tasks generated. Import a tree and link yourself "
             "into it, or run 'roots cluster' first.")
        store.close()
        return 0

    if args.budget:
        from . import budget as budget_mod

        costs = budget_mod.load_costs(args.cost_file) if args.cost_file else None
        led = budget_mod.ledger(store)
        plan = budget_mod.plan_within_budget(tasks, args.budget, costs)
        _rule(f"what to do next, within £{args.budget:.2f}")
        _out(plan.summary())
        if led.spent:
            _out(f"(you have already spent £{led.spent:.2f} on {led.records} records)")
        if plan.subscription_taken:
            repo, fee = plan.subscription_taken
            _out()
            _out(textwrap.fill(
                f"The plan includes one month of {repo} at £{fee:.2f}. That fee "
                "covers every subscription task below, so do them all inside the "
                "same month rather than spreading them out.", 78))
        _out()
        for task, cost in plan.selected:
            price = "free" if cost.kind == "free" else (
                "included" if cost.kind == "subscription" else f"£{cost.amount:.2f}")
            _out(f"{task.line()}\n        cost: {price}"
                 + (f" -- {cost.note}" if cost.note else ""))
            if task.rationale and args.why:
                _out(textwrap.fill(task.rationale, 72,
                                   initial_indent="        why: ",
                                   subsequent_indent="             "))
            _out()
        if plan.deferred:
            _rule("deferred until there is more budget")
            for task, cost in plan.deferred[:10]:
                _out(f"  £{cost.amount:6.2f}  {task.subject}: {task.question}")
        store.close()
        return 0

    _rule("what to look up next")
    _out(textwrap.fill(
        "Ordered by how much each would unlock per pound spent. Costs are "
        "bands, not prices -- see docs/RECORD_SOURCES.md for current figures. "
        "Pass --budget to cost the plan and fit it to a limit.",
        78))
    _out()
    for t in tasks:
        _out(t.line())
        if t.rationale and args.why:
            _out(textwrap.fill(t.rationale, 72, initial_indent="        why: ",
                               subsequent_indent="             "))
        _out()
    store.close()
    return 0


def cmd_import_citations(args) -> int:
    from . import evidence as ev

    store = _open_store(args)
    stats = ev.import_gedcom_sources(store, args.path)
    _out(f"{os.path.basename(args.path)}: {stats['sources']} sources, "
         f"{stats['citations']} citations")
    _out("these record what has already been looked up, which is what stops "
         "you buying the same certificate twice")
    store.close()
    return 0


def cmd_budget(args) -> int:
    from . import budget as budget_mod

    store = _open_store(args)
    if args.set is not None:
        budget_mod.set_cap(store, args.set if args.set > 0 else None)
    led = budget_mod.ledger(store)
    _rule("budget")
    _out(led.summary())
    rows = [[r["id"], (r["title"] or "")[:40], r["repository"] or "-",
             f"{r['cost']:.2f}" if r["cost"] else "0.00"]
            for r in store.db.execute(
                "SELECT * FROM source WHERE cost IS NOT NULL AND cost > 0 "
                "ORDER BY id DESC LIMIT 20")]
    if rows:
        _out()
        _table(rows, ["id", "record", "repository", "cost"])
    _out()
    _out(textwrap.fill(
        f"Cost estimates were last checked on {budget_mod.VERIFIED_ON}; see "
        "docs/RECORD_SOURCES.md. Override them with --cost-file if they have "
        "aged.", 78))
    store.close()
    return 0


def cmd_auto(args) -> int:
    """Run every analysis step that is free and deterministic."""
    from . import budget as budget_mod
    from . import evidence as ev
    from .tree import kinship

    store = _open_store(args)
    gmap = _load_map(store, args)
    done: List[str] = []
    skipped: List[str] = []

    def step(name: str, ok: bool, detail: str = "") -> None:
        (done if ok else skipped).append(f"{name}{': ' + detail if detail else ''}")
        _out(f"  {'ok  ' if ok else 'skip'}  {name}" + (f" -- {detail}" if detail else ""))

    _rule("automated pipeline")
    _out("Every step here is free, offline and repeatable. Nothing that costs")
    _out("money or needs a login happens without you.")
    _out()

    for person in {k.person for k in store.kits() if k.person}:
        kits = [k for k in store.kits_for_person(person) if k.vendor != "merged"]
        merged_label = f"{person}-merged"
        if len(kits) > 1 and not store.kit(merged_label):
            try:
                _new, stats = qc_mod.merge_kits(
                    store, [k.id for k in kits], merged_label, person,
                    build=kits[0].build or "37")
                carried = ingest_mod.carry_matches(
                    store, [k.id for k in kits], _new)
                step(f"merge kits for {person}", True,
                     f"{stats['positions']:,} positions"
                     + (f", {carried['matches']} matches carried"
                        if carried["matches"] else ""))
            except Exception as exc:
                step(f"merge kits for {person}", False, str(exc))
        elif store.kit(merged_label):
            step(f"merge kits for {person}", True, "already merged")

    child = store.kit(args.self_kit) if args.self_kit else _best_kit(store, "self")
    parent = store.kit(args.parent_kit) if args.parent_kit else _best_kit(store, args.role)
    if not child:
        _out("\nno kit for 'self'. Import one with roots import-dna --person self")
        store.close()
        return 1

    if parent:
        cmp = relate_mod.compare(store, child.id, parent.id, gmap)
        verdict = relate_mod.classify_close(cmp, gmap)
        ok = verdict.label == "parent/child"
        step("verify parent-child", ok, f"{cmp.total_cm:.0f} cM, {verdict.label}")
        if ok:
            stats = phase_mod.phase_against_parent(
                store, child.id, parent.id, parent_role=args.role,
                child_sex=child.inferred_sex)
            step("phase", True,
                 f"{stats.phased_fraction:.0%} of sites, "
                 f"{stats.mendelian_rate:.3%} Mendelian errors")
            report = sides_mod.assign_from_parent_list(
                store, child.id, parent.id, parent_role=args.role)
            step("assign sides", bool(report.calls), report.summary())
    else:
        step("verify parent-child", False, f"no kit for '{args.role}'")

    result = cluster_mod.cluster_matches(store, child.id)
    if result.clusters:
        cluster_mod.persist(store, child.id, result)
    step("cluster", bool(result.clusters), result.summary())

    groups = tri_mod.overlap_groups(store, child.id, gmap)
    step("triangulate", bool(groups), f"{len(groups)} shared regions")

    findings = hypo.investigate(store, child.id, gmap, iterations=args.iterations)
    step("investigate", bool(findings), f"{len(findings)} findings")

    idx = kinship.load_tree(store)
    row = store.db.execute(
        "SELECT tree_xref FROM person WHERE tree_xref IS NOT NULL LIMIT 1").fetchone()
    root = row["tree_xref"] if row else None
    tasks = ev.generate_tasks(store, root=root, idx=idx, kit_id=child.id)
    step("plan research", bool(tasks), f"{len(tasks)} tasks")

    if tasks and args.budget:
        plan = budget_mod.plan_within_budget(tasks, args.budget)
        step("fit to budget", True, plan.summary())

    if args.report:
        from . import report as report_mod
        path = report_mod.build_report(
            store, child.id, gmap, args.report, findings=findings,
            cluster_result=result, iterations=args.iterations)
        step("report", True, path)

    _rule("summary")
    _out(f"{len(done)} steps completed, {len(skipped)} skipped")
    if skipped:
        _out()
        _out("skipped steps need data you have not imported yet:")
        for s in skipped:
            _out(f"  - {s}")
    store.close()
    return 0


def _best_kit(store: Store, person: str):
    """Prefer a person's merged kit, falling back to any kit they have."""
    kits = store.kits_for_person(person)
    if not kits:
        return None
    merged = [k for k in kits if k.vendor == "merged"]
    return merged[0] if merged else kits[0]


def cmd_report(args) -> int:
    from . import report as report_mod

    store = _open_store(args)
    gmap = _load_map(store, args)
    kit = _resolve_kit(store, args.kit)
    findings = hypo.investigate(store, kit.id, gmap, iterations=args.iterations,
                                persist=False)
    result = None
    try:
        result = cluster_mod.cluster_matches(store, kit.id)
    except Exception:
        result = None
    out = args.out or f"roots-report-{kit.label}.html"
    path = report_mod.build_report(
        store, kit.id, gmap, out, findings=findings, cluster_result=result,
        iterations=args.iterations,
    )
    _out(f"wrote {path}")
    _out("the report contains identifiable genetic information about you and "
         "your relatives. Treat it accordingly.")
    store.close()
    return 0


def cmd_demo(args) -> int:
    from . import demo as demo_mod

    files = demo_mod.generate(args.out)
    _out(f"wrote a synthetic family to {args.out}:")
    for name, path in sorted(files.items()):
        _out(f"  {name:24s} {os.path.basename(path)}")
    _out()
    _out("This data is simulated, not real. Every parent-child link, segment and "
         "shared-cM figure was produced by actually simulating inheritance, so "
         "the whole pipeline can be exercised against known truth.")
    return 0


def cmd_guide(args) -> int:
    _out(textwrap.dedent("""\
        Recommended order of operations
        ===============================

        You have your own DNA and your mother's. That combination is worth more
        than the sum of its parts, because it splits every match you will ever
        have into maternal and paternal.

        1.  roots init
        2.  roots import-dna you-23andme.txt     --person self   --label self-23
            roots import-dna you-ancestry.txt    --person self   --label self-anc
            roots import-dna mum-23andme.txt     --person mother --label mum-23
        3.  roots merge-kits --person self       (union of both vendors' SNPs)
            roots merge-kits --person mother
        4.  roots compare self-merged mum-merged
              Confirms from the data alone that she is your mother.
        5.  roots phase --child self-merged --parent mum-merged \\
                --export-other father-inferred.txt
              Splits your genome into her half and your father's half, and
              writes the paternal half as an uploadable pseudo-kit.
        6.  Export match lists for BOTH of you from the same site, then:
            roots import-matches you-matches.csv    --kit self-merged
            roots import-matches mum-matches.csv    --kit mum-merged
        7.  roots sides --child self-merged --parent mum-merged
              Everyone on her list is maternal. Everyone absent from it, above
              the reporting floor, is paternal.
        8.  roots import-segments you-segments.csv --kit self-merged   (if available)
            roots import-icw      you-icw.csv      --kit self-merged
        9.  roots cluster --kit self-merged
            roots triangulate --kit self-merged
        10. roots import-tree yourtree.ged
            roots import-citations yourtree.ged
            roots person --label self --tree-xref I1
            roots link --match <id> --xref <xref>     for matches you know
        11. roots investigate --kit self-merged
            roots research --kit self-merged --why
              Turns the walls in your tree into a list of records to order,
              cheapest and highest-yield first.
        12. roots report --kit self-merged --out report.html

        docs/DATA_SOURCES.md   which testing sites to export from, and which
                               of them can triangulate at all.
        docs/RECORD_SOURCES.md which record providers cost what, and in what
                               order to spend money on certificates.
        """))
    return 0


# ---------------------------------------------------------------------------
# argument parsing


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="roots",
        description="Genealogy toolkit: raw DNA, match lists, family trees, and "
                    "the reasoning that connects them. Entirely offline.",
    )
    p.add_argument("--project", default=DEFAULT_PROJECT, help="project database file")
    p.add_argument("--map", help="path to a recombination map (file or directory)")
    p.add_argument("--build", default=None, help="reference assembly, 37 or 38")
    sub = p.add_subparsers(dest="command", required=True)

    sp = sub.add_parser("init", help="create a project")
    sp.add_argument("--build", default="37")
    sp.add_argument("--map")
    sp.set_defaults(func=cmd_init)

    sp = sub.add_parser("guide", help="show the recommended workflow")
    sp.set_defaults(func=cmd_guide)

    sp = sub.add_parser("import-dna", help="import a raw DNA download")
    sp.add_argument("path")
    sp.add_argument("--person", required=True, help="who this kit belongs to")
    sp.add_argument("--label", help="name for this kit (defaults to the filename)")
    sp.add_argument("--vendor", help="override vendor detection")
    sp.add_argument("--build-override", help="override assembly detection")
    sp.set_defaults(func=cmd_import_dna)

    sp = sub.add_parser("kits", help="list imported kits")
    sp.set_defaults(func=cmd_kits)

    sp = sub.add_parser("qc", help="quality-control a kit")
    sp.add_argument("--kit")
    sp.set_defaults(func=cmd_qc)

    sp = sub.add_parser("merge-kits", help="combine one person's kits from several vendors")
    sp.add_argument("--person", required=True)
    sp.add_argument("--label")
    sp.set_defaults(func=cmd_merge_kits)

    sp = sub.add_parser("compare", help="find shared segments between two kits")
    sp.add_argument("kit_a")
    sp.add_argument("kit_b")
    sp.add_argument("--min-cm", type=float, default=7.0)
    sp.add_argument("--min-snps", type=int, default=500)
    sp.add_argument("--save", action="store_true", help="store the segments found")
    sp.set_defaults(func=cmd_compare)

    sp = sub.add_parser("phase", help="split a kit into maternal and paternal halves")
    sp.add_argument("--child", required=True)
    sp.add_argument("--parent", required=True)
    sp.add_argument("--role", choices=["mother", "father"], default="mother")
    sp.add_argument("--export-other", help="write the untested parent's haplotype here")
    sp.add_argument("--export-same", help="write the tested parent's haplotype here")
    sp.add_argument("--min-cm", type=float, default=7.0)
    sp.add_argument("--skip-verify", action="store_true")
    sp.set_defaults(func=cmd_phase)

    sp = sub.add_parser("import-matches", help="import a match list")
    sp.add_argument("path")
    sp.add_argument("--kit", required=True)
    sp.add_argument("--source", help="platform the list came from")
    sp.set_defaults(func=cmd_import_matches)

    sp = sub.add_parser("import-segments", help="import chromosome-browser data")
    sp.add_argument("path")
    sp.add_argument("--kit", required=True)
    sp.add_argument("--source")
    sp.set_defaults(func=cmd_import_segments)

    sp = sub.add_parser("import-icw", help="import shared-match (in-common-with) pairs")
    sp.add_argument("path")
    sp.add_argument("--kit", required=True)
    sp.add_argument("--kind", choices=["icw", "triangulated"], default="icw")
    sp.set_defaults(func=cmd_import_icw)

    sp = sub.add_parser("sides", help="sort matches maternal/paternal using a parent's list")
    sp.add_argument("--child", required=True)
    sp.add_argument("--parent", required=True)
    sp.add_argument("--role", choices=["mother", "father"], default="mother")
    sp.add_argument("--source")
    sp.add_argument("--floor", type=float, default=20.0,
                    help="minimum cM at which absence from the parent list is evidence")
    sp.add_argument("--floor-multiple", type=float, default=2.0,
                    help="safety margin above the parent list's observed cutoff; "
                         "set to 0 if you know the list is complete")
    sp.add_argument("--segments", action="store_true", help="also confirm via segments")
    sp.add_argument("--verbose", action="store_true")
    sp.set_defaults(func=cmd_sides)

    sp = sub.add_parser("cluster", help="group matches into ancestral lines")
    sp.add_argument("--kit", required=True)
    sp.add_argument("--min-cm", type=float, default=40.0)
    sp.add_argument("--max-cm", type=float, default=400.0)
    sp.add_argument("--cutoff", type=float, default=0.80)
    sp.add_argument("--min-size", type=int, default=2)
    sp.add_argument("--source")
    sp.set_defaults(func=cmd_cluster)

    sp = sub.add_parser("triangulate", help="find regions several matches share")
    sp.add_argument("--kit", required=True)
    sp.add_argument("--min-cm", type=float, default=7.0)
    sp.add_argument("--min-members", type=int, default=2)
    sp.add_argument("--max-match-cm", type=float, default=2000.0,
                    help="ignore relatives closer than this; they overlap "
                         "everything and locate nothing")
    sp.add_argument("--limit", type=int, default=30)
    sp.add_argument("--source")
    sp.set_defaults(func=cmd_triangulate)

    sp = sub.add_parser("paint", help="chromosome coverage by side")
    sp.add_argument("--kit", required=True)
    sp.add_argument("--min-cm", type=float, default=7.0)
    sp.set_defaults(func=cmd_paint)

    sp = sub.add_parser("predict", help="rank relationships for an amount of shared DNA")
    sp.add_argument("cm", type=float)
    sp.add_argument("--segments", type=int)
    sp.add_argument("--ibd2", type=float, help="fully-identical cM, if known")
    sp.add_argument("--age-gap", type=float, help="match's birth year minus yours")
    sp.add_argument("--flat-prior", action="store_true")
    sp.add_argument("--min-cm", type=float, default=7.0)
    sp.add_argument("--iterations", type=int, default=1500)
    sp.add_argument("--limit", type=int, default=12)
    sp.set_defaults(func=cmd_predict)

    sp = sub.add_parser("check", help="test one specific relationship hypothesis")
    sp.add_argument("relationship")
    sp.add_argument("--cm", type=float, required=True)
    sp.add_argument("--min-cm", type=float, default=7.0)
    sp.add_argument("--iterations", type=int, default=1500)
    sp.set_defaults(func=cmd_check)

    sp = sub.add_parser("simulate", help="show simulated sharing distributions")
    sp.add_argument("relationship", nargs="?")
    sp.add_argument("--iterations", type=int, default=1500)
    sp.add_argument("--min-cm", type=float, default=7.0)
    sp.set_defaults(func=cmd_simulate)

    sp = sub.add_parser("import-tree", help="import a GEDCOM file")
    sp.add_argument("path")
    sp.add_argument("--replace", action="store_true")
    sp.set_defaults(func=cmd_import_tree)

    sp = sub.add_parser("tree", help="query the family tree")
    sp.add_argument("action", choices=["stats", "find", "relate", "gaps"])
    sp.add_argument("query", nargs="?")
    sp.add_argument("other", nargs="?")
    sp.set_defaults(func=cmd_tree)

    sp = sub.add_parser("person", help="record who a kit belongs to")
    sp.add_argument("--label", required=True)
    sp.add_argument("--sex")
    sp.add_argument("--birth-year", type=int)
    sp.add_argument("--tree-xref")
    sp.set_defaults(func=cmd_person)

    sp = sub.add_parser("link", help="link a match to a person in the tree")
    sp.add_argument("--match", type=int, required=True)
    sp.add_argument("--xref", required=True)
    sp.set_defaults(func=cmd_link)

    sp = sub.add_parser("claim", help="record a believed relationship and test it")
    sp.add_argument("--match", type=int, required=True)
    sp.add_argument("--relationship", required=True)
    sp.add_argument("--min-cm", type=float, default=7.0)
    sp.add_argument("--iterations", type=int, default=1500)
    sp.set_defaults(func=cmd_claim)

    sp = sub.add_parser("matches", help="list imported matches")
    sp.add_argument("--kit", required=True)
    sp.add_argument("--side")
    sp.add_argument("--source")
    sp.add_argument("--limit", type=int, default=40)
    sp.set_defaults(func=cmd_matches)

    sp = sub.add_parser("investigate", help="run every check and report findings")
    sp.add_argument("--kit", required=True)
    sp.add_argument("--iterations", type=int, default=1500)
    sp.add_argument("--min-cm", type=float, default=7.0)
    sp.add_argument("--limit", type=int, default=8)
    sp.set_defaults(func=cmd_investigate)


    sp = sub.add_parser("source", help="record documents you have looked at")
    sp.add_argument("action", choices=["add", "list"], nargs="?", default="list")
    sp.add_argument("--title")
    sp.add_argument("--repository", help="GRO, FindMyPast, ScotlandsPeople, ...")
    sp.add_argument("--type", help="birth | marriage | death | census | parish | will")
    sp.add_argument("--reference", help="index reference, piece number, volume/page")
    sp.add_argument("--url")
    sp.add_argument("--cost", type=float)
    sp.add_argument("--notes")
    sp.set_defaults(func=cmd_source)

    sp = sub.add_parser("evidence", help="what a document asserts about a person")
    sp.add_argument("action", choices=["add", "list", "conflicts"], nargs="?",
                    default="list")
    sp.add_argument("--source", type=int)
    sp.add_argument("--xref")
    sp.add_argument("--type", default="identity")
    sp.add_argument("--claim")
    sp.add_argument("--contradicts", action="store_true",
                    help="the document disagrees with the tree as it stands")
    sp.add_argument("--confidence", default="direct",
                    choices=["direct", "indirect", "negative"])
    sp.add_argument("--notes")
    sp.set_defaults(func=cmd_evidence)

    sp = sub.add_parser("research", help="prioritised list of records to look up")
    sp.add_argument("--root", help="your xref in the tree")
    sp.add_argument("--kit", help="also generate tasks from DNA clusters")
    sp.add_argument("--limit", type=int, default=25)
    sp.add_argument("--why", action="store_true", help="explain each suggestion")
    sp.add_argument("--budget", type=float,
                    help="only propose what fits this many pounds")
    sp.add_argument("--cost-file", help="JSON overriding the built-in price table")
    sp.set_defaults(func=cmd_research)

    sp = sub.add_parser("import-citations",
                        help="pull source citations out of a GEDCOM")
    sp.add_argument("path")
    sp.set_defaults(func=cmd_import_citations)

    sp = sub.add_parser("budget", help="what you have spent, and the cap")
    sp.add_argument("--set", type=float, help="set a spending cap (0 clears it)")
    sp.set_defaults(func=cmd_budget)

    sp = sub.add_parser("auto", help="run every free, deterministic step")
    sp.add_argument("--self-kit", help="kit for you (default: your merged kit)")
    sp.add_argument("--parent-kit", help="kit for the tested parent")
    sp.add_argument("--role", choices=["mother", "father"], default="mother")
    sp.add_argument("--budget", type=float, help="also fit a research plan to a budget")
    sp.add_argument("--report", help="write an HTML report here too")
    sp.add_argument("--iterations", type=int, default=1500)
    sp.set_defaults(func=cmd_auto)

    sp = sub.add_parser("report", help="write an HTML report")
    sp.add_argument("--kit", required=True)
    sp.add_argument("--out")
    sp.add_argument("--iterations", type=int, default=1500)
    sp.set_defaults(func=cmd_report)

    sp = sub.add_parser("demo", help="generate a synthetic family to try the tools on")
    sp.add_argument("--out", default="demo-data")
    sp.set_defaults(func=cmd_demo)

    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except BrokenPipeError:
        return 0
    except KeyboardInterrupt:
        _out("interrupted")
        return 130


if __name__ == "__main__":
    sys.exit(main())

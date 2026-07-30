"""Assemble everything the project knows into one page you can keep.

A genetic genealogy investigation is spread across a dozen half-answers: a
match list, a pile of segments, a clustering run, a tree with holes in it.
None of them mean much alone, and the interesting conclusions live in the
places where two of them disagree.  This module is where they are put side by
side.

Three constraints shape the output, and all three are deliberate.

*It is one file, and it works offline forever.*  The report contains a named
person's genotype statistics, their relatives' names, and a map of which
stretches of their chromosomes came from whom.  Nothing in it may phone home:
no fonts, no CDN, no analytics, not one URL.  Everything -- stylesheet,
chromosome map, cluster matrix -- is inlined, so the file still renders
correctly on a laptop with no network in twenty years' time.  A consequence
worth stating: this file is as sensitive as the raw data it came from, and
the header says so.

*Every string that came from a file is escaped.*  Match names, surnames and
cluster labels arrive from vendor CSV exports and GEDCOM files, which is to
say from arbitrary text nobody validated.  They go through ``html.escape``
without exception, including inside SVG ``<title>`` tooltips.

*Nothing here raises on missing data.*  The commonest moment to want a report
is right after importing the first kit, when there are no matches, no
segments, no clusters and no tree.  Every section degrades to an honest
sentence about what is absent rather than a traceback, because "you have not
imported a match list yet" is itself useful output.

The two charts are the ones that carry information no table can.

The chromosome map answers "where did my DNA come from", and more usefully
its inverse.  Segments are packed into lanes within each chromosome track by
a greedy interval assignment -- the first lane whose previous segment ends
before this one starts -- so overlapping matches sit above one another
instead of hiding each other.  Lanes are capped; anything past the cap is
drawn faintly in a final overflow lane and counted in the caption, because
silently dropping segments from a coverage picture would be a lie.  The gaps
in the picture are the point: they are the DNA no living match explains.

The cluster matrix is the Leeds picture.  Matches are ordered by cluster and
a cell is filled where two of them match each other, so shared ancestry shows
up as solid blocks on the diagonal.  Off-diagonal ink is a cluster boundary
that may not be real.
"""

from __future__ import annotations

import html
import json
import os
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .dna.predict import predict
from .dna.qc import qc_kit
from .genome import AUTOSOMES, CHROM_BP, GeneticMap
from .hypothesis import Finding, research_plan
from .matches.cluster import ClusterResult
from .matches.sides import side_summary
from .matches.triangulate import (
    OverlapGroup,
    PaintedSegment,
    coverage,
    overlap_groups,
    paint,
)
from .store import Store

PROJECT_NAME = "roots"

#: Sides are an identity encoding, not a magnitude, so they take the first
#: three categorical slots (blue / orange / aqua -- the three that clear the
#: all-pairs colour-vision gates in both modes).  "Unassigned" is not a
#: fourth identity: it is the absence of one, so it wears the muted ink and
#: never competes with a real side for attention.
SIDE_ORDER: Tuple[str, ...] = ("maternal", "paternal", "both", "unassigned")
SIDE_LABEL: Dict[str, str] = {
    "maternal": "Maternal",
    "paternal": "Paternal",
    "both": "Both sides",
    "unassigned": "Unassigned",
}

#: Lanes per chromosome track before segments start overlaying each other.
MAX_LANES = 8

#: Largest cluster matrix we will draw in full, in cells per side.
MAX_MATRIX = 120

# ---------------------------------------------------------------------------
# stylesheet
#
# Light and dark are both *selected*, not flipped: each mode names its own
# step of the same hues, chosen against that mode's surface.  Everything below
# is written against roles (--series-maternal, --grid, --ink-2) so the two
# blocks at the top are the only place a colour appears twice.

_CSS = """
:root {
  color-scheme: light;
  --plane: #f9f9f7;
  --surface: #fcfcfb;
  --ink: #0b0b0b;
  --ink-2: #52514e;
  --muted: #898781;
  --grid: #e1e0d9;
  --axis: #c3c2b7;
  --border: rgba(11, 11, 11, 0.10);
  --shadow: rgba(11, 11, 11, 0.05);
  --side-maternal: #2a78d6;
  --side-paternal: #eb6834;
  --side-both: #1baf7a;
  --side-unassigned: #898781;
  --c1: #2a78d6; --c2: #eb6834; --c3: #1baf7a; --c4: #eda100;
  --c5: #e87ba4; --c6: #008300; --c7: #4a3aa7; --c8: #e34948;
  --c-other: #898781;
  --good: #0ca30c;
  --warning: #fab219;
  --serious: #ec835a;
  --critical: #d03b3b;
  --wash: rgba(42, 120, 214, 0.08);
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --plane: #0d0d0d;
    --surface: #1a1a19;
    --ink: #ffffff;
    --ink-2: #c3c2b7;
    --muted: #898781;
    --grid: #2c2c2a;
    --axis: #383835;
    --border: rgba(255, 255, 255, 0.10);
    --shadow: rgba(0, 0, 0, 0.4);
    --side-maternal: #3987e5;
    --side-paternal: #d95926;
    --side-both: #199e70;
    --side-unassigned: #898781;
    --c1: #3987e5; --c2: #d95926; --c3: #199e70; --c4: #c98500;
    --c5: #d55181; --c6: #008300; --c7: #9085e9; --c8: #e66767;
    --c-other: #898781;
    --wash: rgba(57, 135, 229, 0.14);
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --plane: #0d0d0d;
  --surface: #1a1a19;
  --ink: #ffffff;
  --ink-2: #c3c2b7;
  --muted: #898781;
  --grid: #2c2c2a;
  --axis: #383835;
  --border: rgba(255, 255, 255, 0.10);
  --shadow: rgba(0, 0, 0, 0.4);
  --side-maternal: #3987e5;
  --side-paternal: #d95926;
  --side-both: #199e70;
  --side-unassigned: #898781;
  --c1: #3987e5; --c2: #d95926; --c3: #199e70; --c4: #c98500;
  --c5: #d55181; --c6: #008300; --c7: #9085e9; --c8: #e66767;
  --c-other: #898781;
  --wash: rgba(57, 135, 229, 0.14);
}

* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0;
  padding: 0 0 4rem;
  background: var(--plane);
  color: var(--ink);
  font: 15px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif;
  overflow-x: hidden;
}
.wrap { max-width: 68rem; margin: 0 auto; padding: 1.5rem 1rem; }

h1 { font-size: 1.6rem; line-height: 1.2; margin: 0 0 .35rem; font-weight: 650; }
h2 {
  font-size: 1.05rem; font-weight: 650; letter-spacing: .01em;
  margin: 0 0 .2rem;
}
h3 { font-size: .92rem; font-weight: 650; margin: 1.1rem 0 .4rem; }
p { margin: .45rem 0; }
a { color: inherit; }

.sub { color: var(--ink-2); font-size: .9rem; margin: 0; }
.muted { color: var(--muted); }
.small { font-size: .82rem; }
.mono { font-variant-numeric: tabular-nums; }

header.page {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: 12px;
  padding: 1.1rem 1.15rem;
  margin-bottom: 1.1rem;
}
header.page dl {
  display: grid;
  grid-template-columns: auto 1fr;
  gap: .2rem .9rem;
  margin: .8rem 0 0;
  font-size: .88rem;
}
header.page dt { color: var(--muted); }
header.page dd { margin: 0; color: var(--ink-2); }

.privacy {
  margin-top: .9rem;
  border-left: 3px solid var(--critical);
  background: var(--wash);
  border-radius: 0 8px 8px 0;
  padding: .6rem .8rem;
  font-size: .86rem;
  color: var(--ink-2);
}
.privacy strong { color: var(--ink); }

section.card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: 12px;
  padding: 1.05rem 1.15rem 1.15rem;
  margin: 0 0 1.1rem;
}
section.card > .lede { color: var(--ink-2); font-size: .86rem; margin: 0 0 .8rem; }

.scroll { overflow-x: auto; -webkit-overflow-scrolling: touch; max-width: 100%; }

table { border-collapse: collapse; width: 100%; font-size: .86rem; }
caption {
  caption-side: bottom; text-align: left; color: var(--muted);
  font-size: .8rem; padding-top: .5rem;
}
th, td {
  text-align: left; padding: .38rem .55rem; border-bottom: 1px solid var(--grid);
  white-space: nowrap; vertical-align: top;
}
th {
  color: var(--muted); font-weight: 600; font-size: .78rem;
  text-transform: uppercase; letter-spacing: .04em;
}
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
td.wrap-cell { white-space: normal; min-width: 14rem; }
tbody tr:last-child td { border-bottom: none; }

.stats { display: flex; flex-wrap: wrap; gap: .5rem; margin: .2rem 0 .9rem; }
.stat {
  flex: 1 1 8.5rem; min-width: 8.5rem;
  border: 1px solid var(--border); border-radius: 10px;
  padding: .55rem .7rem; background: var(--plane);
}
.stat .k {
  display: flex; align-items: center; gap: .4rem;
  color: var(--muted); font-size: .74rem; text-transform: uppercase;
  letter-spacing: .04em;
}
.stat .v { font-size: 1.3rem; font-weight: 650; line-height: 1.15; }
.stat .n { color: var(--ink-2); font-size: .78rem; }

.dot { width: 9px; height: 9px; border-radius: 50%; flex: 0 0 9px; display: inline-block; }
.swatch {
  width: 11px; height: 11px; border-radius: 3px; display: inline-block;
  vertical-align: -1px; margin-right: .35rem;
}

.legend {
  display: flex; flex-wrap: wrap; gap: .35rem .95rem;
  margin: .55rem 0 .2rem; font-size: .8rem; color: var(--ink-2);
}
.legend span { display: inline-flex; align-items: center; gap: .35rem; }

ul.notes { margin: .4rem 0 0; padding-left: 1.05rem; }
ul.notes li { margin: .2rem 0; color: var(--ink-2); font-size: .86rem; }

ul.caution { list-style: none; margin: .5rem 0 0; padding: 0; }
ul.caution li {
  display: flex; gap: .5rem; align-items: flex-start;
  border-left: 3px solid var(--warning);
  background: var(--plane);
  border-radius: 0 8px 8px 0;
  padding: .4rem .65rem; margin: .3rem 0; font-size: .85rem; color: var(--ink-2);
}
ul.caution li.crit { border-left-color: var(--critical); }
ul.caution .icon { font-weight: 700; color: var(--ink); }

ol.plan { margin: .5rem 0 0; padding-left: 1.3rem; }
ol.plan li { margin: .28rem 0; font-size: .87rem; color: var(--ink-2); }

.finding {
  border: 1px solid var(--border); border-radius: 10px;
  padding: .55rem .7rem; margin: .4rem 0; background: var(--plane);
}
.finding .head {
  display: flex; gap: .5rem; align-items: baseline; flex-wrap: wrap;
}
.finding .score {
  font-variant-numeric: tabular-nums; font-weight: 650; font-size: .8rem;
  color: var(--ink);
}
.finding .subj { font-weight: 600; font-size: .89rem; }
.finding p { margin: .3rem 0 0; font-size: .86rem; color: var(--ink-2); }
.finding ul { margin: .35rem 0 0; padding-left: 1.05rem; }
.finding li { font-size: .83rem; color: var(--ink-2); }

.kindhead {
  display: flex; align-items: baseline; gap: .5rem;
  margin: 1rem 0 .1rem;
}
.kindhead .n { color: var(--muted); font-size: .8rem; }

figure { margin: .6rem 0 0; }
figure svg { display: block; max-width: 100%; }

/* Stacked bar: the 2px gap is the separator, never a stroke. */
.stack { display: flex; gap: 2px; height: 30px; }
.stack .seg {
  border-radius: 4px; min-width: 3px;
  display: flex; align-items: center; justify-content: center;
  color: var(--surface); font-size: .78rem; font-weight: 650;
  white-space: nowrap;
}
.chart-min { min-width: 560px; }
.matrix-min { min-width: 320px; }

.pill {
  display: inline-block; border-radius: 999px; padding: .04rem .45rem;
  font-size: .74rem; border: 1px solid var(--border); color: var(--ink-2);
  background: var(--plane);
}
.pill.warn { border-color: var(--warning); color: var(--ink); }
.pill.crit { border-color: var(--critical); color: var(--ink); }

.empty {
  color: var(--muted); font-size: .87rem; font-style: normal;
  border: 1px dashed var(--border); border-radius: 10px;
  padding: .6rem .75rem; margin: .5rem 0 0;
}

footer.page {
  color: var(--muted); font-size: .8rem; text-align: center; margin-top: 1.5rem;
}

@media print {
  body { background: #fff; }
  section.card { break-inside: avoid; }
}
"""


# ---------------------------------------------------------------------------
# small helpers
#
# Everything user-supplied goes through _e().  There is no second escaping
# path in this module: if a value came from a file, it came through here.


def _e(value: Any) -> str:
    """Escape anything for HTML/SVG text, including None."""
    if value is None:
        return ""
    return html.escape(str(value), quote=True)


def _num(value: Optional[float], places: int = 0, dash: str = "--") -> str:
    if value is None:
        return dash
    try:
        return f"{float(value):,.{places}f}"
    except (TypeError, ValueError):
        return dash


def _pct(value: Optional[float], places: int = 1, dash: str = "--") -> str:
    if value is None:
        return dash
    try:
        return f"{float(value) * 100:.{places}f}%"
    except (TypeError, ValueError):
        return dash


def _side_key(side: Optional[str]) -> str:
    if side in ("maternal", "paternal", "both"):
        return side
    return "unassigned"


def _side_var(side: Optional[str]) -> str:
    return f"var(--side-{_side_key(side)})"


def _cluster_var(position: int) -> str:
    """Colour for the nth cluster, 1-based.

    Eight categorical slots, then a single shared muted grey.  Hues are never
    cycled or generated: a ninth cluster that reused slot 1's blue would be
    indistinguishable from cluster 1 under colour-vision deficiency, so past
    eight the colour stops carrying identity and the printed index carries it
    instead.  Cluster identity is always also in the row label and the table,
    so nothing here is colour-alone.
    """
    if 1 <= position <= 8:
        return f"var(--c{position})"
    return "var(--c-other)"


def _row_get(row: Any, key: str, default: Any = None) -> Any:
    """Read a column from an sqlite3.Row without exploding if it is absent."""
    try:
        value = row[key]
    except (IndexError, KeyError, TypeError):
        return default
    return default if value is None else value


def _empty(message: str) -> str:
    return f'<p class="empty">{_e(message)}</p>'


def _card(title: str, body: str, lede: str = "") -> str:
    parts = [f"<section class=\"card\"><h2>{_e(title)}</h2>"]
    if lede:
        parts.append(f'<p class="lede">{_e(lede)}</p>')
    parts.append(body)
    parts.append("</section>")
    return "".join(parts)


def _scroll(inner: str) -> str:
    return f'<div class="scroll">{inner}</div>'


def _table(
    headers: Sequence[Tuple[str, bool]],
    rows: Sequence[Sequence[str]],
    caption: str = "",
    wrap_cols: Sequence[int] = (),
) -> str:
    """Build a table.  ``headers`` is (label, numeric?); cells are pre-escaped.

    ``wrap_cols`` names the column indices allowed to wrap onto several lines;
    every other column stays on one line and the table scrolls inside its own
    container rather than pushing the page sideways.
    """
    wrapping = set(wrap_cols)
    head = "".join(
        f'<th class="num">{_e(label)}</th>' if numeric else f"<th>{_e(label)}</th>"
        for label, numeric in headers
    )
    body = []
    for row in rows:
        cells = []
        for i, ((_label, numeric), cell) in enumerate(zip(headers, row)):
            if numeric:
                cells.append(f'<td class="num">{cell}</td>')
            elif i in wrapping:
                cells.append(f'<td class="wrap-cell">{cell}</td>')
            else:
                cells.append(f"<td>{cell}</td>")
        body.append(f"<tr>{''.join(cells)}</tr>")
    cap = f"<caption>{_e(caption)}</caption>" if caption else ""
    return (
        f"<table>{cap}<thead><tr>{head}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table>"
    )


def _legend(entries: Sequence[Tuple[str, str]]) -> str:
    """A legend is always present when two or more colours carry meaning."""
    items = "".join(
        f'<span><i class="dot" style="background:{colour}"></i>{_e(label)}</span>'
        for label, colour in entries
    )
    return f'<div class="legend">{items}</div>'


# ---------------------------------------------------------------------------
# 1. header


def _header_html(
    store: Store,
    kit: Any,
    gmap: GeneticMap,
    title: Optional[str],
    generated: datetime,
) -> str:
    project = store.get_meta("project_name") or PROJECT_NAME
    heading = title or f"{project}: DNA research report"
    person = getattr(kit, "person", None) or "unnamed person"
    label = getattr(kit, "label", None) or "no kit"

    approx = ""
    if gmap.is_approximate:
        approx = (
            " Segment sizes are quoted from a uniform-rate approximation, which "
            "is fine for sanity checks but not accurate enough to quote to "
            "another genealogist -- load a real recombination map before "
            "relying on a cM figure here."
        )

    rows = [
        ("Kit", f"{_e(label)} &mdash; {_e(person)}"),
        ("Genetic map", _e(gmap.describe())),
        ("Generated", _e(generated.strftime("%Y-%m-%d %H:%M UTC"))),
        ("Project file", _e(os.path.basename(store.path))),
    ]
    created = store.get_meta("created_at")
    if created:
        rows.append(("Project started", _e(created)))
    dl = "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in rows)

    return (
        '<header class="page">'
        f"<h1>{_e(heading)}</h1>"
        f'<p class="sub">Everything this project currently knows about '
        f"{_e(person)}&rsquo;s DNA, in one offline file.</p>"
        f"<dl>{dl}</dl>"
        '<div class="privacy"><strong>This file is as sensitive as the raw '
        "data it came from.</strong> It contains identifiable genetic "
        "information about the subject and about living relatives who did not "
        "ask to appear in it: names, how much DNA each shares, which side of "
        "the family they fall on, and which stretches of chromosome came from "
        "whom. Genetic information about one person is also information about "
        "their parents, siblings and children. Do not post it, mail it, or "
        "upload it to a site or an assistant; share it deliberately, with one "
        "named person at a time, and say what it contains when you do."
        f"{_e(approx)}</div>"
        "</header>"
    )


# ---------------------------------------------------------------------------
# 2. kits and QC


def _qc_html(store: Store) -> str:
    kits = store.kits()
    if not kits:
        return _card(
            "Kits and quality control",
            _empty("No kits in this project yet. Import a raw data file to start."),
        )

    rows: List[List[str]] = []
    cautions: List[Tuple[str, bool]] = []
    for kit in kits:
        try:
            q = qc_kit(store, kit)
        except Exception as exc:  # a broken kit must not lose the whole report
            rows.append([
                _e(kit.label), _e(kit.vendor), _e(kit.build),
                "--", "--", "--", "--",
            ])
            cautions.append((f"{kit.label}: QC could not be computed ({exc})", True))
            continue

        if q is None or q.total == 0:
            snps = _num(kit.snp_count) if kit.snp_count else "none"
            rows.append([
                _e(kit.label), _e(kit.vendor or "--"), _e(kit.build or "--"),
                snps, "--", "--", _e(kit.inferred_sex or "--"),
            ])
            cautions.append(
                (
                    f"{kit.label}: no genotypes stored, so call rate, "
                    "heterozygosity and sex inference are unavailable. This is "
                    "normal for a kit whose match list was imported without its "
                    "raw data.",
                    False,
                )
            )
            continue

        rows.append([
            _e(kit.label),
            _e(kit.vendor or "--"),
            _e(kit.build or "--"),
            _num(q.total),
            _pct(q.call_rate),
            _pct(q.het_rate),
            _e(q.inferred_sex or "undetermined"),
        ])
        for warning in q.warnings():
            cautions.append((f"{kit.label}: {warning}", True))
        if q.inferred_sex and kit.inferred_sex and q.inferred_sex != kit.inferred_sex:
            cautions.append(
                (
                    f"{kit.label}: sex inferred from the data now ({q.inferred_sex}) "
                    f"disagrees with the value stored at import "
                    f"({kit.inferred_sex}).",
                    True,
                )
            )

    table = _scroll(
        _table(
            [
                ("Kit", False), ("Vendor", False), ("Build", False),
                ("SNPs", True), ("Call rate", True), ("Heterozygosity", True),
                ("Inferred sex", False),
            ],
            rows,
            caption=(
                "Call rate is the share of tested positions with a genotype; "
                "heterozygosity is the share of called positions with two "
                "different alleles."
            ),
        )
    )

    if cautions:
        items = "".join(
            f'<li class="{"crit" if crit else ""}">'
            f'<span class="icon">{"!" if crit else "i"}</span>'
            f"<span>{_e(text)}</span></li>"
            for text, crit in cautions
        )
        caution_block = f"<h3>Cautions</h3><ul class=\"caution\">{items}</ul>"
    else:
        caution_block = (
            "<h3>Cautions</h3>"
            + _empty("No quality problems found in any kit.")
        )

    return _card("Kits and quality control", table + caution_block)


# ---------------------------------------------------------------------------
# 3. sides


def _sides_html(store: Store, kit_id: Optional[int]) -> str:
    summary: Dict[str, Dict[str, float]] = {}
    if kit_id is not None:
        try:
            summary = side_summary(store, kit_id)
        except Exception:
            summary = {}

    if not summary:
        return _card(
            "Sides",
            _empty(
                "No matches on this kit, so there is nothing to sort onto a "
                "side. Testing a parent is the single most effective way to "
                "split a match list in two."
            ),
            lede=(
                "Which parent each match is related through -- the division "
                "that makes every later step easier."
            ),
        )

    ordered = [k for k in SIDE_ORDER if k in summary]
    ordered += [k for k in sorted(summary) if k not in SIDE_ORDER]

    total_matches = sum(int(v.get("matches", 0)) for v in summary.values())
    total_cm = sum(float(v.get("total_cm", 0.0)) for v in summary.values())

    tiles = []
    for key in ordered:
        bucket = summary[key]
        label = SIDE_LABEL.get(key, key.replace("_", " ").title())
        share = (
            float(bucket.get("total_cm", 0.0)) / total_cm if total_cm else 0.0
        )
        tiles.append(
            '<div class="stat">'
            f'<div class="k"><i class="dot" style="background:{_side_var(key)}">'
            f"</i>{_e(label)}</div>"
            f'<div class="v">{_num(bucket.get("matches"))}</div>'
            f'<div class="n">{_num(bucket.get("total_cm"))} cM total '
            f"&middot; {_pct(share, 0)} of shared cM</div>"
            "</div>"
        )

    bar = _sides_bar(summary, ordered, total_cm)

    rows = [
        [
            f'<i class="swatch" style="background:{_side_var(k)}"></i>'
            + _e(SIDE_LABEL.get(k, k)),
            _num(summary[k].get("matches")),
            _num(summary[k].get("total_cm")),
            _num(summary[k].get("largest_cm")),
        ]
        for k in ordered
    ]
    rows.append([
        "<strong>All matches</strong>",
        f"<strong>{_num(total_matches)}</strong>",
        f"<strong>{_num(total_cm)}</strong>",
        "",
    ])

    table = _scroll(
        _table(
            [("Side", False), ("Matches", True), ("Total cM", True),
             ("Largest match, cM", True)],
            rows,
            caption=(
                "Total cM sums each match's whole-genome sharing, so it "
                "double-counts nothing but says nothing about coverage either "
                "-- for that see the chromosome map."
            ),
        )
    )

    return _card(
        "Sides",
        f'<div class="stats">{"".join(tiles)}</div>' + bar + table,
        lede=(
            "Which parent each match is related through -- the division that "
            "makes every later step easier."
        ),
    )


def _sides_bar(
    summary: Dict[str, Dict[str, float]], ordered: Sequence[str], total_cm: float
) -> str:
    """One stacked bar of shared cM by side.

    Laid out in flexbox rather than SVG: the bar has to survive a phone-width
    column without either stretching its text or forcing the page sideways,
    and a stacked bar is one of the few charts that is genuinely simpler as
    boxes.  A 2px gap in the surface colour separates the segments -- no
    strokes are drawn around them.  A segment narrower than about a tenth of
    the bar carries no inline label, because a clipped label is worse than
    none; the legend and the table below carry every value.
    """
    if total_cm <= 0:
        return ""
    parts: List[str] = []
    for key in ordered:
        value = float(summary[key].get("total_cm", 0.0))
        if value <= 0:
            continue
        share = value / total_cm
        label = SIDE_LABEL.get(key, key)
        inline = _pct(share, 0) if share > 0.15 else ""
        parts.append(
            f'<div class="seg" style="flex:{share:.6f} 1 0;'
            f'background:{_side_var(key)}" '
            f'title="{_e(label)}: {_e(_num(value))} cM '
            f'({_e(_pct(share, 0))})">{_e(inline)}</div>'
        )
    legend = _legend(
        [(SIDE_LABEL.get(k, k), _side_var(k)) for k in ordered
         if float(summary[k].get("total_cm", 0.0)) > 0]
    )
    return (
        '<figure><div class="stack" role="img" '
        'aria-label="Shared centimorgans by side">'
        f"{''.join(parts)}</div>{legend}</figure>"
    )


# ---------------------------------------------------------------------------
# 4. chromosome map


def _chromosome_html(
    store: Store, kit_id: Optional[int], gmap: GeneticMap
) -> str:
    painted: Dict[str, List[PaintedSegment]] = {}
    if kit_id is not None:
        try:
            painted = paint(store, kit_id, gmap)
        except Exception:
            painted = {}

    lede = (
        "Where each match's DNA sits on your chromosomes -- and, more usefully, "
        "where none of them do."
    )
    if not painted:
        return _card(
            "Chromosome map",
            _empty(
                "No segment data. Vendor match lists carry segments at "
                "23andMe, MyHeritage, FTDNA and GEDmatch, but Ancestry does "
                "not publish them at all -- with an Ancestry-only project the "
                "clusters below are the substitute."
            ),
            lede=lede,
        )

    svg, overflow, drawn = _chromosome_svg(painted, gmap)
    try:
        cov = coverage(painted, gmap)
    except Exception:
        cov = {}

    caption_bits = [
        f"{drawn:,} segments of at least 7 cM across {len(painted)} "
        "chromosomes.",
        "Segments are packed into lanes within each chromosome so overlapping "
        f"matches sit above one another rather than hiding each other; up to "
        f"{MAX_LANES} lanes are drawn per chromosome.",
    ]
    if overflow:
        caption_bits.append(
            f"{overflow:,} of them needed a ninth lane or beyond and are "
            "overlaid faintly on the last one, so a dense region reads as "
            "dense rather than losing segments."
        )
    caption = " ".join(caption_bits)

    legend = _legend(
        [(SIDE_LABEL[k], f"var(--side-{k})") for k in SIDE_ORDER]
    )

    body = (
        "<figure>"
        + _scroll(f'<div class="chart-min">{svg}</div>')
        + legend
        + f'<figcaption class="small muted">{_e(caption)}</figcaption>'
        + "</figure>"
    )
    return _card("Chromosome map", body + _coverage_html(cov), lede=lede)


def _chromosome_svg(
    painted: Dict[str, List[PaintedSegment]], gmap: GeneticMap
) -> Tuple[str, int, int]:
    """Draw the 22 autosomes as tracks scaled by physical length.

    Returns ``(svg, overflow_count, drawn_count)``.
    """
    bp = CHROM_BP.get(gmap.build) or CHROM_BP["37"]
    longest = max(bp[c] for c in AUTOSOMES)

    gutter = 34.0        # room for the chromosome number
    right = 10.0
    width = 1000.0
    plot = width - gutter - right
    axis_h = 26.0
    lane_h = 6.0
    lane_gap = 2.0       # the surface gap; nothing is stroked
    track_gap = 9.0

    def x_of(pos: int) -> float:
        return gutter + plot * max(0, min(pos, longest)) / longest

    # -- axis: megabases across the top, with recessive solid hairlines that
    # carry down through every track.
    ticks: List[str] = []
    grid: List[str] = []
    step = 25_000_000
    tick = 0
    while tick <= longest:
        x = x_of(tick)
        ticks.append(
            f'<text x="{x:.1f}" y="14" font-size="11" fill="var(--muted)" '
            f'text-anchor="{"start" if tick == 0 else "middle"}">'
            f"{tick // 1_000_000}</text>"
        )
        if tick:
            grid.append(x)
        tick += step
    ticks.append(
        f'<text x="{gutter:.1f}" y="{axis_h - 1:.1f}" font-size="10" '
        f'fill="var(--muted)">Mb</text>'
    )

    # -- tracks
    tracks: List[str] = []
    y = axis_h + 6.0
    overflow_total = 0
    drawn_total = 0
    grid_bottom = y

    for chrom in AUTOSOMES:
        segs = sorted(
            painted.get(chrom, []), key=lambda s: (s.start_bp, s.end_bp)
        )
        lanes, overflow = _pack_lanes(segs, MAX_LANES)
        overflow_total += len(overflow)
        drawn_total += len(segs)
        n_lanes = max(len(lanes), 1)
        track_h = n_lanes * lane_h + (n_lanes - 1) * lane_gap

        rail_y = y + (track_h - lane_h) / 2.0
        tracks.append(
            f'<text x="{gutter - 6:.1f}" y="{y + track_h / 2 + 4:.1f}" '
            f'font-size="11" fill="var(--ink-2)" text-anchor="end">'
            f"{_e(chrom)}</text>"
        )
        tracks.append(
            f'<rect x="{gutter:.1f}" y="{rail_y:.1f}" '
            f'width="{x_of(bp[chrom]) - gutter:.1f}" height="{lane_h:.1f}" '
            f'rx="3" fill="var(--grid)"></rect>'
        )

        for lane_i, lane in enumerate(lanes):
            lane_y = y + lane_i * (lane_h + lane_gap)
            for seg in lane:
                tracks.append(_segment_rect(seg, x_of, lane_y, lane_h, 1.0))
        if overflow:
            lane_y = y + (n_lanes - 1) * (lane_h + lane_gap)
            for seg in overflow:
                tracks.append(_segment_rect(seg, x_of, lane_y, lane_h, 0.35))

        y += track_h + track_gap
        grid_bottom = y

    height = y + 4.0
    gridlines = "".join(
        f'<line x1="{x:.1f}" y1="{axis_h:.1f}" x2="{x:.1f}" '
        f'y2="{grid_bottom - track_gap:.1f}" stroke="var(--grid)" '
        f'stroke-width="1"></line>'
        for x in grid
    )

    svg = (
        f'<svg viewBox="0 0 {width:.0f} {height:.0f}" '
        'style="width:100%;height:auto" role="img" '
        'aria-label="Chromosome map of matching segments by side">'
        f"{gridlines}{''.join(ticks)}{''.join(tracks)}</svg>"
    )
    return svg, overflow_total, drawn_total


def _pack_lanes(
    segs: Sequence[PaintedSegment], max_lanes: int
) -> Tuple[List[List[PaintedSegment]], List[PaintedSegment]]:
    """Greedy interval packing: first lane whose last segment has ended.

    Segments arrive sorted by start position, so the first lane that is free
    is also the one that keeps the picture compact.  Anything that would need
    a lane past the cap is returned separately and drawn faintly over the last
    lane -- visible as density, counted in the caption, never dropped.
    """
    lanes: List[List[PaintedSegment]] = []
    lane_end: List[int] = []
    overflow: List[PaintedSegment] = []
    # 1 Mb of padding keeps two nearly-touching segments visually separate.
    pad = 1_000_000
    for seg in segs:
        placed = False
        for i, end in enumerate(lane_end):
            if seg.start_bp > end + pad:
                lanes[i].append(seg)
                lane_end[i] = seg.end_bp
                placed = True
                break
        if placed:
            continue
        if len(lanes) < max_lanes:
            lanes.append([seg])
            lane_end.append(seg.end_bp)
        else:
            overflow.append(seg)
    return lanes, overflow


def _segment_rect(
    seg: PaintedSegment, x_of, y: float, height: float, opacity: float
) -> str:
    x0 = x_of(seg.start_bp)
    x1 = x_of(seg.end_bp)
    w = max(x1 - x0, 1.2)
    side = SIDE_LABEL[_side_key(seg.side)]
    tip = (
        f"{seg.name} -- chr{seg.chrom}:"
        f"{seg.start_bp / 1e6:.1f}-{seg.end_bp / 1e6:.1f} Mb, "
        f"{seg.cm:.1f} cM, {side.lower()}"
    )
    op = "" if opacity >= 1.0 else f' opacity="{opacity}"'
    return (
        f"<g{op}><title>{_e(tip)}</title>"
        f'<rect x="{x0:.1f}" y="{y:.1f}" width="{w:.1f}" '
        f'height="{height:.1f}" rx="2" fill="{_side_var(seg.side)}"></rect></g>'
    )


def _coverage_html(cov: Dict[str, float]) -> str:
    if not cov:
        return ""
    genome = float(cov.get("genome") or 0.0)
    any_cm = float(cov.get("any") or 0.0)
    unattributed = float(cov.get("unattributed") or 0.0)

    def share(value: float) -> str:
        return _pct(value / genome, 0) if genome else "--"

    tiles = [
        ("Maternal", cov.get("maternal"), "var(--side-maternal)"),
        ("Paternal", cov.get("paternal"), "var(--side-paternal)"),
        ("Covered by anyone", any_cm, "var(--side-both)"),
        ("Explained by nobody", unattributed, "var(--side-unassigned)"),
    ]
    cards = "".join(
        '<div class="stat">'
        f'<div class="k"><i class="dot" style="background:{colour}"></i>'
        f"{_e(label)}</div>"
        f'<div class="v">{_num(value)}<span class="n"> cM</span></div>'
        f'<div class="n">{share(float(value or 0.0))} of {_num(genome)} cM</div>'
        "</div>"
        for label, value, colour in tiles
    )
    return (
        "<h3>Coverage</h3>"
        f'<div class="stats">{cards}</div>'
        '<p class="small muted">Coverage counts a stretch once however many '
        "matches cover it, so it does not add up to the total shared cM above. "
        "The unattributed figure is the genealogically interesting one: it is "
        "the part of the genome no current match explains, and it is where the "
        "undiscovered branches are.</p>"
    )


# ---------------------------------------------------------------------------
# 5. clusters


def _clusters_html(
    store: Store, kit_id: Optional[int], result: Optional[ClusterResult]
) -> str:
    lede = (
        "The Leeds picture: matches ordered by cluster, with a cell wherever "
        "two of them match each other. Shared ancestry shows up as a block on "
        "the diagonal."
    )
    if result is None:
        return _card(
            "Clusters",
            _empty(
                "No clustering run was supplied. Clustering needs an "
                "in-common-with export rather than just the match list, and it "
                "is the main analytical tool when segment data is unavailable."
            ),
            lede=lede,
        )
    if not result.clusters:
        notes = "".join(f"<li>{_e(w)}</li>" for w in result.warnings)
        body = _empty("The clustering run produced no clusters.")
        if notes:
            body += f'<ul class="notes">{notes}</ul>'
        return _card("Clusters", body, lede=lede)

    ids = _edge_index_ids(store, kit_id, result)
    pairs = set()
    for edge in result.edges:
        try:
            i, j = edge
        except (TypeError, ValueError):
            continue
        if 0 <= i < len(ids) and 0 <= j < len(ids) and i != j:
            a, b = ids[i], ids[j]
            pairs.add((min(a, b), max(a, b)))

    display, cluster_of, truncated, shown_clusters = _matrix_order(result)
    matrix = _matrix_svg(display, cluster_of, pairs, result)

    caption = (
        f"{len(display)} matches shown"
        + (
            f", from the {shown_clusters} largest of "
            f"{len(result.clusters)} clusters -- the full matrix would exceed "
            f"{MAX_MATRIX}x{MAX_MATRIX} cells."
            if truncated
            else f" across {len(result.clusters)} clusters."
        )
    )

    rows = []
    for c in result.clusters:
        position = result.clusters.index(c) + 1
        side = c.side or "unassigned"
        side_cell = _e(side)
        if side == "mixed":
            side_cell = f'<span class="pill warn">{_e(side)}</span>'
        rows.append([
            f'<i class="swatch" style="background:{_cluster_var(position)}">'
            f"</i>cluster {_e(c.index)}",
            _num(c.size),
            side_cell,
            _pct(c.cohesion, 0),
            _num(c.top_cm),
            _num(c.total_cm),
            _e(", ".join(c.names[:3]) + ("..." if len(c.names) > 3 else "")),
        ])

    table = _scroll(
        _table(
            [
                ("Cluster", False), ("Size", True), ("Side", False),
                ("Cohesion", True), ("Strongest, cM", True),
                ("Total cM", True), ("Members", False),
            ],
            rows,
            caption=(
                "Cohesion is the share of possible within-cluster shared-match "
                "links that are actually present. A cluster near 100% descends "
                "from one couple; a loose one has probably swallowed two "
                "branches and should be split by lowering the distance cutoff."
            ),
            wrap_cols=(6,),
        )
    )

    warnings_html = ""
    if result.warnings:
        items = "".join(
            f'<li><span class="icon">i</span><span>{_e(w)}</span></li>'
            for w in result.warnings
        )
        warnings_html = f'<ul class="caution">{items}</ul>'

    band = (
        f'<p class="small muted">Band {result.band[0]:.0f}-{result.band[1]:.0f} '
        f"cM. {_e(result.summary())}</p>"
    )

    body = (
        band
        + "<figure>"
        + _scroll(f'<div class="matrix-min">{matrix}</div>')
        + f'<figcaption class="small muted">{_e(caption)} Cells inside a '
        "cluster take that cluster&rsquo;s colour; a cell between two clusters "
        "is drawn in grey, and a lot of grey means the boundary between them "
        "may not be real. Colour repeats past the eighth cluster, so the "
        "printed index is what identifies a cluster, not its hue."
        "</figcaption></figure>"
        + table
        + warnings_html
    )
    return _card("Clusters", body, lede=lede)


def _edge_index_ids(
    store: Store, kit_id: Optional[int], result: ClusterResult
) -> List[int]:
    """Work out which match ids the edge indices refer to.

    ``ClusterResult.edges`` holds index pairs, but the index space is the
    clusterer's input order (matches sorted by shared cM, descending) rather
    than the cluster-grouped ``order``.  Rather than trust either reading, we
    reconstruct the input order from the store and score both against the
    shared-match table, keeping whichever actually reproduces real pairs.  If
    there is nothing to score against, we take the reconstruction, which is
    what the current clusterer produces.
    """
    order = list(result.order)
    if not order:
        return []
    in_order = set(order)

    reconstructed = order
    if kit_id is not None:
        try:
            by_cm = [m["id"] for m in store.matches(kit_id) if m["id"] in in_order]
            if len(by_cm) == len(order):
                reconstructed = by_cm
        except Exception:
            pass

    truth = set()
    if kit_id is not None:
        try:
            for row in store.shared_matches(kit_id):
                a, b = row["a_id"], row["b_id"]
                truth.add((min(a, b), max(a, b)))
        except Exception:
            truth = set()
    if not truth:
        return reconstructed

    best, best_score = reconstructed, -1
    for candidate in (reconstructed, order):
        score = 0
        for edge in result.edges:
            try:
                i, j = edge
            except (TypeError, ValueError):
                continue
            if 0 <= i < len(candidate) and 0 <= j < len(candidate):
                a, b = candidate[i], candidate[j]
                if (min(a, b), max(a, b)) in truth:
                    score += 1
        if score > best_score:
            best, best_score = candidate, score
    return best


def _matrix_order(
    result: ClusterResult,
) -> Tuple[List[int], Dict[int, int], bool, int]:
    """Pick the matches to draw, largest clusters first, under the cell cap."""
    display: List[int] = []
    cluster_of: Dict[int, int] = {}
    positions = sorted(
        range(len(result.clusters)),
        key=lambda i: -result.clusters[i].size,
    )
    shown = 0
    truncated = False
    for i in positions:
        c = result.clusters[i]
        if len(display) + c.size > MAX_MATRIX and display:
            truncated = True
            continue
        for mid in c.match_ids:
            cluster_of[mid] = i + 1
        display.extend(c.match_ids[:MAX_MATRIX])
        shown += 1
    if len(display) > MAX_MATRIX:
        display = display[:MAX_MATRIX]
        truncated = True
    # Restore cluster-index order so the blocks run down the diagonal in the
    # same sequence the table lists them.
    display.sort(key=lambda mid: (cluster_of.get(mid, 10 ** 6), 0))
    remaining = MAX_MATRIX - len(display)
    if remaining > 0 and result.unclustered:
        extra = result.unclustered[:remaining]
        display.extend(extra)
        if len(extra) < len(result.unclustered):
            truncated = True
    elif result.unclustered:
        truncated = True
    return display, cluster_of, truncated, shown


def _matrix_svg(
    display: Sequence[int],
    cluster_of: Dict[int, int],
    pairs: set,
    result: ClusterResult,
) -> str:
    n = len(display)
    if n == 0:
        return ""
    cell = 8.0
    gutter = 58.0
    pad = 6.0
    size = n * cell
    width = gutter + size + pad
    height = size + pad * 2

    parts: List[str] = []

    # Cluster blocks: a wash behind each cluster's rows and columns, so the
    # eye finds the block before it reads any individual cell.
    start = 0
    while start < n:
        cid = cluster_of.get(display[start])
        end = start
        while end + 1 < n and cluster_of.get(display[end + 1]) == cid:
            end += 1
        if cid is not None:
            x = gutter + start * cell
            span = (end - start + 1) * cell
            parts.append(
                f'<rect x="{x:.1f}" y="{pad + start * cell:.1f}" '
                f'width="{span:.1f}" height="{span:.1f}" rx="2" '
                f'fill="{_cluster_var(cid)}" opacity="0.14"></rect>'
            )
            label = result.clusters[cid - 1].index if cid <= len(result.clusters) else cid
            parts.append(
                f'<text x="{gutter - 6:.1f}" '
                f'y="{pad + start * cell + span / 2 + 3.5:.1f}" font-size="10" '
                f'fill="var(--ink-2)" text-anchor="end">cluster '
                f"{_e(label)}</text>"
            )
        start = end + 1

    index_of = {mid: i for i, mid in enumerate(display)}
    for i, a in enumerate(display):
        # The diagonal: a match always shares itself, drawn recessively so it
        # reads as the frame of the block rather than as evidence.
        parts.append(
            f'<rect x="{gutter + i * cell + 1:.1f}" y="{pad + i * cell + 1:.1f}" '
            f'width="{cell - 2:.1f}" height="{cell - 2:.1f}" rx="1.5" '
            f'fill="var(--axis)" opacity="0.5"></rect>'
        )
        for b in display[i + 1:]:
            if (min(a, b), max(a, b)) not in pairs:
                continue
            j = index_of[b]
            same = (
                cluster_of.get(a) is not None
                and cluster_of.get(a) == cluster_of.get(b)
            )
            fill = _cluster_var(cluster_of[a]) if same else "var(--muted)"
            opacity = "1" if same else "0.5"
            name_a = _e(result.labels.get(a, a))
            name_b = _e(result.labels.get(b, b))
            tip = f"<title>{name_a} &amp; {name_b}</title>"
            for (r, c) in ((i, j), (j, i)):
                parts.append(
                    f"<g>{tip}<rect "
                    f'x="{gutter + c * cell + 1:.1f}" '
                    f'y="{pad + r * cell + 1:.1f}" '
                    f'width="{cell - 2:.1f}" height="{cell - 2:.1f}" rx="1.5" '
                    f'fill="{fill}" opacity="{opacity}"></rect></g>'
                )

    return (
        f'<svg viewBox="0 0 {width:.0f} {height:.0f}" '
        f'width="{width:.0f}" height="{height:.0f}" '
        'style="max-width:100%;height:auto" role="img" '
        'aria-label="Cluster matrix of matches that match each other">'
        f"{''.join(parts)}</svg>"
    )


# ---------------------------------------------------------------------------
# 6. top matches


def _matches_html(
    store: Store,
    kit_id: Optional[int],
    gmap: GeneticMap,
    predict_top: int,
    iterations: int,
) -> str:
    lede = (
        "The strongest matches, with the vendor's guess beside this project's "
        "own -- which is a ranked set of possibilities, because adjacent "
        "relationships overlap heavily in shared cM."
    )
    rows_in: List[Any] = []
    total = 0
    if kit_id is not None:
        try:
            rows_in = store.matches(kit_id)
            total = len(rows_in)
        except Exception:
            rows_in = []
    if not rows_in:
        return _card(
            "Top matches",
            _empty("No matches imported for this kit."),
            lede=lede,
        )

    rows_in = sorted(rows_in, key=lambda m: -(m["total_cm"] or 0.0))
    shown = rows_in[:max(predict_top, 0)]

    memo: Dict[Tuple[int, Optional[int]], str] = {}
    failed = False
    out_rows: List[List[str]] = []
    for m in shown:
        cm = float(m["total_cm"] or 0.0)
        segs = _row_get(m, "seg_count")
        key = (int(round(cm)), int(segs) if segs else None)
        if key in memo:
            prediction = memo[key]
        else:
            prediction, ok = _predict_cell(store, cm, segs, gmap, iterations)
            failed = failed or not ok
            memo[key] = prediction

        side = _side_key(_row_get(m, "side"))
        out_rows.append([
            _e(m["name"] or m["remote_id"] or f"match {m['id']}"),
            _num(cm),
            _num(_row_get(m, "seg_count")),
            f'<i class="swatch" style="background:{_side_var(side)}"></i>'
            + _e(SIDE_LABEL[side]),
            _e(_row_get(m, "predicted") or "--"),
            prediction,
        ])

    caption = (
        f"Showing the {len(shown)} strongest of {total:,} matches; each row "
        f"costs a simulation lookup, so the table is capped at {predict_top}."
    )
    if failed:
        caption += (
            " Some predictions could not be computed and are shown as "
            "unavailable."
        )

    table = _scroll(
        _table(
            [
                ("Match", False), ("Shared cM", True), ("Segments", True),
                ("Side", False), ("Vendor's guess", False),
                ("This project's top three", False),
            ],
            out_rows,
            caption=caption,
            wrap_cols=(5,),
        )
    )
    return _card("Top matches", table, lede=lede)


def _predict_cell(
    store: Store,
    cm: float,
    segments: Optional[int],
    gmap: GeneticMap,
    iterations: int,
) -> Tuple[str, bool]:
    if cm <= 0:
        return '<span class="muted">no sharing recorded</span>', True
    try:
        pred = predict(
            store, cm, gmap=gmap, iterations=iterations,
            segments=int(segments) if segments else None,
        )
        top = pred.top(3)
    except Exception:
        return '<span class="muted">prediction unavailable</span>', False
    if not top:
        return '<span class="muted">no candidate relationship fits</span>', True
    bits = []
    for c in top:
        rng = c.range_text()
        extra = f' <span class="muted">{_e(rng)}</span>' if rng else ""
        bits.append(
            f'<div><span class="mono">{_pct(c.posterior, 0)}</span> '
            f"{_e(c.name)}{extra}</div>"
        )
    return "".join(bits), True


# ---------------------------------------------------------------------------
# 7. triangulation groups


def _triangulation_html(
    store: Store, kit_id: Optional[int], gmap: GeneticMap, limit: int = 25
) -> str:
    lede = (
        "Regions where several matches share the same stretch with you. "
        "Overlapping is not triangulating: you have two copies of every "
        "chromosome, so matches on opposite copies can overlap while being "
        "entirely unrelated to each other."
    )
    groups: List[OverlapGroup] = []
    if kit_id is not None:
        try:
            groups = overlap_groups(store, kit_id, gmap)
        except Exception:
            groups = []
    if not groups:
        return _card(
            "Triangulation groups",
            _empty(
                "No region has two or more matches sharing at least 7 cM with "
                "you. This needs segment data and at least a couple of matches "
                "in the same place."
            ),
            lede=lede,
        )

    conflicts = sum(1 for g in groups if g.conflict)
    rows = []
    for g in groups[:limit]:
        side = g.side()
        if side == "conflicting":
            side_cell = '<span class="pill crit">conflicting sides</span>'
        elif side:
            side_cell = (
                f'<i class="swatch" style="background:{_side_var(side)}"></i>'
                + _e(SIDE_LABEL.get(side, side))
            )
        else:
            side_cell = f'<span class="muted">{_e(SIDE_LABEL["unassigned"])}</span>'
        support = _e(g.support)
        if g.support == "triangulated":
            support = f'<span class="pill">{support}</span>'
        rows.append([
            _e(f"chr{g.chrom}:{g.start_bp / 1e6:.1f}-{g.end_bp / 1e6:.1f} Mb"),
            _num(g.cm, 1),
            _num(g.size),
            support,
            _pct(g.icw_fraction, 0),
            side_cell,
            _e(", ".join(g.names[:4]) + ("..." if len(g.names) > 4 else "")),
        ])

    table = _scroll(
        _table(
            [
                ("Location", False), ("cM", True), ("Matches", True),
                ("Support", False), ("Pairs in common", True),
                ("Side", False), ("Members", False),
            ],
            rows,
            caption=(
                f"Top {len(rows)} of {len(groups)} groups, largest first. "
                "Support says which standard the group actually meets: "
                "'overlap only' means the members share the region with you "
                "but are not known to match each other, which is the weakest "
                "reading."
            ),
            wrap_cols=(6,),
        )
    )

    flag = ""
    if conflicts:
        flag = (
            '<ul class="caution"><li class="crit"><span class="icon">!</span>'
            f"<span>{conflicts} group(s) contain matches from both sides, so "
            "those segments sit on opposite copies of the chromosome. Treating "
            "such a group as one ancestral line would merge two unrelated "
            "families.</span></li></ul>"
        )
    return _card("Triangulation groups", table + flag, lede=lede)


# ---------------------------------------------------------------------------
# 8. findings


def _findings_html(
    store: Store, findings: Optional[Sequence[Finding]]
) -> str:
    lede = "Where the data disagrees with itself, and what to do about it."
    items = list(findings) if findings else []
    if not items:
        items = _findings_from_store(store)
    if not items:
        return _card(
            "Findings",
            _empty(
                "No findings. Run the hypothesis engine once there are matches "
                "and a tree to read against each other."
            ),
            lede=lede,
        )

    ordered = sorted(items, key=lambda f: -float(f.score or 0.0))
    # Group by kind, but keep the groups themselves in priority order so the
    # most important thing is still the first thing on the page.
    kind_rank: Dict[str, float] = {}
    for f in ordered:
        kind_rank.setdefault(f.kind or "other", float(f.score or 0.0))

    blocks: List[str] = []
    for kind in sorted(kind_rank, key=lambda k: -kind_rank[k]):
        group = [f for f in ordered if (f.kind or "other") == kind]
        label = (kind or "other").replace("-", " ").replace("_", " ")
        blocks.append(
            f'<div class="kindhead"><h3>{_e(label)}</h3>'
            f'<span class="n">{len(group)}</span></div>'
        )
        for f in group:
            actions = ""
            if f.actions:
                actions = "<ul>" + "".join(
                    f"<li>{_e(a)}</li>" for a in f.actions
                ) + "</ul>"
            blocks.append(
                '<div class="finding"><div class="head">'
                f'<span class="score">{float(f.score or 0.0):.2f}</span>'
                f'<span class="subj">{_e(f.subject)}</span></div>'
                f"<p>{_e(f.summary)}</p>{actions}</div>"
            )

    plan = research_plan(ordered)
    if plan:
        plan_html = (
            "<h3>Research plan</h3><ol class=\"plan\">"
            + "".join(f"<li>{_e(step)}</li>" for step in plan)
            + "</ol>"
        )
    else:
        plan_html = "<h3>Research plan</h3>" + _empty(
            "No findings carried a suggested action."
        )

    return _card("Findings", "".join(blocks) + plan_html, lede=lede)


def _findings_from_store(store: Store) -> List[Finding]:
    """Fall back to whatever the hypothesis engine last persisted."""
    out: List[Finding] = []
    try:
        rows = store.hypotheses()
    except Exception:
        return out
    for row in rows:
        detail: Dict[str, Any] = {}
        raw = _row_get(row, "detail")
        if raw:
            try:
                loaded = json.loads(raw)
                if isinstance(loaded, dict):
                    detail = loaded
            except (ValueError, TypeError):
                detail = {}
        out.append(
            Finding(
                kind=_row_get(row, "kind", "other") or "other",
                subject=_row_get(row, "subject", "") or "",
                summary=_row_get(row, "summary", "") or "",
                score=float(_row_get(row, "score", 0.0) or 0.0),
                detail=detail,
            )
        )
    return out


# ---------------------------------------------------------------------------
# 9. tree summary


def _tree_html(store: Store) -> str:
    lede = "What the paper trail currently covers."
    try:
        individuals = store.individuals()
        families = store.families()
    except Exception:
        individuals, families = [], []

    if not individuals and not families:
        return _card(
            "Tree summary",
            _empty(
                "No tree imported. A GEDCOM makes every match with a tree link "
                "testable against simulated sharing, which is where most "
                "findings come from."
            ),
            lede=lede,
        )

    surnames: Dict[str, int] = {}
    with_birth = 0
    for row in individuals:
        name = (_row_get(row, "surname", "") or "").strip()
        if name:
            surnames[name] = surnames.get(name, 0) + 1
        if _row_get(row, "birth_year"):
            with_birth += 1

    tiles = [
        ("Individuals", len(individuals), ""),
        ("Families", len(families), ""),
        ("Distinct surnames", len(surnames), ""),
        ("With a birth year", with_birth,
         _pct(with_birth / len(individuals), 0) if individuals else ""),
    ]
    cards = "".join(
        '<div class="stat">'
        f'<div class="k">{_e(label)}</div>'
        f'<div class="v">{_num(value)}</div>'
        f'<div class="n">{_e(note)}</div></div>'
        for label, value, note in tiles
    )

    top = sorted(surnames.items(), key=lambda kv: (-kv[1], kv[0]))[:15]
    if top:
        biggest = top[0][1]
        rows = []
        for name, count in top:
            # A single series, so slot 1 for every bar and no legend box: the
            # column header already says what is plotted.
            bar_w = 140.0 * count / biggest if biggest else 0.0
            bar = (
                '<span style="display:inline-block;vertical-align:middle;'
                f"width:{max(bar_w, 3.0):.0f}px;height:8px;border-radius:4px;"
                'background:var(--c1)"></span>'
            )
            rows.append([_e(name), _num(count), bar])
        table = _scroll(
            _table(
                [("Surname", False), ("Individuals", True), ("", False)],
                rows,
                caption=(
                    "The 15 commonest surnames in the tree. A surname that "
                    "recurs across an unidentified cluster is usually the "
                    "fastest way into it."
                ),
            )
        )
    else:
        table = _empty("No surnames recorded on any individual.")

    return _card(
        "Tree summary", f'<div class="stats">{cards}</div>' + table, lede=lede
    )


# ---------------------------------------------------------------------------
# assembly


def build_report(
    store: Store,
    kit_id: Any,
    gmap: GeneticMap,
    out_path: str,
    findings: Optional[Sequence[Finding]] = None,
    cluster_result: Optional[ClusterResult] = None,
    title: Optional[str] = None,
    predict_top: int = 15,
    iterations: int = 1500,
) -> str:
    """Write the whole investigation to one self-contained HTML file.

    Returns the path written.  Every section tolerates its inputs being
    missing: a project with one imported kit and nothing else produces a short
    but valid report rather than an exception.
    """
    kit = None
    try:
        kit = store.kit(kit_id)
    except Exception:
        kit = None
    if kit is not None:
        kid: Optional[int] = kit.id
    elif isinstance(kit_id, int):
        kid = kit_id
    elif str(kit_id).isdigit():
        kid = int(kit_id)
    else:
        kid = None

    generated = datetime.now(timezone.utc)
    sections = [
        _header_html(store, kit, gmap, title, generated),
        _qc_html(store),
        _sides_html(store, kid),
        _chromosome_html(store, kid, gmap),
        _clusters_html(store, kid, cluster_result),
        _matches_html(store, kid, gmap, predict_top, iterations),
        _triangulation_html(store, kid, gmap),
        _findings_html(store, findings),
        _tree_html(store),
    ]

    doc_title = title or (
        f"{PROJECT_NAME}: {kit.label} report" if kit else f"{PROJECT_NAME} report"
    )
    footer = (
        '<footer class="page">Generated by '
        f"{_e(PROJECT_NAME)} on "
        f"{_e(generated.strftime('%Y-%m-%d %H:%M UTC'))}. "
        "This file is self-contained and makes no network requests. "
        "Keep it where you would keep the raw data.</footer>"
    )

    html_doc = (
        "<!DOCTYPE html>\n"
        '<html lang="en">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<meta name="robots" content="noindex, noarchive, noimageindex">\n'
        '<meta name="referrer" content="no-referrer">\n'
        '<meta name="color-scheme" content="light dark">\n'
        f"<title>{_e(doc_title)}</title>\n"
        f"<style>{_CSS}</style>\n"
        "</head>\n<body>\n"
        f'<div class="wrap">{"".join(sections)}{footer}</div>\n'
        "</body>\n</html>\n"
    )

    parent = os.path.dirname(os.path.abspath(out_path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html_doc)
    return out_path


# ---------------------------------------------------------------------------
# plain text
#
# The same story at terminal width.  No ANSI: this gets piped into files and
# emails as often as it gets read on screen.


def build_text_summary(
    store: Store,
    kit_id: Any,
    gmap: GeneticMap,
    findings: Optional[Sequence[Finding]] = None,
) -> str:
    """A compact plain-text version of the report, for terminal output."""
    try:
        kit = store.kit(kit_id)
    except Exception:
        kit = None
    if kit is not None:
        kid: Optional[int] = kit.id
    elif isinstance(kit_id, int):
        kid = kit_id
    elif str(kit_id).isdigit():
        kid = int(kit_id)
    else:
        kid = None

    lines: List[str] = []

    def head(text: str) -> None:
        lines.append("")
        lines.append(text)
        lines.append("-" * len(text))

    label = getattr(kit, "label", None) or str(kit_id)
    person = getattr(kit, "person", None) or "unnamed person"
    lines.append(f"{PROJECT_NAME}: report for {label} ({person})")
    lines.append(gmap.describe())
    lines.append(
        datetime.now(timezone.utc).strftime("generated %Y-%m-%d %H:%M UTC")
    )
    lines.append(
        "Contains identifiable genetic information about the subject and "
        "their relatives -- do not share casually."
    )

    # kits
    head("Kits")
    kits = store.kits()
    if not kits:
        lines.append("  no kits imported")
    for k in kits:
        try:
            q = qc_kit(store, k)
        except Exception:
            q = None
        if q is None or q.total == 0:
            lines.append(
                f"  {k.label}: {k.vendor or '?'} build {k.build or '?'}, "
                "no genotypes stored"
            )
            continue
        lines.append(
            f"  {k.label}: {k.vendor or '?'} build {k.build or '?'}, "
            f"{q.total:,} SNPs, call rate {q.call_rate:.1%}, "
            f"het {q.het_rate:.1%}, sex {q.inferred_sex or '?'}"
        )
        for w in q.warnings():
            lines.append(f"    ! {w}")

    # sides
    head("Sides")
    summary: Dict[str, Dict[str, float]] = {}
    if kid is not None:
        try:
            summary = side_summary(store, kid)
        except Exception:
            summary = {}
    if not summary:
        lines.append("  no matches on this kit")
    else:
        keys = [k for k in SIDE_ORDER if k in summary]
        keys += [k for k in sorted(summary) if k not in SIDE_ORDER]
        for key in keys:
            b = summary[key]
            lines.append(
                f"  {SIDE_LABEL.get(key, key):<12} "
                f"{int(b.get('matches', 0)):>5} matches  "
                f"{b.get('total_cm', 0.0):>9,.0f} cM total  "
                f"largest {b.get('largest_cm', 0.0):,.0f} cM"
            )

    # coverage
    head("Chromosome coverage")
    painted: Dict[str, List[PaintedSegment]] = {}
    if kid is not None:
        try:
            painted = paint(store, kid, gmap)
        except Exception:
            painted = {}
    if not painted:
        lines.append("  no segment data")
    else:
        cov = coverage(painted, gmap)
        genome = cov.get("genome", 0.0) or 1.0
        for key in ("maternal", "paternal", "any", "unattributed"):
            value = cov.get(key, 0.0)
            name = {"any": "covered", "unattributed": "unexplained"}.get(key, key)
            lines.append(
                f"  {name:<12} {value:>9,.0f} cM  ({value / genome:>5.1%})"
            )
        lines.append(
            "  the unexplained portion is where the undiscovered branches are"
        )

    # matches
    head("Strongest matches")
    rows: List[Any] = []
    if kid is not None:
        try:
            rows = sorted(
                store.matches(kid), key=lambda m: -(m["total_cm"] or 0.0)
            )
        except Exception:
            rows = []
    if not rows:
        lines.append("  none")
    for m in rows[:10]:
        name = m["name"] or m["remote_id"] or f"match {m['id']}"
        lines.append(
            f"  {str(name)[:32]:<32} {float(m['total_cm'] or 0.0):>8,.0f} cM  "
            f"{_side_key(_row_get(m, 'side')):<10} "
            f"{_row_get(m, 'predicted', '') or ''}"
        )
    if len(rows) > 10:
        lines.append(f"  ... and {len(rows) - 10:,} more")

    # triangulation
    head("Triangulation groups")
    groups: List[OverlapGroup] = []
    if kid is not None:
        try:
            groups = overlap_groups(store, kid, gmap)
        except Exception:
            groups = []
    if not groups:
        lines.append("  none")
    for g in groups[:8]:
        flag = "  [CONFLICTING SIDES]" if g.conflict else ""
        lines.append(f"  {g.describe()}{flag}")

    # findings
    head("Findings")
    items = list(findings) if findings else _findings_from_store(store)
    if not items:
        lines.append("  none")
    for f in sorted(items, key=lambda f: -float(f.score or 0.0))[:12]:
        lines.append(f"  {f.line()}")
    plan = research_plan(sorted(items, key=lambda f: -float(f.score or 0.0)))
    if plan:
        head("Research plan")
        for i, step in enumerate(plan, start=1):
            lines.append(f"  {i}. {step}")

    # tree
    head("Tree")
    try:
        individuals = store.individuals()
        families = store.families()
    except Exception:
        individuals, families = [], []
    if not individuals and not families:
        lines.append("  no tree imported")
    else:
        surnames: Dict[str, int] = {}
        for row in individuals:
            name = (_row_get(row, "surname", "") or "").strip()
            if name:
                surnames[name] = surnames.get(name, 0) + 1
        lines.append(
            f"  {len(individuals):,} individuals, {len(families):,} families, "
            f"{len(surnames):,} distinct surnames"
        )
        top = sorted(surnames.items(), key=lambda kv: (-kv[1], kv[0]))[:15]
        if top:
            lines.append(
                "  commonest: "
                + ", ".join(f"{n} ({c})" for n, c in top)
            )

    lines.append("")
    return "\n".join(lines)

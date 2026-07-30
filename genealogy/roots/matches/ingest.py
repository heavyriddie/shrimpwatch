"""Import DNA match lists, segment data, and shared-match data.

There is no standard format here.  Every testing site exports something
different, most of them export nothing officially and users rely on browser
extensions, and the column names drift between years.  So rather than
hard-coding one layout per vendor, this module recognises columns by meaning
and tolerates whatever order and spelling it finds.

What we import, and why each matters:

*Match lists* -- who matches you and by how much.  The raw material.

*Segment data* -- where on the chromosomes each match sits.  Only FTDNA,
MyHeritage, GEDmatch and 23andMe expose this; Ancestry does not, which is
why Ancestry work leans on clustering instead of chromosome mapping.

*Shared-match ("in common with") data* -- which of your matches also match
each other.  This is what makes clustering possible, and clustering is what
lets you sort an Ancestry match list into ancestral lines with no segment
data at all.

One conversion deserves a warning.  23andMe reports a percentage rather than
centimorgans, and the percentage counts both chromosome copies, so a parent
shows 50% and a full sibling also shows about 50% -- even though their
half-identical totals are 3485 cM and 2613 cM respectively.  Converting a
percentage to cM by simple scaling therefore misstates full siblings and
closer.  We do the conversion when there is no cM column, but we flag it.
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import zipfile
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from ..genome import normalize_chrom
from ..store import Store

#: 23andMe percentages count both chromosome copies, so 100% would be
#: two genomes' worth.  Used only when no centimorgan column exists.
PERCENT_TO_CM = 6970.0

_COLUMN_ALIASES: Dict[str, Tuple[str, ...]] = {
    "name": (
        "match name", "name", "display name", "full name", "matchname",
        "tester name", "relative name", "match", "person",
    ),
    "remote_id": (
        "match id", "id", "profile id", "guid", "kit number", "kit no", "kit",
        "testid", "test id", "human id", "profile",
    ),
    "total_cm": (
        "shared cm", "sharedcm", "total cm", "shared dna cm", "centimorgans",
        "shared_cm", "total shared cm", "shared dna", "cm", "total centimorgans",
        "shared length",
    ),
    "percent": ("shared dna percent", "shared dna %", "percent", "% shared dna", "percentage"),
    "seg_count": (
        "shared segments", "segments", "num segments", "number of shared segments",
        "segment count", "shared segment count",
    ),
    "largest_cm": (
        "longest cm", "largest segment", "longest segment", "largest cm",
        "longest block", "largest seg",
    ),
    "predicted": (
        "predicted relationship", "relationship", "relationship range",
        "estimated relationship", "predicted", "range",
    ),
    "side": ("side", "parent side", "maternal or paternal"),
    "sex": ("sex", "gender"),
    "birth_year": ("birth year", "birthyear", "born", "year of birth"),
    "notes": ("notes", "note", "comment", "comments"),
    "shared_x_cm": ("x cm", "shared x", "x-dna", "shared x cm"),
    "true_relationship": ("true relationship", "known relationship", "actual relationship"),
    # segment files
    "chrom": ("chromosome", "chr", "chromosome number"),
    "start_bp": ("start position", "start location", "start", "start point", "b37 start"),
    "end_bp": ("end position", "end location", "end", "end point", "b37 end"),
    "seg_cm": ("centimorgans", "cm", "shared cm", "length cm", "genetic distance"),
    "snps": ("matching snps", "snps", "num snps", "number of snps", "snp count"),
    # shared-match files
    "name_a": ("match name a", "name a", "match a", "person a", "source match"),
    "name_b": ("match name b", "name b", "match b", "person b", "shared match"),
}


@dataclass
class ImportReport:
    path: str
    source: str = ""
    rows_read: int = 0
    imported: int = 0
    skipped: int = 0
    unmatched_names: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    columns: Dict[str, str] = field(default_factory=dict)

    def summary(self) -> str:
        bits = [f"{self.imported} imported", f"{self.rows_read} rows read"]
        if self.skipped:
            bits.append(f"{self.skipped} skipped")
        if self.unmatched_names:
            bits.append(f"{len(self.unmatched_names)} names not found in the match list")
        return ", ".join(bits)


def normalize_name(name: str) -> str:
    """Collapse a display name for joining across files.

    Match lists and segment exports come from different tools and disagree
    about punctuation, case, and whitespace, so joins have to be fuzzy at
    this level or nothing lines up.
    """
    if not name:
        return ""
    n = name.strip().lower()
    n = re.sub(r"[’']", "", n)
    n = re.sub(r"[^a-z0-9 ]+", " ", n)
    n = re.sub(r"\s+", " ", n).strip()
    return n


def _canonical(header_cell: str) -> Optional[str]:
    cell = re.sub(r"[^a-z0-9 %]+", " ", header_cell.strip().lower())
    cell = re.sub(r"\s+", " ", cell).strip()
    for field_name, aliases in _COLUMN_ALIASES.items():
        for alias in aliases:
            if cell == alias:
                return field_name
    for field_name, aliases in _COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in cell:
                return field_name
    return None


def _open_rows(path: str) -> Iterator[Dict[str, str]]:
    """Yield dict rows from a CSV/TSV file, possibly inside a zip."""
    if path.endswith(".zip") or zipfile.is_zipfile(path):
        zf = zipfile.ZipFile(path)
        members = [n for n in zf.namelist() if n.lower().endswith((".csv", ".tsv", ".txt"))]
        if not members:
            raise ValueError(f"{path}: no delimited file inside the archive")
        handle: Any = io.TextIOWrapper(zf.open(members[0]), encoding="utf-8-sig", errors="replace")
    else:
        handle = open(path, "r", encoding="utf-8-sig", errors="replace", newline="")
    with handle as fh:
        sample = fh.read(8192)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        for row in reader:
            yield {(k or ""): (v or "") for k, v in row.items()}


def _map_columns(row: Dict[str, str]) -> Dict[str, str]:
    """Map canonical field -> actual column name for this file."""
    out: Dict[str, str] = {}
    for col in row:
        canon = _canonical(col)
        if canon and canon not in out:
            out[canon] = col
    return out


def _num(value: str) -> Optional[float]:
    if value is None:
        return None
    v = str(value).strip().replace(",", "").replace("%", "")
    v = re.sub(r"\s*cM\s*$", "", v, flags=re.I)
    if not v or v in ("-", "n/a", "na", "none"):
        return None
    try:
        return float(v)
    except ValueError:
        m = re.search(r"-?\d+(?:\.\d+)?", v)
        return float(m.group()) if m else None


def _int(value: str) -> Optional[int]:
    f = _num(value)
    return int(f) if f is not None else None


# ---------------------------------------------------------------------------
# match lists


def import_matches(
    store: Store,
    kit_id: int,
    path: str,
    source: Optional[str] = None,
) -> ImportReport:
    """Load a match list for one kit.

    ``source`` names the platform the list came from.  It is not cosmetic:
    absence of a person from a list only means something when you compare
    lists from the *same* platform, so side inference refuses to work across
    sources.
    """
    src = source or _guess_source(path)
    report = ImportReport(path=path, source=src)
    rows = _open_rows(path)
    first = next(rows, None)
    if first is None:
        report.warnings.append("file contained no data rows")
        return report
    cols = _map_columns(first)
    report.columns = dict(cols)
    if "name" not in cols and "remote_id" not in cols:
        report.warnings.append(
            "could not find a name or id column; check the header row"
        )
        return report
    if "total_cm" not in cols and "percent" in cols:
        report.warnings.append(
            "no centimorgan column found; converting from percentage. "
            "Percentages count both chromosome copies, so anything at the "
            "sibling level or closer will be overstated -- prefer a cM export"
        )

    for row in _chain(first, rows):
        report.rows_read += 1
        name = (row.get(cols.get("name", ""), "") or "").strip()
        remote = (row.get(cols.get("remote_id", ""), "") or "").strip()
        if not name and not remote:
            report.skipped += 1
            continue
        remote = remote or normalize_name(name)

        total_cm = _num(row.get(cols.get("total_cm", ""), "")) if "total_cm" in cols else None
        if total_cm is None and "percent" in cols:
            pct = _num(row.get(cols["percent"], ""))
            if pct is not None:
                total_cm = pct / 100.0 * PERCENT_TO_CM
        if total_cm is None:
            report.skipped += 1
            continue

        extra: Dict[str, Any] = {}
        if "true_relationship" in cols:
            extra["true_relationship"] = row.get(cols["true_relationship"], "")
        for col, value in row.items():
            if _canonical(col) is None and value:
                extra[col] = value

        store.upsert_match(
            kit_id,
            src,
            remote,
            name=name or remote,
            total_cm=total_cm,
            seg_count=_int(row.get(cols.get("seg_count", ""), "")),
            largest_cm=_num(row.get(cols.get("largest_cm", ""), "")),
            shared_x_cm=_num(row.get(cols.get("shared_x_cm", ""), "")),
            predicted=(row.get(cols.get("predicted", ""), "") or "").strip() or None,
            side=(row.get(cols.get("side", ""), "") or "").strip().lower() or None,
            sex=(row.get(cols.get("sex", ""), "") or "").strip()[:1].upper() or None,
            birth_year=_int(row.get(cols.get("birth_year", ""), "")),
            notes=(row.get(cols.get("notes", ""), "") or "").strip() or None,
            extra=extra or None,
        )
        report.imported += 1
    store.commit()
    return report


def _chain(first: Dict[str, str], rest: Iterable[Dict[str, str]]) -> Iterator[Dict[str, str]]:
    yield first
    yield from rest


def _guess_source(path: str) -> str:
    base = os.path.basename(path).lower()
    for needle, name in (
        ("23andme", "23andMe"), ("ancestry", "Ancestry"), ("ftdna", "FTDNA"),
        ("familytree", "FTDNA"), ("myheritage", "MyHeritage"), ("gedmatch", "GEDmatch"),
        ("livingdna", "LivingDNA"),
    ):
        if needle in base:
            return name
    return "unknown"


# ---------------------------------------------------------------------------
# segments


def import_segments(
    store: Store,
    kit_id: int,
    path: str,
    source: Optional[str] = None,
    create_missing: bool = True,
) -> ImportReport:
    """Load a chromosome-browser style segment export.

    Rows are joined to already-imported matches by normalised name.  Names
    that do not resolve are reported rather than silently dropped, because a
    systematic mismatch (a different export naming convention) would
    otherwise look like a match with no segments.
    """
    src = source or _guess_source(path)
    report = ImportReport(path=path, source=src)
    by_name: Dict[str, int] = {}
    for m in store.matches(kit_id):
        by_name[normalize_name(m["name"] or "")] = m["id"]
        if m["remote_id"]:
            by_name.setdefault(normalize_name(m["remote_id"]), m["id"])

    pending: Dict[int, List[Tuple[Any, ...]]] = {}
    rows = _open_rows(path)
    first = next(rows, None)
    if first is None:
        report.warnings.append("file contained no data rows")
        return report
    cols = _map_columns(first)
    report.columns = dict(cols)
    needed = {"chrom", "start_bp", "end_bp"}
    if not needed.issubset(cols):
        report.warnings.append(
            f"missing required columns {sorted(needed - set(cols))}; got {sorted(cols)}"
        )
        return report

    unmatched: set = set()
    for row in _chain(first, rows):
        report.rows_read += 1
        name = (row.get(cols.get("name", ""), "") or "").strip()
        key = normalize_name(name)
        mid = by_name.get(key)
        if mid is None:
            if not create_missing or not name:
                unmatched.add(name or "(blank)")
                report.skipped += 1
                continue
            mid = store.upsert_match(kit_id, src, key, name=name, total_cm=0.0)
            by_name[key] = mid
        chrom = normalize_chrom(row.get(cols["chrom"], ""))
        start = _int(row.get(cols["start_bp"], ""))
        end = _int(row.get(cols["end_bp"], ""))
        if chrom is None or start is None or end is None:
            report.skipped += 1
            continue
        if start > end:
            start, end = end, start
        pending.setdefault(mid, []).append(
            (
                chrom, start, end,
                _num(row.get(cols.get("seg_cm", ""), "")),
                _int(row.get(cols.get("snps", ""), "")),
            )
        )
        report.imported += 1

    for mid, segs in pending.items():
        store.replace_match_segments(mid, segs)
        # Where the match list gave no total, derive one from the segments.
        row = store.match_by_id(mid)
        if row and not row["total_cm"]:
            total = sum(s[3] or 0.0 for s in segs)
            largest = max((s[3] or 0.0 for s in segs), default=0.0)
            store.db.execute(
                "UPDATE match SET total_cm=?, seg_count=?, largest_cm=? WHERE id=?",
                (total, len(segs), largest, mid),
            )
    report.unmatched_names = sorted(unmatched)[:50]
    store.commit()
    return report


# ---------------------------------------------------------------------------
# shared matches


def import_shared_matches(
    store: Store,
    kit_id: int,
    path: str,
    kind: str = "icw",
) -> ImportReport:
    """Load in-common-with pairs, the input to clustering.

    ``kind`` distinguishes plain shared-match data ("these two both match
    you") from true triangulation ("these two also match *each other* on the
    same segment").  The difference is genealogically load-bearing: two
    people can both match you without being related to each other at all,
    once through your mother and once through your father.
    """
    report = ImportReport(path=path, source=kind)
    by_name: Dict[str, int] = {}
    for m in store.matches(kit_id):
        by_name[normalize_name(m["name"] or "")] = m["id"]
        if m["remote_id"]:
            by_name.setdefault(normalize_name(m["remote_id"]), m["id"])

    rows = _open_rows(path)
    first = next(rows, None)
    if first is None:
        report.warnings.append("file contained no data rows")
        return report
    cols = _map_columns(first)
    report.columns = dict(cols)
    if "name_a" not in cols or "name_b" not in cols:
        keys = list(first.keys())
        if len(keys) >= 2:
            cols["name_a"], cols["name_b"] = keys[0], keys[1]
            report.warnings.append(
                f"no recognised name columns; assuming the first two are the pair "
                f"({keys[0]!r}, {keys[1]!r})"
            )
        else:
            report.warnings.append("need two name columns for shared-match data")
            return report

    unmatched: set = set()
    for row in _chain(first, rows):
        report.rows_read += 1
        a_name = (row.get(cols["name_a"], "") or "").strip()
        b_name = (row.get(cols["name_b"], "") or "").strip()
        a = by_name.get(normalize_name(a_name))
        b = by_name.get(normalize_name(b_name))
        if a is None or b is None or a == b:
            for nm, mid in ((a_name, a), (b_name, b)):
                if mid is None and nm:
                    unmatched.add(nm)
            report.skipped += 1
            continue
        store.add_shared_match(kit_id, a, b, _num(row.get(cols.get("seg_cm", ""), "")), kind)
        report.imported += 1
    report.unmatched_names = sorted(unmatched)[:50]
    store.commit()
    return report

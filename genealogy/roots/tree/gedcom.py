"""GEDCOM 5.5.1 reading and writing.

GEDCOM (GEnealogical Data COMmunication) is the lineage-linked interchange
format every genealogy program can read and write.  A file is a flat list of
records -- one per individual (INDI), one per family (FAM) -- joined by
pointers: an individual says which families it is a *child* in (FAMC) and
which it is a *spouse* in (FAMS); a family names its HUSB, WIFE and CHIL.
There is no "parent of" link; parenthood is always mediated by a family
record, which is what lets the format express a remarriage, a child of an
unknown mother, or an adoption without inventing new relationship types.

Version 5.5.1 dates from 1999 and is still the near-universal export target.
Version 7.0 exists and is better specified, but as of today Ancestry,
MyHeritage, FamilySearch, Gramps and RootsMagic all *export* 5.5.1 by
default, so that is what lands in a user's downloads folder and that is what
we read and write.  Everything here is deliberately lenient on input and
conservative on output: we accept anything we can make sense of, and emit
only structures the spec blesses.

The real-world messiness this parser is built to survive, all of it observed
in files from the sites above:

* Text values split across CONC (join with no separator) and CONT (join with
  a newline) sub-records.  Notes and long place names always arrive this way,
  because the spec caps a physical line at 255 characters.
* Encodings.  The header declares CHAR, but the declaration lies often enough
  that we simply try UTF-8, then UTF-16 if there is a byte-order mark, then
  fall back to latin-1, which cannot fail.  A file never refuses to load.
  ANSEL -- the pre-Unicode library encoding the 5.5.1 spec actually mandates
  -- we detect but cannot decode faithfully, so we warn: accented characters
  in such a file will be wrong, and the fix is to re-export as UTF-8.
* Names as ``John Henry /Smith/``, where the slashes delimit the surname, but
  also bare ``John Smith`` with no slashes at all, and GIVN/SURN sub-tags
  that contradict the NAME line (the sub-tags win -- they were entered as
  structured fields, the NAME line is a rendering of them).
* Dates that are not dates: ``ABT 1832``, ``BEF 1900``, ``BET 1830 AND
  1840``, ``EST 1700``, dual-dated ``1845/46`` from the Julian-to-Gregorian
  changeover.  We keep the raw string verbatim, because the qualifier *is*
  the evidence, and separately extract a bare year for arithmetic.
* One-directional links.  Some exporters write CHIL on the family but omit
  FAMC on the child, or vice versa.  We cross-link after parsing so both
  directions are populated, and synthesise stub records for pointers that go
  nowhere rather than dropping the relationship.
* Records with no xref, and xrefs used twice.  Both get a warning; the first
  definition of a repeated xref wins.

Two representational choices worth knowing.  Inside this module an xref is
the literal delimited token from the file, ``@I1@``, because pointer values
are then matched by plain string equality.  In the store the delimiters are
stripped, so the column reads as ``I1`` and can be prefixed with a
source-file slug when two uploads collide (they always do -- every exporter
starts numbering at @I1@).  `export_gedcom` puts the delimiters back.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

#: Physical line length at which a value is folded onto a CONC line.  The
#: spec's own limit is 255 including level and tag; 200 leaves room for both
#: and matches what most exporters do in practice.
CONC_LIMIT = 200

#: Years outside this range are parse artefacts (a street number, a page
#: reference, a mistyped date), not birth years.
MIN_YEAR = 1000
MAX_YEAR = 2200

_LINE_RE = re.compile(
    r"""^\s*
        (\d+)\s+                 # level
        (?:(@[^@\s]*@)\s+)?      # optional xref of the record being defined
        ([A-Za-z0-9_]+)          # tag
        (?:\s(.*))?$             # optional value, one space delimiter
    """,
    re.VERBOSE,
)

_YEAR_RE = re.compile(r"\d{4}")
_POINTER_RE = re.compile(r"^@[^@\s]+@$")
_PUNCT_RE = re.compile(r"[^a-z0-9 ]+")

#: Placeholders genealogists type when a name is unknown.  Treating these as
#: empty stops "Unknown Smith" matching "Unknown Jones" on given name.
_NULL_NAMES = {"", "unknown", "unk", "nn", "nn nn", "n n", "none", "no name", "living"}

#: Tags we read off an individual or family and model as first-class fields.
#: Anything else with a scalar value lands in `extra` instead of being lost.
_INDI_KNOWN = {
    "NAME", "SEX", "BIRT", "CHR", "DEAT", "BURI", "FAMC", "FAMS", "NOTE",
}
_FAM_KNOWN = {"HUSB", "WIFE", "CHIL", "MARR", "NOTE"}

#: Top-level record types we understand.  Everything else is counted, not
#: parsed -- an unrecognised SOUR or OBJE record must never stop an import.
_TOP_KNOWN = {"HEAD", "INDI", "FAM", "NOTE", "TRLR", "SUBM", "SUBN"}


@dataclass
class GedcomIndividual:
    """One INDI record.

    `birth_date` is the raw GEDCOM string ("ABT 1832"); `birth_year` is the
    integer a caller can do arithmetic with, or None when the date carries no
    usable year.  Keeping both means a report can quote the qualifier while a
    filter still works.
    """

    xref: str
    given: str = ""
    surname: str = ""
    sex: Optional[str] = None
    birth_date: Optional[str] = None
    birth_place: Optional[str] = None
    death_date: Optional[str] = None
    death_place: Optional[str] = None
    birth_year: Optional[int] = None
    death_year: Optional[int] = None
    notes: str = ""
    famc: List[str] = field(default_factory=list)
    fams: List[str] = field(default_factory=list)
    extra: Dict[str, Any] = field(default_factory=dict)

    @property
    def display_name(self) -> str:
        parts = [p for p in (self.given.strip(), self.surname.strip()) if p]
        return " ".join(parts) if parts else "(unknown)"


@dataclass
class GedcomFamily:
    """One FAM record.

    5.5.1 has exactly two spouse slots, HUSB and WIFE, and no way to say
    "two mothers"; a same-sex couple round-trips only by putting one partner
    in each slot.  That is a limitation of the format, not a modelling choice
    made here, and it is why `extra` may carry a `spouse_slot_guessed` flag
    when we had to place a spouse without knowing their sex.
    """

    xref: str
    husb: Optional[str] = None
    wife: Optional[str] = None
    children: List[str] = field(default_factory=list)
    marr_date: Optional[str] = None
    marr_place: Optional[str] = None
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class GedcomFile:
    """The parsed contents of one .ged file."""

    path: str
    individuals: Dict[str, GedcomIndividual] = field(default_factory=dict)
    families: Dict[str, GedcomFamily] = field(default_factory=dict)
    charset: Optional[str] = None
    source_software: Optional[str] = None
    warnings: List[str] = field(default_factory=list)

    def summary(self) -> str:
        years = [i.birth_year for i in self.individuals.values() if i.birth_year]
        span = f"{min(years)}-{max(years)}" if years else "no dated births"
        named = sum(1 for i in self.individuals.values() if i.surname.strip())
        return (
            f"{os.path.basename(self.path)}: {len(self.individuals)} individuals "
            f"({named} with a surname), {len(self.families)} families, "
            f"births {span}, charset {self.charset or 'undeclared'}, "
            f"software {self.source_software or 'unknown'}, "
            f"{len(self.warnings)} warning(s)"
        )


@dataclass
class _Node:
    """A raw GEDCOM line plus its sub-lines, before interpretation."""

    level: int
    tag: str
    value: str = ""
    xref: Optional[str] = None
    children: List["_Node"] = field(default_factory=list)

    def first(self, tag: str) -> Optional["_Node"]:
        for c in self.children:
            if c.tag == tag:
                return c
        return None

    def value_of(self, tag: str) -> Optional[str]:
        n = self.first(tag)
        return n.value.strip() if n is not None and n.value.strip() else None

    def all_of(self, tag: str) -> List["_Node"]:
        return [c for c in self.children if c.tag == tag]


# -- reading ------------------------------------------------------------


def _read_text(path: str, warnings: List[str]) -> str:
    """Decode a GEDCOM file, never raising on encoding.

    Order matters: a byte-order mark is authoritative, UTF-8 is what modern
    exporters emit, and latin-1 decodes any byte sequence at all, so it is
    the guaranteed terminal fallback.  A file that is really ANSEL will come
    through as latin-1 mojibake in the accented characters only; the header
    check in `parse_gedcom` warns about that separately.
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw.startswith(b"\xef\xbb\xbf"):
        raw = raw[3:]
    elif raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        try:
            return raw.decode("utf-16")
        except UnicodeDecodeError:
            warnings.append("file has a UTF-16 BOM but is not valid UTF-16; read as latin-1")
            return raw.decode("latin-1")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        warnings.append(
            "file is not valid UTF-8; decoded as latin-1, so non-ASCII "
            "characters may be wrong"
        )
        return raw.decode("latin-1")


def _parse_nodes(text: str, warnings: List[str]) -> List[_Node]:
    """Turn physical lines into a forest of level-0 records.

    CONC and CONT never become nodes: they are folded straight into the value
    of the line above them, which is the only sane representation for a
    consumer.  A CONC joins with nothing, a CONT joins with a newline.
    """
    roots: List[_Node] = []
    stack: List[_Node] = []
    bad = 0
    for lineno, raw_line in enumerate(text.replace("\r\n", "\n").replace("\r", "\n").split("\n"), 1):
        if not raw_line.strip():
            continue
        m = _LINE_RE.match(raw_line)
        if not m:
            bad += 1
            continue
        level = int(m.group(1))
        xref = m.group(2)
        tag = m.group(3).upper()
        value = m.group(4) or ""
        if tag in ("CONC", "CONT"):
            if not stack:
                bad += 1
                continue
            target = stack[min(level, len(stack)) - 1]
            target.value += ("\n" + value) if tag == "CONT" else value
            continue
        del stack[level:]
        node = _Node(level=level, tag=tag, value=value, xref=xref)
        if stack:
            stack[-1].children.append(node)
        else:
            roots.append(node)
            if level != 0:
                # A record whose level is not 0 with nothing above it: the
                # file is damaged, but the record itself may still be usable.
                bad += 1
        stack.append(node)
    if bad:
        warnings.append(f"{bad} line(s) could not be parsed and were skipped")
    return roots


def parse_year(raw: Optional[str]) -> Optional[int]:
    """Pull a usable year out of a GEDCOM date phrase.

    The first four-digit number in range wins, which gives the right answer
    for every common form: ``ABT 1832``, ``12 MAR 1845``, ``BEF 1900``, and
    ``BET 1830 AND 1840`` (the earliest bound, the conservative reading for a
    birth).  Dual dates such as ``1845/46`` yield 1845, the Julian year as
    written; the Gregorian alternative is only ever the following year, so
    nothing downstream is misled by more than twelve months.
    """
    if not raw:
        return None
    for m in _YEAR_RE.finditer(raw):
        year = int(m.group(0))
        if MIN_YEAR <= year <= MAX_YEAR:
            return year
    return None


def split_name(value: str) -> Tuple[str, str, str]:
    """Split a NAME value into (given, surname, suffix).

    The slash-delimited form is unambiguous.  Without slashes we take the
    last whitespace-delimited token as the surname, which is right for
    Western naming order and wrong for the minority of files that use family
    name first -- there is no way to tell from the value alone, and guessing
    Western order recovers far more names than it corrupts.  A lone token is
    treated as a given name, because a person recorded by one name is usually
    recorded by their first.
    """
    value = value.strip()
    if not value:
        return "", "", ""
    m = re.match(r"^(.*?)/([^/]*)/(.*)$", value)
    if m:
        return _tidy(m.group(1)), _tidy(m.group(2)), _tidy(m.group(3))
    tokens = value.split()
    if len(tokens) == 1:
        return tokens[0], "", ""
    return " ".join(tokens[:-1]), tokens[-1], ""


def _tidy(s: str) -> str:
    return " ".join(s.split()).strip().strip(",")


def _event(node: _Node) -> Tuple[Optional[str], Optional[str]]:
    return node.value_of("DATE"), node.value_of("PLAC")


def _collect_unmodelled(node: _Node, known: set, extra: Dict[str, Any]) -> None:
    """Keep scalar values from tags we do not model, rather than losing them.

    Only the value is kept, never the sub-structure: a full SOUR citation
    with its PAGE, DATA and QUAY children is a research artefact we have no
    schema for, and half of it in a JSON blob would be worse than an honest
    pointer back to the original file.
    """
    for child in node.children:
        if child.tag in known or not child.value.strip():
            continue
        key = child.tag.lower()
        val = child.value.strip()
        if key in extra:
            if isinstance(extra[key], list):
                extra[key].append(val)
            else:
                extra[key] = [extra[key], val]
        else:
            extra[key] = val


def _build_individual(node: _Node, xref: str, warnings: List[str]) -> GedcomIndividual:
    indi = GedcomIndividual(xref=xref)
    note_refs: List[str] = []
    birth: Dict[str, Optional[str]] = {"date": None, "place": None}
    death: Dict[str, Optional[str]] = {"date": None, "place": None}
    fallback_birth: Dict[str, Optional[str]] = {"date": None, "place": None}
    fallback_death: Dict[str, Optional[str]] = {"date": None, "place": None}
    seen_name = False

    for child in node.children:
        tag = child.tag
        if tag == "NAME":
            given, surname, suffix = split_name(child.value)
            structured_given = child.value_of("GIVN")
            structured_surname = child.value_of("SURN")
            if structured_given:
                given = _tidy(structured_given)
            if structured_surname:
                surname = _tidy(structured_surname)
            if not seen_name:
                indi.given, indi.surname = given, surname
                if suffix:
                    indi.extra["name_suffix"] = suffix
                seen_name = True
            else:
                indi.extra.setdefault("alt_names", []).append(
                    " ".join(p for p in (given, surname) if p)
                )
        elif tag == "SEX":
            sex = child.value.strip().upper()[:1]
            indi.sex = sex if sex in ("M", "F") else None
        elif tag in ("BIRT", "CHR", "DEAT", "BURI"):
            date, place = _event(child)
            target = {
                "BIRT": birth, "CHR": fallback_birth,
                "DEAT": death, "BURI": fallback_death,
            }[tag]
            if date and not target["date"]:
                target["date"] = date
            if place and not target["place"]:
                target["place"] = place
        elif tag == "FAMC":
            fam = child.value.strip()
            if fam:
                if fam not in indi.famc:
                    indi.famc.append(fam)
                pedi = child.value_of("PEDI")
                if pedi:
                    indi.extra.setdefault("pedi", {})[fam] = pedi.lower()
        elif tag == "FAMS":
            fam = child.value.strip()
            if fam and fam not in indi.fams:
                indi.fams.append(fam)
        elif tag == "NOTE":
            note_refs.append(child.value)

    # A christening or burial stands in for a missing birth or death.  Both
    # are usually within days or weeks of the event they proxy, which is far
    # better than nothing for sorting generations -- but the substitution is
    # recorded so a report can say "christened 1832" rather than asserting a
    # birth year the source never gave.
    if not birth["date"] and fallback_birth["date"]:
        birth["date"] = fallback_birth["date"]
        indi.extra["birth_date_from"] = "CHR"
    if not birth["place"] and fallback_birth["place"]:
        birth["place"] = fallback_birth["place"]
    if not death["date"] and fallback_death["date"]:
        death["date"] = fallback_death["date"]
        indi.extra["death_date_from"] = "BURI"
    if not death["place"] and fallback_death["place"]:
        death["place"] = fallback_death["place"]

    indi.birth_date, indi.birth_place = birth["date"], birth["place"]
    indi.death_date, indi.death_place = death["date"], death["place"]
    indi.birth_year = parse_year(indi.birth_date)
    indi.death_year = parse_year(indi.death_date)
    indi.extra["_note_refs"] = note_refs
    _collect_unmodelled(node, _INDI_KNOWN, indi.extra)
    return indi


def _build_family(node: _Node, xref: str, warnings: List[str]) -> GedcomFamily:
    fam = GedcomFamily(xref=xref)
    note_refs: List[str] = []
    for child in node.children:
        tag = child.tag
        if tag == "HUSB":
            if fam.husb is None:
                fam.husb = child.value.strip() or None
            else:
                warnings.append(f"{xref}: more than one HUSB, keeping {fam.husb}")
        elif tag == "WIFE":
            if fam.wife is None:
                fam.wife = child.value.strip() or None
            else:
                warnings.append(f"{xref}: more than one WIFE, keeping {fam.wife}")
        elif tag == "CHIL":
            ptr = child.value.strip()
            if ptr and ptr not in fam.children:
                fam.children.append(ptr)
        elif tag == "MARR":
            fam.marr_date, fam.marr_place = _event(child)
        elif tag == "NOTE":
            note_refs.append(child.value)
    fam.extra["_note_refs"] = note_refs
    _collect_unmodelled(node, _FAM_KNOWN, fam.extra)
    return fam


def _resolve_notes(
    holder_extra: Dict[str, Any], notes: Dict[str, str], warnings: List[str]
) -> str:
    """Turn NOTE values -- inline text or pointers -- into one string.

    Unresolvable pointers are dropped rather than kept as ``@N12@``: a
    dangling pointer is a broken export, and leaving the token in a note
    field only pollutes anything that later displays or searches it.
    """
    refs = holder_extra.pop("_note_refs", [])
    out: List[str] = []
    for ref in refs:
        text = ref.strip()
        if _POINTER_RE.match(text):
            resolved = notes.get(text)
            if resolved is None:
                warnings.append(f"note pointer {text} does not resolve")
                continue
            text = resolved
        if text:
            out.append(text)
    return "\n\n".join(out)


def _cross_link(gf: GedcomFile) -> None:
    """Make FAMC/FAMS and CHIL/HUSB/WIFE agree in both directions.

    Exports routinely state a relationship once.  Ancestry omits FAMC on
    children; some hand-built files list only FAMC and never CHIL.  A
    consumer that trusts one direction silently loses half the tree, so we
    take the union.  Pointers to records that do not exist become stubs: an
    unnamed person is still evidence that a person was there.
    """
    for indi in list(gf.individuals.values()):
        for fx in indi.famc + indi.fams:
            if fx not in gf.families:
                gf.families[fx] = GedcomFamily(xref=fx)
                gf.warnings.append(f"{indi.xref} points at missing family {fx}; stub created")
    for fam in list(gf.families.values()):
        for ix in [p for p in (fam.husb, fam.wife) if p] + list(fam.children):
            if ix not in gf.individuals:
                gf.individuals[ix] = GedcomIndividual(xref=ix)
                gf.warnings.append(f"{fam.xref} points at missing individual {ix}; stub created")

    for fam in gf.families.values():
        for ix in fam.children:
            indi = gf.individuals[ix]
            if fam.xref not in indi.famc:
                indi.famc.append(fam.xref)
        for ix in (fam.husb, fam.wife):
            if ix:
                indi = gf.individuals[ix]
                if fam.xref not in indi.fams:
                    indi.fams.append(fam.xref)

    for indi in gf.individuals.values():
        for fx in indi.famc:
            fam = gf.families[fx]
            if indi.xref not in fam.children:
                fam.children.append(indi.xref)
        for fx in indi.fams:
            fam = gf.families[fx]
            if indi.xref in (fam.husb, fam.wife):
                continue
            if fam.husb is None and indi.sex != "F":
                fam.husb = indi.xref
                if indi.sex is None:
                    fam.extra["spouse_slot_guessed"] = True
            elif fam.wife is None:
                fam.wife = indi.xref
                if indi.sex is None:
                    fam.extra["spouse_slot_guessed"] = True
            else:
                gf.warnings.append(
                    f"{indi.xref} claims FAMS {fx} but both spouse slots are taken"
                )


def parse_gedcom(path: str) -> GedcomFile:
    """Parse a GEDCOM file into individuals and families.

    Never raises on malformed content: unparseable lines, unknown record
    types, missing xrefs and dangling pointers all become entries in
    `GedcomFile.warnings` so that a caller can report what was lost while
    still importing everything that was intelligible.
    """
    warnings: List[str] = []
    text = _read_text(path, warnings)
    roots = _parse_nodes(text, warnings)
    gf = GedcomFile(path=path, warnings=warnings)

    notes: Dict[str, str] = {}
    unknown_types: Dict[str, int] = {}
    auto_indi = 0
    auto_fam = 0
    auto_note = 0

    for node in roots:
        tag = node.tag
        if tag == "HEAD":
            gf.charset = node.value_of("CHAR")
            sour = node.first("SOUR")
            if sour is not None:
                name = sour.value_of("NAME") or sour.value.strip()
                vers = sour.value_of("VERS")
                if name:
                    gf.source_software = f"{name} {vers}" if vers else name
            if gf.charset and gf.charset.upper().startswith("ANSEL"):
                warnings.append(
                    "header declares CHAR ANSEL, a legacy library encoding we cannot "
                    "decode faithfully; accented characters will be wrong -- re-export "
                    "the file as UTF-8 if names look mangled"
                )
        elif tag == "INDI":
            xref = node.xref
            if not xref:
                auto_indi += 1
                xref = f"@I_auto_{auto_indi}@"
                warnings.append(f"INDI record with no xref; assigned {xref}")
            if xref in gf.individuals:
                warnings.append(f"duplicate individual xref {xref}; keeping the first")
                continue
            gf.individuals[xref] = _build_individual(node, xref, warnings)
        elif tag == "FAM":
            xref = node.xref
            if not xref:
                auto_fam += 1
                xref = f"@F_auto_{auto_fam}@"
                warnings.append(f"FAM record with no xref; assigned {xref}")
            if xref in gf.families:
                warnings.append(f"duplicate family xref {xref}; keeping the first")
                continue
            gf.families[xref] = _build_family(node, xref, warnings)
        elif tag == "NOTE":
            xref = node.xref
            if not xref:
                auto_note += 1
                xref = f"@N_auto_{auto_note}@"
            notes[xref] = node.value.strip()
        elif tag not in _TOP_KNOWN:
            unknown_types[tag] = unknown_types.get(tag, 0) + 1

    if unknown_types:
        listed = ", ".join(f"{t} x{n}" for t, n in sorted(unknown_types.items()))
        warnings.append(f"ignored unrecognised top-level records: {listed}")

    for indi in gf.individuals.values():
        indi.notes = _resolve_notes(indi.extra, notes, warnings)
    for fam in gf.families.values():
        fam_notes = _resolve_notes(fam.extra, notes, warnings)
        if fam_notes:
            fam.extra["notes"] = fam_notes

    _cross_link(gf)
    return gf


# -- import into the store ----------------------------------------------


def _bare(xref: str) -> str:
    return xref.strip("@")


def _slug(filename: str) -> str:
    """Short, file-derived prefix used to disambiguate colliding xrefs."""
    stem = os.path.splitext(os.path.basename(filename))[0].lower()
    cleaned = re.sub(r"[^a-z0-9]+", "", stem)
    return (cleaned or "src")[:8]


def _rename_map(
    xrefs: Sequence[str], existing: Dict[str, Optional[str]], slug: str, source: str
) -> Dict[str, str]:
    """Map file xrefs onto store xrefs, renaming only on a real collision.

    Re-importing the same file keeps the same ids, so an import is
    idempotent.  A different file claiming an id the store already holds gets
    the slug prefix, because @I1@ from two uploads is two different people
    far more often than it is one.
    """
    out: Dict[str, str] = {}
    taken = set(existing)
    for xref in xrefs:
        bare = _bare(xref)
        if bare in existing and existing[bare] != source:
            candidate = f"{slug}_{bare}"
            n = 2
            while candidate in taken:
                candidate = f"{slug}{n}_{bare}"
                n += 1
        else:
            candidate = bare
        out[xref] = candidate
        taken.add(candidate)
    return out


def import_gedcom(store: Any, path: str, replace: bool = False) -> Dict[str, int]:
    """Parse a GEDCOM file and write it into the store's tree tables.

    With ``replace`` the existing tree is dropped first; otherwise the file is
    merged alongside whatever is already there, which is the usual case --
    a case built from several people's uploaded trees.  Merging means xref
    collisions, so ids that clash with a *different* source file are
    namespaced by a slug of this file's name, and every FAMC/FAMS/CHIL link
    is rewritten to match.  Nothing is deduplicated here: two records for one
    ancestor stay two records until a human confirms otherwise (see
    `merge_duplicates`).
    """
    gf = parse_gedcom(path)
    if replace:
        store.clear_tree()
    source = os.path.basename(path)
    slug = _slug(source)

    existing_indi = {row["xref"]: row["source_file"] for row in store.individuals()}
    existing_fam = {row["xref"]: row["source_file"] for row in store.families()}
    imap = _rename_map(list(gf.individuals), existing_indi, slug, source)
    fmap = _rename_map(list(gf.families), existing_fam, slug, source)

    for xref, indi in gf.individuals.items():
        extra = dict(indi.extra)
        pedi = extra.get("pedi")
        if pedi:
            extra["pedi"] = {fmap.get(k, _bare(k)): v for k, v in pedi.items()}
        store.upsert_individual(
            xref=imap[xref],
            given=indi.given or None,
            surname=indi.surname or None,
            sex=indi.sex,
            birth_date=indi.birth_date,
            birth_place=indi.birth_place,
            death_date=indi.death_date,
            death_place=indi.death_place,
            birth_year=indi.birth_year,
            death_year=indi.death_year,
            notes=indi.notes or None,
            source_file=source,
            extra=json.dumps(extra, default=str) if extra else None,
        )

    for xref, fam in gf.families.items():
        store.upsert_family(
            xref=fmap[xref],
            husb=imap.get(fam.husb) if fam.husb else None,
            wife=imap.get(fam.wife) if fam.wife else None,
            marr_date=fam.marr_date,
            marr_place=fam.marr_place,
            source_file=source,
        )

    child_rows = 0
    for xref, fam in gf.families.items():
        for ix in fam.children:
            pedi = gf.individuals[ix].extra.get("pedi") or {}
            store.add_child(fmap[xref], imap[ix], pedi.get(xref, "birth"))
            child_rows += 1

    store.commit()
    return {
        "individuals": len(gf.individuals),
        "families": len(gf.families),
        "children": child_rows,
        "warnings": len(gf.warnings),
    }


# -- export from the store ----------------------------------------------


def _emit(out: List[str], level: int, tag: str, value: Optional[str] = None) -> None:
    """Append one logical value as one or more physical GEDCOM lines.

    Embedded newlines become CONT sub-lines and over-long runs become CONC
    sub-lines, which is exactly what `parse_gedcom` folds back together, so
    long notes survive a round trip unchanged.
    """
    if value is None or value == "":
        out.append(f"{level} {tag}")
        return
    text = str(value).replace("\r\n", "\n").replace("\r", "\n")
    for i, physical in enumerate(text.split("\n")):
        chunks = [physical[j:j + CONC_LIMIT] for j in range(0, len(physical), CONC_LIMIT)] or [""]
        if i == 0:
            out.append(f"{level} {tag} {chunks[0]}" if chunks[0] else f"{level} {tag}")
        else:
            out.append(f"{level + 1} CONT {chunks[0]}" if chunks[0] else f"{level + 1} CONT")
        for chunk in chunks[1:]:
            out.append(f"{level + 1} CONC {chunk}")


def _file_xrefs(xrefs: Sequence[str], prefix: str) -> Dict[str, str]:
    """Sanitise store xrefs into pointer tokens old readers will accept.

    5.5.1 allows a limited character set and 20 characters; a slug-prefixed
    id can exceed both, and truncation can collide, so anything unsafe is
    replaced by a sequential id rather than silently mangled.
    """
    out: Dict[str, str] = {}
    used: set = set()
    for n, xref in enumerate(xrefs, 1):
        candidate = re.sub(r"[^A-Za-z0-9_]", "_", str(xref))[:20]
        if not candidate or candidate in used:
            candidate = f"{prefix}{n}"
        while candidate in used:
            candidate = f"{prefix}{n}_{len(used)}"
        used.add(candidate)
        out[str(xref)] = candidate
    return out


def export_gedcom(store: Any, path: str, submitter: str = "roots") -> Dict[str, int]:
    """Write the store's tree as a GEDCOM 5.5.1 file.

    The output is deliberately plain: no vendor extensions, no _UID tags, no
    embedded media.  It is meant to be readable by anything, including
    twenty-year-old desktop software, so lines are CRLF-terminated and every
    value is folded at `CONC_LIMIT`.  Data we hold that 5.5.1 has no place
    for -- the contents of the `extra` column beyond adoption pedigree --
    is not written; a GEDCOM export is a lossy view of the project, and the
    project file remains the record of truth.
    """
    individuals = store.individuals()
    families = store.families()
    children = store.children()

    imap = _file_xrefs([r["xref"] for r in individuals], "I")
    fmap = _file_xrefs([r["xref"] for r in families], "F")

    fam_children: Dict[str, List[Tuple[str, Optional[str]]]] = {}
    indi_famc: Dict[str, List[Tuple[str, Optional[str]]]] = {}
    for row in children:
        fam_children.setdefault(row["fam_xref"], []).append((row["indi_xref"], row["rel"]))
        indi_famc.setdefault(row["indi_xref"], []).append((row["fam_xref"], row["rel"]))
    indi_fams: Dict[str, List[str]] = {}
    for row in families:
        for slot in ("husb", "wife"):
            if row[slot]:
                indi_fams.setdefault(row[slot], []).append(row["xref"])

    out: List[str] = []
    out.append("0 HEAD")
    _emit(out, 1, "SOUR", submitter)
    _emit(out, 2, "NAME", "roots genealogy toolkit")
    out.append("1 DEST ANY")
    out.append("1 GEDC")
    out.append("2 VERS 5.5.1")
    out.append("2 FORM LINEAGE-LINKED")
    out.append("1 CHAR UTF-8")
    _emit(out, 1, "SUBM", "@SUBM1@")
    out.append("0 @SUBM1@ SUBM")
    _emit(out, 1, "NAME", submitter)

    for row in individuals:
        xid = imap[row["xref"]]
        out.append(f"0 @{xid}@ INDI")
        given = (row["given"] or "").strip()
        surname = (row["surname"] or "").strip()
        _emit(out, 1, "NAME", f"{given} /{surname}/".strip())
        if given:
            _emit(out, 2, "GIVN", given)
        if surname:
            _emit(out, 2, "SURN", surname)
        if row["sex"]:
            _emit(out, 1, "SEX", row["sex"])
        for tag, date, place in (
            ("BIRT", row["birth_date"], row["birth_place"]),
            ("DEAT", row["death_date"], row["death_place"]),
        ):
            if date or place:
                out.append(f"1 {tag}")
                if date:
                    _emit(out, 2, "DATE", date)
                if place:
                    _emit(out, 2, "PLAC", place)
        if row["notes"]:
            _emit(out, 1, "NOTE", row["notes"])
        for fam_xref, rel in indi_famc.get(row["xref"], []):
            if fam_xref not in fmap:
                continue
            _emit(out, 1, "FAMC", f"@{fmap[fam_xref]}@")
            if rel and rel != "birth":
                _emit(out, 2, "PEDI", rel)
        for fam_xref in indi_fams.get(row["xref"], []):
            _emit(out, 1, "FAMS", f"@{fmap[fam_xref]}@")

    child_rows = 0
    for row in families:
        fid = fmap[row["xref"]]
        out.append(f"0 @{fid}@ FAM")
        if row["husb"] and row["husb"] in imap:
            _emit(out, 1, "HUSB", f"@{imap[row['husb']]}@")
        if row["wife"] and row["wife"] in imap:
            _emit(out, 1, "WIFE", f"@{imap[row['wife']]}@")
        if row["marr_date"] or row["marr_place"]:
            out.append("1 MARR")
            if row["marr_date"]:
                _emit(out, 2, "DATE", row["marr_date"])
            if row["marr_place"]:
                _emit(out, 2, "PLAC", row["marr_place"])
        for indi_xref, _rel in fam_children.get(row["xref"], []):
            if indi_xref not in imap:
                continue
            _emit(out, 1, "CHIL", f"@{imap[indi_xref]}@")
            child_rows += 1
    out.append("0 TRLR")

    with open(path, "w", encoding="utf-8", newline="\r\n") as fh:
        fh.write("\n".join(out) + "\n")

    return {
        "individuals": len(individuals),
        "families": len(families),
        "children": child_rows,
        "lines": len(out),
    }


# -- duplicate detection ------------------------------------------------

_W_SURNAME = 0.45
_W_GIVEN = 0.30
_W_YEAR = 0.20
_W_SEX = 0.05


def _normalise(name: Optional[str]) -> str:
    if not name:
        return ""
    text = _PUNCT_RE.sub(" ", str(name).lower())
    text = " ".join(text.split())
    return "" if text in _NULL_NAMES else text


def _edit_distance(a: str, b: str) -> int:
    """Levenshtein distance, two rows of state.

    Written out rather than imported because the whole toolkit is standard
    library only, and because the alternative libraries all normalise
    differently, which would make scores here incomparable between installs.
    """
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        current = [i]
        for j, cb in enumerate(b, 1):
            current.append(
                min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (ca != cb))
            )
        previous = current
    return previous[-1]


def _similarity(a: str, b: str) -> float:
    if not a and not b:
        return 0.0
    longest = max(len(a), len(b))
    return 1.0 - _edit_distance(a, b) / longest


def _token_similarity(a: str, b: str) -> float:
    """Compare two single name tokens, treating an initial as a wildcard.

    "J" against "John" is the same evidence as "John" against "John" minus a
    little, because an initial is what a census taker wrote down, not a
    different name.
    """
    if a == b:
        return 1.0
    if len(a) == 1 or len(b) == 1:
        return 0.9 if a[0] == b[0] else 0.0
    return _similarity(a, b)


def _given_similarity(a: str, b: str) -> float:
    """Token-overlap score for given names, ignoring order and extra tokens.

    A middle name recorded in one tree and not the other is the normal state
    of affairs and is not evidence that the two are different people, so
    unmatched surplus tokens cost nothing.
    """
    ta, tb = a.split(), b.split()
    if not ta and not tb:
        return 0.5
    if not ta or not tb:
        return 0.3
    short, long_ = (ta, tb) if len(ta) <= len(tb) else (tb, ta)
    used: set = set()
    total = 0.0
    for token in short:
        best, best_j = 0.0, None
        for j, other in enumerate(long_):
            if j in used:
                continue
            score = _token_similarity(token, other)
            if score > best:
                best, best_j = score, j
        if best_j is not None:
            used.add(best_j)
        total += best
    return total / len(short)


def _year_component(a: Optional[int], b: Optional[int]) -> Optional[float]:
    """Score a birth-year pair, or None when the pair is disqualified."""
    if a is None or b is None:
        return 0.5
    gap = abs(a - b)
    if gap == 0:
        return 1.0
    if gap == 1:
        return 0.85
    if gap == 2:
        return 0.7
    if gap <= 5:
        return 0.3
    return None


def _score_pair(a: Any, b: Any) -> float:
    """Blend the evidence into 0..1, or 0.0 when the pair is disqualified."""
    sex_a, sex_b = (a["sex"] or "").upper()[:1], (b["sex"] or "").upper()[:1]
    if sex_a in ("M", "F") and sex_b in ("M", "F") and sex_a != sex_b:
        return 0.0
    year = _year_component(a["birth_year"], b["birth_year"])
    if year is None:
        return 0.0

    sa, sb = _normalise(a["surname"]), _normalise(b["surname"])
    if not sa and not sb:
        surname = 0.5
    elif not sa or not sb:
        surname = 0.2
    elif sa == sb:
        surname = 1.0
    else:
        sim = _similarity(sa, sb)
        surname = sim if sim >= 0.7 else 0.0

    given = _given_similarity(_normalise(a["given"]), _normalise(b["given"]))
    sex = 1.0 if (sex_a and sex_a == sex_b) else 0.5

    return _W_SURNAME * surname + _W_GIVEN * given + _W_YEAR * year + _W_SEX * sex


def merge_duplicates(store: Any, threshold: float = 0.85) -> List[Tuple[str, str, float]]:
    """Propose pairs of individuals that may be the same person.

    When several people upload their trees into one project, a shared
    ancestor appears once per tree; finding those repetitions is what makes
    the trees usable together.  What this function will not do is merge them.

    Automatic merging of genealogical records is a bad idea, and not merely a
    risky one.  A merge is destructive and effectively irreversible once
    later research is built on top of it; the classic failure -- two cousins
    named for the same grandfather, born in the same parish four years apart
    -- produces a person who is their own uncle and a DNA hypothesis that
    can never be reconciled with the segment data.  The cost is asymmetric:
    a missed duplicate leaves two nodes a researcher can still join later,
    while a wrong merge silently corrupts every descendant line and every
    relationship prediction computed from it.  Scores here are a triage
    queue, ordered so a human reviews the strongest evidence first, and the
    threshold is a display convenience rather than a decision boundary.

    Scoring weights surname equality most heavily (a shared surname is the
    single strongest cheap signal), then given-name similarity computed by
    token overlap with initials treated as wildcards, then birth years within
    two years, then agreeing sex.  Contradicting sex or birth years more than
    five years apart disqualify a pair outright.  Individuals with no surname
    in either record cannot reach a default threshold by construction; there
    is not enough there to justify a researcher's time.

    Candidates within a single source file are included, but treat them with
    extra scepticism: repeated names inside one tree are usually a father and
    son, not a duplicate.
    """
    rows = store.individuals()
    blocks: Dict[str, List[Any]] = {}
    for row in rows:
        key = (_normalise(row["surname"]) or " ")[:1]
        blocks.setdefault(key, []).append(row)

    pairs: List[Tuple[str, str, float]] = []
    for block in blocks.values():
        for i in range(len(block)):
            for j in range(i + 1, len(block)):
                score = _score_pair(block[i], block[j])
                if score >= threshold:
                    a, b = sorted((str(block[i]["xref"]), str(block[j]["xref"])))
                    pairs.append((a, b, round(score, 4)))
    pairs.sort(key=lambda p: (-p[2], p[0], p[1]))
    return pairs


__all__ = [
    "GedcomIndividual",
    "GedcomFamily",
    "GedcomFile",
    "parse_gedcom",
    "import_gedcom",
    "export_gedcom",
    "merge_duplicates",
    "parse_year",
    "split_name",
]

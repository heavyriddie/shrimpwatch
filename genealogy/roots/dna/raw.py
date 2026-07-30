"""Readers for consumer-DNA raw genotype downloads.

Supported out of the box: 23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA
(Family Finder), LivingDNA, and the generic "rsid/chromosome/position/
genotype" layout that most other tools emit.  Files may be plain text,
``.gz``, or ``.zip`` (the vendors all hand you a zip).

Two normalisations happen here and matter downstream:

* Chromosome labels become canonical strings ('1'..'22', 'X', 'Y', 'MT').
  Ancestry's numeric 23/24/25/26 coding is translated, with the
  pseudo-autosomal code 25 folded into X.

* Genotypes become allele-sorted uppercase strings: both ``GA`` and ``AG``
  are stored as ``AG``, so a genotype comparison is a string comparison.
  No-calls of every vendor's spelling ('--', '00', 'NN', '', 'I'/'D' pairs
  with a missing side) collapse to ``--``.

A word on strand: all of these vendors report the plus strand of the
reference assembly, so genotypes are directly comparable between vendors
without flipping.  That is an assumption, not a guarantee, and
`roots.dna.qc.concordance` exists to check it -- if two kits from the same
person disagree at more than a fraction of a percent of shared SNPs, treat
the merge as suspect rather than silently reconciling it.
"""

from __future__ import annotations

import gzip
import io
import os
import re
import zipfile
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, TextIO, Tuple

from ..genome import CHROM_BP, normalize_chrom

NO_CALL = "--"

Record = Tuple[str, int, Optional[str], str]  # chrom, pos, rsid, genotype


@dataclass
class RawFileInfo:
    """What we learned about a raw file before reading its body."""

    path: str
    vendor: str = "unknown"
    build: str = "37"
    header_lines: List[str] = field(default_factory=list)
    delimiter: str = "\t"
    columns: List[str] = field(default_factory=list)
    layout: str = "genotype"  # 'genotype' (one column) or 'alleles' (two)
    inner_name: Optional[str] = None  # member name, when inside a zip

    def summary(self) -> str:
        inner = f" ({self.inner_name})" if self.inner_name else ""
        return f"{os.path.basename(self.path)}{inner}: {self.vendor}, build {self.build}"


@dataclass
class KitStats:
    """Counters accumulated while streaming a file, used for QC and sexing."""

    total: int = 0
    called: int = 0
    no_call: int = 0
    het: int = 0
    indel: int = 0
    per_chrom: Dict[str, int] = field(default_factory=dict)
    y_called: int = 0
    y_total: int = 0
    x_called: int = 0
    x_het: int = 0
    duplicates: int = 0
    off_reference: int = 0

    @property
    def call_rate(self) -> float:
        return self.called / self.total if self.total else 0.0

    def infer_sex(self) -> Optional[str]:
        """Infer sex from Y call rate and X heterozygosity.

        Both signals are needed.  Some vendors omit the Y entirely for
        female kits (no rows at all), which looks identical to "Y rows
        present but all no-call"; X heterozygosity disambiguates, since a
        male has essentially no heterozygous X calls outside the
        pseudo-autosomal regions.
        """
        y_rate = self.y_called / self.y_total if self.y_total else 0.0
        x_het_rate = self.x_het / self.x_called if self.x_called else 0.0
        if self.y_total == 0 and self.x_called == 0:
            return None
        if x_het_rate > 0.10:
            return "F"
        if y_rate > 0.30 or (self.x_called > 1000 and x_het_rate < 0.02):
            return "M"
        if self.y_total > 100 and y_rate < 0.10:
            return "F"
        return None


# ---------------------------------------------------------------------------
# opening


def _open_text(path: str) -> Tuple[TextIO, Optional[str]]:
    """Open a raw file as text, transparently handling .gz and .zip."""
    if path.endswith(".zip") or zipfile.is_zipfile(path):
        zf = zipfile.ZipFile(path)
        members = [
            n for n in zf.namelist()
            if not n.endswith("/") and not n.startswith("__MACOSX")
        ]
        if not members:
            raise ValueError(f"{path}: zip archive is empty")
        # Prefer the largest text member; vendor zips often carry a readme.
        members.sort(key=lambda n: zf.getinfo(n).file_size, reverse=True)
        name = members[0]
        raw = zf.open(name)
        return io.TextIOWrapper(raw, encoding="utf-8", errors="replace"), name
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace"), None
    return open(path, "r", encoding="utf-8", errors="replace"), None


# ---------------------------------------------------------------------------
# sniffing


def inspect(path: str) -> RawFileInfo:
    """Read just the header of a raw file to identify vendor and layout."""
    info = RawFileInfo(path=path)
    fh, inner = _open_text(path)
    info.inner_name = inner
    try:
        header_blob: List[str] = []
        first_data: Optional[str] = None
        column_line: Optional[str] = None
        for line in fh:
            if line.startswith("#"):
                header_blob.append(line.rstrip("\n"))
                continue
            stripped = line.strip()
            if not stripped:
                continue
            if column_line is None and re.search(r"rsid|chromosome", stripped, re.I):
                column_line = stripped
                continue
            first_data = stripped
            break
    finally:
        fh.close()

    info.header_lines = header_blob[:40]
    blob = "\n".join(header_blob).lower()

    if "23andme" in blob:
        info.vendor = "23andMe"
    elif "ancestry" in blob:
        info.vendor = "AncestryDNA"
    elif "myheritage" in blob:
        info.vendor = "MyHeritage"
    elif "living dna" in blob or "livingdna" in blob:
        info.vendor = "LivingDNA"
    elif "ftdna" in blob or "family tree dna" in blob:
        info.vendor = "FTDNA"

    sample = column_line or first_data or ""
    if sample.count(",") >= 3:
        info.delimiter = ","
    elif "\t" in sample:
        info.delimiter = "\t"
    else:
        info.delimiter = None  # any run of whitespace

    if column_line:
        cols = [c.strip().strip('"').lower() for c in _split(column_line, info.delimiter)]
        info.columns = cols
        if "allele1" in cols or "allele 1" in cols:
            info.layout = "alleles"
        if info.vendor == "unknown":
            if "result" in cols:
                # FTDNA and MyHeritage both use RESULT; MyHeritage quotes.
                info.vendor = "MyHeritage" if '"' in column_line else "FTDNA"
            elif "allele1" in cols:
                info.vendor = "AncestryDNA"
            elif "genotype" in cols:
                info.vendor = "23andMe"
    elif first_data:
        n = len(_split(first_data, info.delimiter))
        info.layout = "alleles" if n >= 5 else "genotype"

    m = re.search(r"(grch|build|reference human assembly)\D{0,12}(3[678])", blob)
    if m:
        info.build = "38" if m.group(2) == "38" else "37"
    else:
        info.build = ""  # unknown; caller may sniff from coordinates
    return info


def _split(line: str, delimiter: Optional[str]) -> List[str]:
    if delimiter is None:
        return line.split()
    return line.split(delimiter)


# ---------------------------------------------------------------------------
# parsing


def normalize_genotype(tokens: List[str]) -> str:
    """Turn vendor allele tokens into a canonical genotype string."""
    alleles: List[str] = []
    for tok in tokens:
        t = tok.strip().strip('"').upper()
        for ch in t:
            if ch in "ACGTID":
                alleles.append(ch)
            elif ch in "-0N":
                alleles.append("-")
    if not alleles:
        return NO_CALL
    if any(a == "-" for a in alleles):
        # Partial calls are not usable for IBD; treat the whole site as missing.
        return NO_CALL
    if len(alleles) == 1:
        return alleles[0]  # hemizygous (male X/Y, mitochondrial)
    return "".join(sorted(alleles[:2]))


def read(path: str, info: Optional[RawFileInfo] = None) -> Tuple[RawFileInfo, Iterator[Record]]:
    """Open a raw file and return ``(info, record_iterator)``.

    The iterator is lazy and single-pass; consume it before the caller
    inspects mutable stats.  Use `read_with_stats` when you want both.
    """
    info = info or inspect(path)
    fh, _ = _open_text(path)

    def gen() -> Iterator[Record]:
        try:
            delim = info.delimiter
            for line in fh:
                if not line or line[0] == "#":
                    continue
                s = line.strip()
                if not s:
                    continue
                parts = [p.strip().strip('"') for p in _split(s, delim)]
                if len(parts) < 4:
                    continue
                rsid = parts[0]
                if rsid.lower() in ("rsid", "rs id", "snp"):
                    continue  # column header
                chrom = normalize_chrom(parts[1])
                if chrom is None:
                    continue
                try:
                    pos = int(parts[2])
                except ValueError:
                    continue
                gt = normalize_genotype(parts[3:5] if len(parts) >= 5 else parts[3:4])
                # Mitochondrial and Y calls are hemizygous by definition, but
                # two-column formats have nowhere to say so and write the
                # allele twice. Collapse them so the same site read from two
                # vendors compares equal instead of looking like a conflict.
                if chrom in ("Y", "MT") and len(gt) == 2 and gt[0] == gt[1]:
                    gt = gt[0]
                yield chrom, pos, (rsid or None), gt
        finally:
            fh.close()

    return info, gen()


def read_with_stats(path: str, info: Optional[RawFileInfo] = None) -> Tuple[RawFileInfo, Iterator[Record], KitStats]:
    """Like `read`, but the returned `KitStats` fills in as you consume."""
    info, base = read(path, info)
    stats = KitStats()
    seen: set = set()

    def gen() -> Iterator[Record]:
        for chrom, pos, rsid, gt in base:
            stats.total += 1
            key = (chrom, pos)
            if key in seen:
                stats.duplicates += 1
                continue
            seen.add(key)
            stats.per_chrom[chrom] = stats.per_chrom.get(chrom, 0) + 1
            if chrom == "Y":
                stats.y_total += 1
            if gt == NO_CALL:
                stats.no_call += 1
            else:
                stats.called += 1
                if chrom == "Y":
                    stats.y_called += 1
                if chrom == "X":
                    stats.x_called += 1
                    if len(gt) == 2 and gt[0] != gt[1]:
                        stats.x_het += 1
                if len(gt) == 2 and gt[0] != gt[1]:
                    stats.het += 1
                if "I" in gt or "D" in gt:
                    stats.indel += 1
            yield chrom, pos, rsid, gt

    return info, gen(), stats


def resolve_build(info: RawFileInfo, sample: List[Tuple[str, int]]) -> str:
    """Settle on an assembly, preferring the file header over coordinates."""
    if info.build in ("37", "38"):
        return info.build
    from ..genome import detect_build

    return detect_build(sample)


def off_reference_count(records: List[Record], build: str) -> int:
    """How many positions fall past the end of their chromosome.

    A non-trivial count means the build guess is wrong.
    """
    lengths = CHROM_BP[build]
    bad = 0
    for chrom, pos, _rsid, _gt in records:
        limit = lengths.get(chrom)
        if limit and pos > limit:
            bad += 1
    return bad

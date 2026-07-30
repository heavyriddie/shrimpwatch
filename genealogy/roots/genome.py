"""Reference-genome geometry and genetic (recombination) maps.

Everything downstream that talks about "centimorgans" comes through here.
Two things matter:

1.  Physical geometry (base pairs) differs between GRCh37 (build 37 / hg19)
    and GRCh38.  Every consumer testing company currently ships raw data on
    build 37, but files do exist on 38, so positions carry a build tag.

2.  Genetic distance is *not* proportional to physical distance.  Real
    recombination rates vary by an order of magnitude along a chromosome
    (cold centromeres, hot subtelomeric regions).  If you have a real
    recombination map, load it -- see `GeneticMap.load`.  The built-in
    fallback is a uniform-rate approximation which is fine for sanity checks
    and for the meiosis simulator (which only needs total genetic lengths),
    but is *not* accurate enough to quote segment sizes to a genealogist.
    `GeneticMap.is_approximate` tells you which one you are holding, and the
    CLI prints a warning whenever an approximate map produces a number that
    ends up in a report.
"""

from __future__ import annotations

import bisect
import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

AUTOSOMES: Tuple[str, ...] = tuple(str(i) for i in range(1, 23))
ALL_CHROMS: Tuple[str, ...] = AUTOSOMES + ("X", "Y", "MT")

# Physical lengths, in base pairs.
CHROM_BP: Dict[str, Dict[str, int]] = {
    "37": {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
        "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
        "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
        "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
        "MT": 16569,
    },
    "38": {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
        "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
        "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
        "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
        "MT": 16569,
    },
}

# Sex-averaged genetic lengths in centimorgans, rounded from the standard
# published linkage maps.  These are approximations good to ~1%; they set the
# scale for the simulator and for the fallback map.  Autosomal total is about
# 3545 cM here, against the ~3485 cM that the testing companies quote for
# their reporting region (they exclude some low-confidence terminal and
# centromeric stretches), so simulated totals run a shade high.
CHROM_CM: Dict[str, float] = {
    "1": 281.5, "2": 263.7, "3": 224.2, "4": 214.5, "5": 209.4, "6": 194.1,
    "7": 187.0, "8": 169.0, "9": 166.4, "10": 181.1, "11": 158.2, "12": 174.7,
    "13": 125.7, "14": 120.2, "15": 141.9, "16": 134.0, "17": 128.4,
    "18": 117.7, "19": 107.7, "20": 108.3, "21": 62.8, "22": 74.1,
    # X recombines only in female meiosis; this is the female map length.
    "X": 180.0,
}

#: What the consumer testing industry calls "the whole genome" when it reports
#: a shared-cM total.  Used only for reporting ratios, never for arithmetic on
#: actual segments.
INDUSTRY_AUTOSOMAL_CM = 3485.0

AUTOSOMAL_CM_TOTAL = sum(CHROM_CM[c] for c in AUTOSOMES)


def normalize_chrom(raw: str) -> Optional[str]:
    """Map any vendor's chromosome spelling onto our canonical labels.

    Ancestry uses numeric codes through 26; 23andMe and MyHeritage use
    letters; some tools prefix ``chr``.  Returns None for anything we do not
    model (pseudo-autosomal regions are folded into X).
    """
    if raw is None:
        return None
    c = str(raw).strip().upper()
    if c.startswith("CHR"):
        c = c[3:]
    if c in ("23", "X"):
        return "X"
    if c in ("24", "Y"):
        return "Y"
    if c in ("25", "XY", "PAR"):
        # Pseudo-autosomal: physically on X and Y, reported separately by
        # Ancestry.  We fold it into X, which is where its coordinates live.
        return "X"
    if c in ("26", "MT", "M", "MITO"):
        return "MT"
    if c.isdigit() and 1 <= int(c) <= 22:
        return str(int(c))
    return None


@dataclass
class _ChromMap:
    """Interpolation table for one chromosome: sorted bp -> cumulative cM."""

    positions: List[int] = field(default_factory=list)
    cms: List[float] = field(default_factory=list)

    def cm_at(self, pos: int) -> float:
        if not self.positions:
            return 0.0
        i = bisect.bisect_left(self.positions, pos)
        if i == 0:
            return self.cms[0]
        if i >= len(self.positions):
            return self.cms[-1]
        p0, p1 = self.positions[i - 1], self.positions[i]
        c0, c1 = self.cms[i - 1], self.cms[i]
        if p1 == p0:
            return c1
        return c0 + (c1 - c0) * (pos - p0) / (p1 - p0)


class GeneticMap:
    """bp <-> cM conversion, either interpolated from a real map or uniform."""

    def __init__(self, build: str = "37", source: str = "linear-approximation"):
        self.build = build
        self.source = source
        self._maps: Dict[str, _ChromMap] = {}

    # -- construction ---------------------------------------------------

    @classmethod
    def linear(cls, build: str = "37") -> "GeneticMap":
        """Uniform recombination rate within each chromosome."""
        gm = cls(build=build, source="linear-approximation")
        for chrom, cm_len in CHROM_CM.items():
            bp_len = CHROM_BP[build].get(chrom)
            if not bp_len:
                continue
            m = _ChromMap(positions=[1, bp_len], cms=[0.0, cm_len])
            gm._maps[chrom] = m
        return gm

    @classmethod
    def load(cls, path: str, build: str = "37") -> "GeneticMap":
        """Load a real recombination map.

        Accepts either a directory of per-chromosome files (the usual HapMap
        layout, ``genetic_map_GRCh37_chr1.txt`` and friends) or a single file
        carrying a chromosome column.  Recognised column layouts:

        * HapMap: ``Chromosome Position(bp) Rate(cM/Mb) Map(cM)``
        * SHAPEIT/Impute: ``position COMBINED_rate(cM/Mb) Genetic_Map(cM)``
        * Bare two-column ``position cM``
        """
        gm = cls(build=build, source=os.path.abspath(path))
        if os.path.isdir(path):
            for name in sorted(os.listdir(path)):
                m = re.search(r"chr([0-9XYMT]+)", name, re.I)
                if not m:
                    continue
                chrom = normalize_chrom(m.group(1))
                if not chrom:
                    continue
                gm._load_file(os.path.join(path, name), default_chrom=chrom)
        else:
            gm._load_file(path, default_chrom=None)
        if not gm._maps:
            raise ValueError(f"no usable genetic map data found in {path!r}")
        return gm

    def _load_file(self, path: str, default_chrom: Optional[str]) -> None:
        rows: Dict[str, List[Tuple[int, float]]] = {}
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            header: Optional[List[str]] = None
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = re.split(r"[\t,\s]+", line)
                if header is None and not _looks_numeric(parts[0]):
                    header = [p.lower() for p in parts]
                    continue
                chrom, pos, cm = self._parse_map_row(parts, header, default_chrom)
                if chrom is None or pos is None or cm is None:
                    continue
                rows.setdefault(chrom, []).append((pos, cm))
        for chrom, pairs in rows.items():
            pairs.sort()
            m = self._maps.setdefault(chrom, _ChromMap())
            # Merge, keeping the map monotone; duplicate positions collapse.
            merged = dict(zip(m.positions, m.cms))
            merged.update(dict(pairs))
            keys = sorted(merged)
            m.positions = keys
            m.cms = [merged[k] for k in keys]

    @staticmethod
    def _parse_map_row(
        parts: Sequence[str],
        header: Optional[Sequence[str]],
        default_chrom: Optional[str],
    ) -> Tuple[Optional[str], Optional[int], Optional[float]]:
        def idx(*names: str) -> Optional[int]:
            if not header:
                return None
            for i, h in enumerate(header):
                for n in names:
                    if n in h:
                        return i
            return None

        c_i = idx("chr")
        p_i = idx("position", "pos", "bp")
        m_i = idx("map(cm)", "genetic_map", "map", "cm")
        try:
            if p_i is not None and m_i is not None:
                chrom = normalize_chrom(parts[c_i]) if c_i is not None else default_chrom
                return chrom, int(float(parts[p_i])), float(parts[m_i])
            # Positional fallbacks by column count.
            if len(parts) >= 4:
                return normalize_chrom(parts[0]) or default_chrom, int(float(parts[1])), float(parts[3])
            if len(parts) == 3:
                return default_chrom, int(float(parts[0])), float(parts[2])
            if len(parts) == 2:
                return default_chrom, int(float(parts[0])), float(parts[1])
        except (ValueError, IndexError):
            return None, None, None
        return None, None, None

    # -- queries --------------------------------------------------------

    @property
    def is_approximate(self) -> bool:
        return self.source == "linear-approximation"

    def has(self, chrom: str) -> bool:
        return chrom in self._maps

    def cm_at(self, chrom: str, pos: int) -> float:
        m = self._maps.get(chrom)
        if m is None:
            return 0.0
        return m.cm_at(pos)

    def length_cm(self, chrom: str, start_bp: int, end_bp: int) -> float:
        """Genetic length of a physical interval, never negative."""
        if end_bp < start_bp:
            start_bp, end_bp = end_bp, start_bp
        return max(0.0, self.cm_at(chrom, end_bp) - self.cm_at(chrom, start_bp))

    def chrom_cm(self, chrom: str) -> float:
        m = self._maps.get(chrom)
        if m is None or not m.cms:
            return CHROM_CM.get(chrom, 0.0)
        return m.cms[-1] - m.cms[0]

    def bp_at_cm(self, chrom: str, cm: float) -> int:
        """Inverse lookup, for turning simulated segments back into coordinates."""
        m = self._maps.get(chrom)
        if m is None or not m.cms:
            return 0
        i = bisect.bisect_left(m.cms, cm)
        if i == 0:
            return m.positions[0]
        if i >= len(m.cms):
            return m.positions[-1]
        c0, c1 = m.cms[i - 1], m.cms[i]
        p0, p1 = m.positions[i - 1], m.positions[i]
        if c1 == c0:
            return p1
        return int(p0 + (p1 - p0) * (cm - c0) / (c1 - c0))

    def autosomal_total_cm(self) -> float:
        return sum(self.chrom_cm(c) for c in AUTOSOMES)

    def describe(self) -> str:
        kind = "approximate (uniform rate)" if self.is_approximate else self.source
        return f"build {self.build}, {kind}, autosomal total {self.autosomal_total_cm():.0f} cM"


def _looks_numeric(tok: str) -> bool:
    try:
        float(tok)
        return True
    except ValueError:
        return False


def detect_build(positions: Iterable[Tuple[str, int]]) -> str:
    """Guess the assembly from the largest coordinate seen per chromosome.

    Raw files usually say which build they are on in a header comment; this is
    the backstop for when they do not.  A position past the GRCh38 end of a
    chromosome can only be build 37, and vice versa.  We tally votes and
    return the winner, defaulting to 37 because that is what every consumer
    vendor currently ships.
    """
    max_pos: Dict[str, int] = {}
    for chrom, pos in positions:
        if chrom in max_pos:
            if pos > max_pos[chrom]:
                max_pos[chrom] = pos
        else:
            max_pos[chrom] = pos
    votes = {"37": 0, "38": 0}
    for chrom, pos in max_pos.items():
        b37 = CHROM_BP["37"].get(chrom)
        b38 = CHROM_BP["38"].get(chrom)
        if not b37 or not b38:
            continue
        if pos > b38 and pos <= b37:
            votes["37"] += 1
        elif pos > b37 and pos <= b38:
            votes["38"] += 1
    if votes["38"] > votes["37"]:
        return "38"
    return "37"

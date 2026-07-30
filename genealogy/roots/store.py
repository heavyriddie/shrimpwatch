"""SQLite-backed project store.

One file holds everything about an investigation: raw genotypes from every
kit, imported match lists, computed IBD segments, phasing results, the family
tree, and cached simulator output.  Nothing here talks to the network, and
nothing leaves the file.

The genotype table is the only large one -- roughly 600k to 1.8M rows per kit
depending on chip version.  It is declared WITHOUT ROWID with a
``(kit_id, chrom, pos)`` primary key so that the two-kit joins that drive
IBD detection and phasing walk a single clustered index in position order.
"""

from __future__ import annotations

import json
import os
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

SCHEMA_VERSION = 1

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT
);

-- A "person" is a human; a "kit" is one test of one person at one vendor.
-- Two kits can and should map to the same person (23andMe + Ancestry).
CREATE TABLE IF NOT EXISTS person (
    id     INTEGER PRIMARY KEY,
    label  TEXT UNIQUE NOT NULL,
    sex    TEXT,
    birth_year INTEGER,
    tree_xref  TEXT,
    notes  TEXT
);

CREATE TABLE IF NOT EXISTS kit (
    id          INTEGER PRIMARY KEY,
    label       TEXT UNIQUE NOT NULL,
    person_id   INTEGER REFERENCES person(id),
    vendor      TEXT,
    build       TEXT,
    source_path TEXT,
    imported_at TEXT,
    snp_count   INTEGER,
    called      INTEGER,
    inferred_sex TEXT,
    notes       TEXT
);

CREATE TABLE IF NOT EXISTS genotype (
    kit_id INTEGER NOT NULL,
    chrom  TEXT NOT NULL,
    pos    INTEGER NOT NULL,
    rsid   TEXT,
    gt     TEXT NOT NULL,
    PRIMARY KEY (kit_id, chrom, pos)
) WITHOUT ROWID;

-- Phasing output for a kit whose parent is also in the project.  `mat` and
-- `pat` are single alleles; NULL means the site could not be resolved.
CREATE TABLE IF NOT EXISTS phased (
    kit_id INTEGER NOT NULL,
    chrom  TEXT NOT NULL,
    pos    INTEGER NOT NULL,
    mat    TEXT,
    pat    TEXT,
    method TEXT,
    PRIMARY KEY (kit_id, chrom, pos)
) WITHOUT ROWID;

-- IBD segments we computed ourselves from two kits' genotypes.
CREATE TABLE IF NOT EXISTS ibd (
    id       INTEGER PRIMARY KEY,
    kit_a    INTEGER NOT NULL,
    kit_b    INTEGER NOT NULL,
    chrom    TEXT NOT NULL,
    start_bp INTEGER NOT NULL,
    end_bp   INTEGER NOT NULL,
    cm       REAL,
    snps     INTEGER,
    kind     TEXT   -- 'HIR' (half-identical) or 'FIR' (fully identical)
);
CREATE INDEX IF NOT EXISTS ibd_pair ON ibd (kit_a, kit_b, chrom);

-- A DNA relative, as reported by a testing site, relative to one of our kits.
CREATE TABLE IF NOT EXISTS match (
    id            INTEGER PRIMARY KEY,
    kit_id        INTEGER NOT NULL REFERENCES kit(id),
    source        TEXT NOT NULL,
    remote_id     TEXT,
    name          TEXT,
    total_cm      REAL,
    seg_count     INTEGER,
    largest_cm    REAL,
    shared_x_cm   REAL,
    predicted     TEXT,    -- the vendor's own guess, kept for comparison
    side          TEXT,    -- 'maternal' | 'paternal' | 'both' | NULL
    side_evidence TEXT,
    tree_xref     TEXT,    -- link into the individual table, once identified
    sex           TEXT,
    birth_year    INTEGER,
    notes         TEXT,
    extra         TEXT,
    UNIQUE (kit_id, source, remote_id)
);
CREATE INDEX IF NOT EXISTS match_kit ON match (kit_id);
CREATE INDEX IF NOT EXISTS match_name ON match (name);

CREATE TABLE IF NOT EXISTS match_segment (
    id       INTEGER PRIMARY KEY,
    match_id INTEGER NOT NULL REFERENCES match(id) ON DELETE CASCADE,
    chrom    TEXT NOT NULL,
    start_bp INTEGER NOT NULL,
    end_bp   INTEGER NOT NULL,
    cm       REAL,
    snps     INTEGER
);
CREATE INDEX IF NOT EXISTS seg_match ON match_segment (match_id);
CREATE INDEX IF NOT EXISTS seg_locus ON match_segment (chrom, start_bp, end_bp);

-- "In common with" / shared-match edges, and true triangulation where the
-- platform actually reports B-C sharing rather than just co-membership.
CREATE TABLE IF NOT EXISTS shared_match (
    id       INTEGER PRIMARY KEY,
    kit_id   INTEGER NOT NULL,
    a_id     INTEGER NOT NULL REFERENCES match(id) ON DELETE CASCADE,
    b_id     INTEGER NOT NULL REFERENCES match(id) ON DELETE CASCADE,
    cm       REAL,
    kind     TEXT,   -- 'icw' | 'triangulated'
    UNIQUE (a_id, b_id, kind)
);

CREATE TABLE IF NOT EXISTS cluster (
    id      INTEGER PRIMARY KEY,
    kit_id  INTEGER NOT NULL,
    run     TEXT,
    name    TEXT,
    side    TEXT,
    ancestral_hint TEXT
);

CREATE TABLE IF NOT EXISTS cluster_member (
    cluster_id INTEGER NOT NULL REFERENCES cluster(id) ON DELETE CASCADE,
    match_id   INTEGER NOT NULL REFERENCES match(id) ON DELETE CASCADE,
    PRIMARY KEY (cluster_id, match_id)
);

-- Family tree, mirroring GEDCOM structure closely enough to round-trip.
CREATE TABLE IF NOT EXISTS individual (
    xref        TEXT PRIMARY KEY,
    given       TEXT,
    surname     TEXT,
    sex         TEXT,
    birth_date  TEXT,
    birth_place TEXT,
    death_date  TEXT,
    death_place TEXT,
    birth_year  INTEGER,
    death_year  INTEGER,
    notes       TEXT,
    source_file TEXT,
    extra       TEXT
);

CREATE TABLE IF NOT EXISTS family (
    xref       TEXT PRIMARY KEY,
    husb       TEXT,
    wife       TEXT,
    marr_date  TEXT,
    marr_place TEXT,
    source_file TEXT
);

CREATE TABLE IF NOT EXISTS child (
    fam_xref  TEXT NOT NULL REFERENCES family(xref) ON DELETE CASCADE,
    indi_xref TEXT NOT NULL REFERENCES individual(xref) ON DELETE CASCADE,
    rel       TEXT,  -- birth | adopted | foster, per GEDCOM PEDI
    PRIMARY KEY (fam_xref, indi_xref)
);
CREATE INDEX IF NOT EXISTS child_indi ON child (indi_xref);

-- Findings the hypothesis engine produced, so they survive between runs and
-- can be marked accepted or rejected by a human.
CREATE TABLE IF NOT EXISTS hypothesis (
    id        INTEGER PRIMARY KEY,
    kind      TEXT,
    subject   TEXT,
    summary   TEXT,
    score     REAL,
    status    TEXT DEFAULT 'open',
    detail    TEXT,
    created_at TEXT
);

CREATE TABLE IF NOT EXISTS sim_cache (
    key     TEXT PRIMARY KEY,
    payload TEXT
);
"""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass
class Kit:
    id: int
    label: str
    person_id: Optional[int]
    vendor: Optional[str]
    build: Optional[str]
    snp_count: Optional[int]
    inferred_sex: Optional[str]
    person: Optional[str] = None

    @property
    def display(self) -> str:
        who = self.person or "?"
        return f"{self.label} [{who}/{self.vendor or '?'}]"


class Store:
    def __init__(self, path: str):
        self.path = path
        first_time = not os.path.exists(path)
        parent = os.path.dirname(os.path.abspath(path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        self.db = sqlite3.connect(path)
        self.db.row_factory = sqlite3.Row
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.execute("PRAGMA synchronous=NORMAL")
        self.db.execute("PRAGMA foreign_keys=ON")
        self.db.execute("PRAGMA cache_size=-64000")  # ~64 MB page cache
        self.db.executescript(SCHEMA)
        if first_time:
            self.set_meta("schema_version", str(SCHEMA_VERSION))
            self.set_meta("created_at", _now())
        self.db.commit()

    # -- lifecycle ------------------------------------------------------

    def close(self) -> None:
        self.db.commit()
        self.db.close()

    def __enter__(self) -> "Store":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    # -- meta -----------------------------------------------------------

    def set_meta(self, key: str, value: str) -> None:
        self.db.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            (key, value),
        )
        self.db.commit()

    def get_meta(self, key: str, default: Optional[str] = None) -> Optional[str]:
        row = self.db.execute("SELECT value FROM meta WHERE key=?", (key,)).fetchone()
        return row["value"] if row else default

    # -- people and kits ------------------------------------------------

    def upsert_person(
        self,
        label: str,
        sex: Optional[str] = None,
        birth_year: Optional[int] = None,
        tree_xref: Optional[str] = None,
        notes: Optional[str] = None,
    ) -> int:
        row = self.db.execute("SELECT id FROM person WHERE label=?", (label,)).fetchone()
        if row:
            pid = row["id"]
            sets, vals = [], []
            for col, val in (
                ("sex", sex),
                ("birth_year", birth_year),
                ("tree_xref", tree_xref),
                ("notes", notes),
            ):
                if val is not None:
                    sets.append(f"{col}=?")
                    vals.append(val)
            if sets:
                vals.append(pid)
                self.db.execute(f"UPDATE person SET {', '.join(sets)} WHERE id=?", vals)
            self.db.commit()
            return pid
        cur = self.db.execute(
            "INSERT INTO person(label, sex, birth_year, tree_xref, notes) VALUES(?,?,?,?,?)",
            (label, sex, birth_year, tree_xref, notes),
        )
        self.db.commit()
        return int(cur.lastrowid)

    def person_id(self, label: str) -> Optional[int]:
        row = self.db.execute("SELECT id FROM person WHERE label=?", (label,)).fetchone()
        return row["id"] if row else None

    def people(self) -> List[sqlite3.Row]:
        return list(self.db.execute("SELECT * FROM person ORDER BY id"))

    def create_kit(
        self,
        label: str,
        person: str,
        vendor: str,
        build: str,
        source_path: str,
    ) -> int:
        pid = self.upsert_person(person)
        self.db.execute("DELETE FROM kit WHERE label=?", (label,))
        cur = self.db.execute(
            "INSERT INTO kit(label, person_id, vendor, build, source_path, imported_at) "
            "VALUES(?,?,?,?,?,?)",
            (label, pid, vendor, build, source_path, _now()),
        )
        self.db.commit()
        return int(cur.lastrowid)

    def kit(self, label_or_id: Any) -> Optional[Kit]:
        if isinstance(label_or_id, int) or str(label_or_id).isdigit():
            row = self.db.execute(
                "SELECT k.*, p.label AS person FROM kit k LEFT JOIN person p ON p.id=k.person_id "
                "WHERE k.id=?",
                (int(label_or_id),),
            ).fetchone()
        else:
            row = self.db.execute(
                "SELECT k.*, p.label AS person FROM kit k LEFT JOIN person p ON p.id=k.person_id "
                "WHERE k.label=?",
                (label_or_id,),
            ).fetchone()
        return self._kit_from_row(row) if row else None

    def kits(self) -> List[Kit]:
        rows = self.db.execute(
            "SELECT k.*, p.label AS person FROM kit k LEFT JOIN person p ON p.id=k.person_id "
            "ORDER BY k.id"
        )
        return [self._kit_from_row(r) for r in rows]

    def kits_for_person(self, person: str) -> List[Kit]:
        return [k for k in self.kits() if k.person == person]

    @staticmethod
    def _kit_from_row(row: sqlite3.Row) -> Kit:
        return Kit(
            id=row["id"],
            label=row["label"],
            person_id=row["person_id"],
            vendor=row["vendor"],
            build=row["build"],
            snp_count=row["snp_count"],
            inferred_sex=row["inferred_sex"],
            person=row["person"] if "person" in row.keys() else None,
        )

    def delete_kit(self, kit_id: int) -> None:
        self.db.execute("DELETE FROM genotype WHERE kit_id=?", (kit_id,))
        self.db.execute("DELETE FROM phased WHERE kit_id=?", (kit_id,))
        self.db.execute("DELETE FROM kit WHERE id=?", (kit_id,))
        self.db.commit()

    # -- genotypes ------------------------------------------------------

    def insert_genotypes(
        self, kit_id: int, rows: Iterable[Tuple[str, int, Optional[str], str]], batch: int = 50_000
    ) -> int:
        """Bulk-load ``(chrom, pos, rsid, gt)`` tuples.  Returns rows written."""
        total = 0
        buf: List[Tuple[int, str, int, Optional[str], str]] = []
        sql = (
            "INSERT OR REPLACE INTO genotype(kit_id, chrom, pos, rsid, gt) VALUES(?,?,?,?,?)"
        )
        for chrom, pos, rsid, gt in rows:
            buf.append((kit_id, chrom, pos, rsid, gt))
            if len(buf) >= batch:
                self.db.executemany(sql, buf)
                total += len(buf)
                buf.clear()
        if buf:
            self.db.executemany(sql, buf)
            total += len(buf)
        self.db.commit()
        return total

    def finalize_kit(
        self, kit_id: int, snp_count: int, called: int, inferred_sex: Optional[str]
    ) -> None:
        self.db.execute(
            "UPDATE kit SET snp_count=?, called=?, inferred_sex=? WHERE id=?",
            (snp_count, called, inferred_sex, kit_id),
        )
        self.db.commit()

    def genotypes(self, kit_id: int, chrom: Optional[str] = None) -> Iterator[sqlite3.Row]:
        if chrom:
            yield from self.db.execute(
                "SELECT chrom, pos, rsid, gt FROM genotype WHERE kit_id=? AND chrom=? ORDER BY pos",
                (kit_id, chrom),
            )
        else:
            from .genome import ALL_CHROMS

            for c in ALL_CHROMS:
                yield from self.db.execute(
                    "SELECT chrom, pos, rsid, gt FROM genotype "
                    "WHERE kit_id=? AND chrom=? ORDER BY pos",
                    (kit_id, c),
                )

    def paired_genotypes(
        self, kit_a: int, kit_b: int, chrom: str
    ) -> Iterator[Tuple[int, str, str]]:
        """Walk positions typed on both kits, in ascending order.

        Yields ``(pos, gt_a, gt_b)``.  This is the hot loop for IBD detection
        and phasing, so it is a single indexed join rather than two scans
        merged in Python.
        """
        cur = self.db.execute(
            "SELECT a.pos, a.gt, b.gt FROM genotype a "
            "JOIN genotype b ON b.kit_id=? AND b.chrom=a.chrom AND b.pos=a.pos "
            "WHERE a.kit_id=? AND a.chrom=? ORDER BY a.pos",
            (kit_b, kit_a, chrom),
        )
        for row in cur:
            yield row[0], row[1], row[2]

    def chroms_present(self, kit_id: int) -> List[str]:
        from .genome import ALL_CHROMS

        rows = self.db.execute(
            "SELECT DISTINCT chrom FROM genotype WHERE kit_id=?", (kit_id,)
        ).fetchall()
        present = {r[0] for r in rows}
        return [c for c in ALL_CHROMS if c in present]

    # -- phasing --------------------------------------------------------

    def clear_phase(self, kit_id: int) -> None:
        self.db.execute("DELETE FROM phased WHERE kit_id=?", (kit_id,))
        self.db.commit()

    def insert_phase(
        self,
        kit_id: int,
        rows: Iterable[Tuple[str, int, Optional[str], Optional[str], str]],
        batch: int = 50_000,
    ) -> int:
        total = 0
        buf: List[Tuple[Any, ...]] = []
        sql = (
            "INSERT OR REPLACE INTO phased(kit_id, chrom, pos, mat, pat, method) "
            "VALUES(?,?,?,?,?,?)"
        )
        for chrom, pos, mat, pat, method in rows:
            buf.append((kit_id, chrom, pos, mat, pat, method))
            if len(buf) >= batch:
                self.db.executemany(sql, buf)
                total += len(buf)
                buf.clear()
        if buf:
            self.db.executemany(sql, buf)
            total += len(buf)
        self.db.commit()
        return total

    # -- IBD ------------------------------------------------------------

    def clear_ibd(self, kit_a: int, kit_b: int) -> None:
        self.db.execute("DELETE FROM ibd WHERE kit_a=? AND kit_b=?", (kit_a, kit_b))
        self.db.commit()

    def insert_ibd(self, rows: Sequence[Tuple[Any, ...]]) -> None:
        self.db.executemany(
            "INSERT INTO ibd(kit_a, kit_b, chrom, start_bp, end_bp, cm, snps, kind) "
            "VALUES(?,?,?,?,?,?,?,?)",
            rows,
        )
        self.db.commit()

    def ibd_segments(self, kit_a: int, kit_b: int, kind: str = "HIR") -> List[sqlite3.Row]:
        return list(
            self.db.execute(
                "SELECT * FROM ibd WHERE kit_a=? AND kit_b=? AND kind=? ORDER BY chrom, start_bp",
                (kit_a, kit_b, kind),
            )
        )

    # -- matches --------------------------------------------------------

    def upsert_match(self, kit_id: int, source: str, remote_id: str, **fields: Any) -> int:
        cols = [
            "name", "total_cm", "seg_count", "largest_cm", "shared_x_cm", "predicted",
            "side", "side_evidence", "tree_xref", "sex", "birth_year", "notes", "extra",
        ]
        row = self.db.execute(
            "SELECT id FROM match WHERE kit_id=? AND source=? AND remote_id=?",
            (kit_id, source, remote_id),
        ).fetchone()
        payload = {c: fields.get(c) for c in cols if fields.get(c) is not None}
        if isinstance(payload.get("extra"), (dict, list)):
            payload["extra"] = json.dumps(payload["extra"])
        if row:
            mid = row["id"]
            if payload:
                sets = ", ".join(f"{c}=?" for c in payload)
                self.db.execute(
                    f"UPDATE match SET {sets} WHERE id=?", list(payload.values()) + [mid]
                )
            return mid
        keys = ["kit_id", "source", "remote_id"] + list(payload)
        vals = [kit_id, source, remote_id] + list(payload.values())
        cur = self.db.execute(
            f"INSERT INTO match({', '.join(keys)}) VALUES({', '.join('?' * len(keys))})", vals
        )
        return int(cur.lastrowid)

    def set_match_side(self, match_id: int, side: Optional[str], evidence: str) -> None:
        self.db.execute(
            "UPDATE match SET side=?, side_evidence=? WHERE id=?", (side, evidence, match_id)
        )

    def matches(self, kit_id: Optional[int] = None, source: Optional[str] = None) -> List[sqlite3.Row]:
        sql = "SELECT * FROM match"
        where, args = [], []
        if kit_id is not None:
            where.append("kit_id=?")
            args.append(kit_id)
        if source:
            where.append("source=?")
            args.append(source)
        if where:
            sql += " WHERE " + " AND ".join(where)
        sql += " ORDER BY total_cm DESC"
        return list(self.db.execute(sql, args))

    def match_by_id(self, match_id: int) -> Optional[sqlite3.Row]:
        return self.db.execute("SELECT * FROM match WHERE id=?", (match_id,)).fetchone()

    def replace_match_segments(self, match_id: int, segs: Sequence[Tuple[Any, ...]]) -> None:
        self.db.execute("DELETE FROM match_segment WHERE match_id=?", (match_id,))
        self.db.executemany(
            "INSERT INTO match_segment(match_id, chrom, start_bp, end_bp, cm, snps) "
            "VALUES(?,?,?,?,?,?)",
            [(match_id, *s) for s in segs],
        )

    def match_segments(self, match_id: Optional[int] = None) -> List[sqlite3.Row]:
        if match_id is None:
            return list(
                self.db.execute(
                    "SELECT s.*, m.name, m.kit_id FROM match_segment s "
                    "JOIN match m ON m.id=s.match_id ORDER BY s.chrom, s.start_bp"
                )
            )
        return list(
            self.db.execute(
                "SELECT * FROM match_segment WHERE match_id=? ORDER BY chrom, start_bp",
                (match_id,),
            )
        )

    def add_shared_match(
        self, kit_id: int, a_id: int, b_id: int,
        cm: Optional[float] = None, kind: str = "icw",
    ) -> None:
        lo, hi = sorted((a_id, b_id))
        self.db.execute(
            "INSERT OR IGNORE INTO shared_match(kit_id, a_id, b_id, cm, kind) VALUES(?,?,?,?,?)",
            (kit_id, lo, hi, cm, kind),
        )

    def shared_matches(self, kit_id: int, kind: Optional[str] = None) -> List[sqlite3.Row]:
        sql = "SELECT * FROM shared_match WHERE kit_id=?"
        args: List[Any] = [kit_id]
        if kind:
            sql += " AND kind=?"
            args.append(kind)
        return list(self.db.execute(sql, args))

    # -- tree -----------------------------------------------------------

    def upsert_individual(self, **fields: Any) -> None:
        cols = [
            "xref", "given", "surname", "sex", "birth_date", "birth_place",
            "death_date", "death_place", "birth_year", "death_year", "notes",
            "source_file", "extra",
        ]
        vals = [fields.get(c) for c in cols]
        self.db.execute(
            f"INSERT OR REPLACE INTO individual({', '.join(cols)}) "
            f"VALUES({', '.join('?' * len(cols))})",
            vals,
        )

    def upsert_family(self, **fields: Any) -> None:
        cols = ["xref", "husb", "wife", "marr_date", "marr_place", "source_file"]
        vals = [fields.get(c) for c in cols]
        self.db.execute(
            f"INSERT OR REPLACE INTO family({', '.join(cols)}) "
            f"VALUES({', '.join('?' * len(cols))})",
            vals,
        )

    def add_child(self, fam_xref: str, indi_xref: str, rel: Optional[str] = None) -> None:
        self.db.execute(
            "INSERT OR REPLACE INTO child(fam_xref, indi_xref, rel) VALUES(?,?,?)",
            (fam_xref, indi_xref, rel),
        )

    def individuals(self) -> List[sqlite3.Row]:
        return list(self.db.execute("SELECT * FROM individual ORDER BY xref"))

    def families(self) -> List[sqlite3.Row]:
        return list(self.db.execute("SELECT * FROM family ORDER BY xref"))

    def children(self) -> List[sqlite3.Row]:
        return list(self.db.execute("SELECT * FROM child"))

    def clear_tree(self) -> None:
        self.db.execute("DELETE FROM child")
        self.db.execute("DELETE FROM family")
        self.db.execute("DELETE FROM individual")
        self.db.commit()

    # -- hypotheses -----------------------------------------------------

    def add_hypothesis(
        self, kind: str, subject: str, summary: str, score: float, detail: Dict[str, Any]
    ) -> int:
        cur = self.db.execute(
            "INSERT INTO hypothesis(kind, subject, summary, score, detail, created_at) "
            "VALUES(?,?,?,?,?,?)",
            (kind, subject, summary, score, json.dumps(detail, default=str), _now()),
        )
        return int(cur.lastrowid)

    def clear_hypotheses(self, kind: Optional[str] = None) -> None:
        if kind:
            self.db.execute("DELETE FROM hypothesis WHERE kind=? AND status='open'", (kind,))
        else:
            self.db.execute("DELETE FROM hypothesis WHERE status='open'")
        self.db.commit()

    def hypotheses(self, kind: Optional[str] = None) -> List[sqlite3.Row]:
        if kind:
            return list(
                self.db.execute(
                    "SELECT * FROM hypothesis WHERE kind=? ORDER BY score DESC", (kind,)
                )
            )
        return list(self.db.execute("SELECT * FROM hypothesis ORDER BY score DESC"))

    # -- simulator cache ------------------------------------------------

    def sim_get(self, key: str) -> Optional[Any]:
        row = self.db.execute("SELECT payload FROM sim_cache WHERE key=?", (key,)).fetchone()
        return json.loads(row["payload"]) if row else None

    def sim_put(self, key: str, payload: Any) -> None:
        self.db.execute(
            "INSERT OR REPLACE INTO sim_cache(key, payload) VALUES(?,?)",
            (key, json.dumps(payload)),
        )
        self.db.commit()

    def commit(self) -> None:
        self.db.commit()

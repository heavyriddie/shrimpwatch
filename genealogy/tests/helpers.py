"""Shared fixtures: small synthetic kits built by simulating real inheritance.

The tests need ground truth, and the only way to have ground truth about a
parent-child pair is to construct one properly -- draw founder haplotypes,
recombine them, and transmit.  Genotypes faked without a meiosis model
produce Mendelian violations everywhere and would make the phasing tests
pass or fail for the wrong reasons.
"""

from __future__ import annotations

import os
import random
import tempfile
from typing import Dict, List, Tuple

from roots.genome import AUTOSOMES, CHROM_BP, CHROM_CM, GeneticMap
from roots.store import Store

ALLELES = ("A", "C", "G", "T")


def make_positions(per_chrom: int = 1200, seed: int = 1) -> Dict[str, List[int]]:
    rng = random.Random(seed)
    out: Dict[str, List[int]] = {}
    for chrom in AUTOSOMES:
        limit = CHROM_BP["37"][chrom]
        out[chrom] = sorted(rng.sample(range(1000, limit - 1000), per_chrom))
    return out


def make_founder(
    positions: Dict[str, List[int]], rng: random.Random
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Two haplotypes for one untransmitted individual."""
    h1: Dict[str, List[str]] = {}
    h2: Dict[str, List[str]] = {}
    for chrom, pos in positions.items():
        a1 = [rng.choice(ALLELES) for _ in pos]
        # The second allele agrees most of the time, which is what makes a
        # real genome mostly homozygous at any given site.
        a2 = [a if rng.random() < 0.7 else rng.choice(ALLELES) for a in a1]
        h1[chrom] = a1
        h2[chrom] = a2
    return h1, h2


def transmit(
    h1: Dict[str, List[str]],
    h2: Dict[str, List[str]],
    positions: Dict[str, List[int]],
    gmap: GeneticMap,
    rng: random.Random,
) -> Dict[str, List[str]]:
    """One gamete: a recombined mosaic of the two input haplotypes."""
    out: Dict[str, List[str]] = {}
    for chrom, pos in positions.items():
        length = CHROM_CM[chrom]
        breaks: List[float] = []
        x = rng.expovariate(0.01)
        while x < length:
            breaks.append(x)
            x += rng.expovariate(0.01)
        break_bp = [gmap.bp_at_cm(chrom, b) for b in breaks]
        current = rng.getrandbits(1)
        alleles: List[str] = []
        bi = 0
        for i, p in enumerate(pos):
            while bi < len(break_bp) and p >= break_bp[bi]:
                current ^= 1
                bi += 1
            alleles.append((h1 if current == 0 else h2)[chrom][i])
        out[chrom] = alleles
    return out


def genotypes(
    hap_a: Dict[str, List[str]],
    hap_b: Dict[str, List[str]],
    positions: Dict[str, List[int]],
    rng: random.Random,
    error_rate: float = 0.0,
    no_call_rate: float = 0.0,
):
    for chrom, pos in positions.items():
        for i, p in enumerate(pos):
            a, b = hap_a[chrom][i], hap_b[chrom][i]
            if no_call_rate and rng.random() < no_call_rate:
                yield chrom, p, f"rs{chrom}_{i}", "--"
                continue
            if error_rate and rng.random() < error_rate:
                a = rng.choice(ALLELES)
            yield chrom, p, f"rs{chrom}_{i}", "".join(sorted((a, b)))


def parent_child_store(
    per_chrom: int = 1200,
    seed: int = 7,
    error_rate: float = 0.002,
    no_call_rate: float = 0.003,
) -> Tuple[Store, int, int, GeneticMap, str]:
    """A store holding a true mother and her child.

    Returns ``(store, child_kit_id, mother_kit_id, gmap, db_path)``.
    """
    rng = random.Random(seed)
    gmap = GeneticMap.linear("37")
    positions = make_positions(per_chrom, seed)

    m1, m2 = make_founder(positions, rng)
    f1, f2 = make_founder(positions, rng)
    from_mother = transmit(m1, m2, positions, gmap, rng)
    from_father = transmit(f1, f2, positions, gmap, rng)

    fd, path = tempfile.mkstemp(suffix=".db")
    os.close(fd)
    os.unlink(path)
    store = Store(path)

    mother_kit = store.create_kit("mother", "mother", "test", "37", "synthetic")
    store.insert_genotypes(
        mother_kit, genotypes(m1, m2, positions, rng, error_rate, no_call_rate)
    )
    child_kit = store.create_kit("child", "self", "test", "37", "synthetic")
    store.insert_genotypes(
        child_kit,
        genotypes(from_mother, from_father, positions, rng, error_rate, no_call_rate),
    )
    return store, child_kit, mother_kit, gmap, path


def unrelated_store(per_chrom: int = 1200, seed: int = 11):
    """Two people with no relationship, as a negative control."""
    rng = random.Random(seed)
    gmap = GeneticMap.linear("37")
    positions = make_positions(per_chrom, seed)
    fd, path = tempfile.mkstemp(suffix=".db")
    os.close(fd)
    os.unlink(path)
    store = Store(path)
    ids = []
    for label in ("a", "b"):
        h1, h2 = make_founder(positions, rng)
        kit = store.create_kit(label, label, "test", "37", "synthetic")
        store.insert_genotypes(kit, genotypes(h1, h2, positions, rng))
        ids.append(kit)
    return store, ids[0], ids[1], gmap, path


def cleanup(store: Store, path: str) -> None:
    try:
        store.close()
    except Exception:
        pass
    for suffix in ("", "-wal", "-shm"):
        try:
            os.unlink(path + suffix)
        except OSError:
            pass

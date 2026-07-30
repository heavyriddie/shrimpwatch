"""Tests for match ingestion, side inference, clustering, and the tree."""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from roots.genome import GeneticMap
from roots.matches import cluster as cluster_mod
from roots.matches import ingest as ingest_mod
from roots.matches import sides as sides_mod
from roots.matches import triangulate as tri_mod
from roots.store import Store
from roots.tree import gedcom as gedcom_mod
from roots.tree import kinship
from tests import helpers

GEDCOM = """\
0 HEAD
1 SOUR test
1 GEDC
2 VERS 5.5.1
1 CHAR UTF-8
0 @I1@ INDI
1 NAME Alfred /Whitlock/
1 SEX M
1 BIRT
2 DATE ABT 1901
2 PLAC Leeds, England
1 FAMS @F1@
0 @I2@ INDI
1 NAME Edith /Ramsay/
1 SEX F
1 BIRT
2 DATE 1903
1 FAMS @F1@
0 @I3@ INDI
1 NAME Harold /Whitlock/
1 SEX M
1 BIRT
2 DATE 12 MAR 1930
1 FAMC @F1@
1 FAMS @F2@
0 @I4@ INDI
1 NAME Ruth /Whitlock/
1 SEX F
1 BIRT
2 DATE BET 1932 AND 1934
1 FAMC @F1@
1 FAMS @F3@
0 @I5@ INDI
1 NAME Diane /Whitlock/
1 SEX F
1 FAMC @F2@
0 @I6@ INDI
1 NAME Colin /Barrow/
1 SEX M
1 FAMC @F3@
0 @F1@ FAM
1 HUSB @I1@
1 WIFE @I2@
1 CHIL @I3@
1 CHIL @I4@
0 @F2@ FAM
1 HUSB @I3@
1 CHIL @I5@
0 @F3@ FAM
1 WIFE @I4@
1 CHIL @I6@
0 TRLR
"""

MATCHES_CSV = """\
Match Name,Shared cM,Shared Segments,Longest cM,Predicted Relationship
Mum Herself,3545,23,281,Parent
Alice Adams,1750,28,120,Close family
Bob Brown,880,32,88,1st cousin
Carol Clark,210,14,42,2nd cousin
Dan Davis,95,7,30,3rd cousin
Erin Ellis,45,4,22,4th cousin
"""

MOTHER_CSV = """\
Match Name,Shared cM,Shared Segments,Longest cM
Alice Adams,3540,23,281
Carol Clark,430,20,55
Erin Ellis,92,8,28
"""

SEGMENTS_CSV = """\
Match Name,Chromosome,Start Position,End Position,Centimorgans,Matching SNPs
Bob Brown,1,10000000,40000000,35.2,4200
Carol Clark,1,15000000,38000000,28.0,3600
Dan Davis,1,12000000,36000000,26.5,3300
Erin Ellis,7,5000000,25000000,21.0,2900
"""

ICW_CSV = """\
Match Name A,Match Name B
Bob Brown,Carol Clark
Bob Brown,Dan Davis
Carol Clark,Dan Davis
"""


def _temp_store():
    fd, path = tempfile.mkstemp(suffix=".db")
    os.close(fd)
    os.unlink(path)
    return Store(path), path


def _write(text: str, suffix: str = ".csv") -> str:
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    return path


class TestMatchIngest(unittest.TestCase):
    def setUp(self):
        self.store, self.path = _temp_store()
        self.kit = self.store.create_kit("me", "self", "test", "37", "-")

    def tearDown(self):
        helpers.cleanup(self.store, self.path)

    def test_columns_are_recognised_by_meaning(self):
        csv_path = _write(MATCHES_CSV)
        try:
            report = ingest_mod.import_matches(self.store, self.kit, csv_path, "Ancestry")
            self.assertEqual(report.imported, 6)
            self.assertEqual(report.columns["total_cm"], "Shared cM")
            self.assertEqual(report.columns["name"], "Match Name")
        finally:
            os.unlink(csv_path)

    def test_unusual_spellings_still_map(self):
        odd = "Display Name;Total cM;# Segments\nZoe Zane;300;12\n"
        csv_path = _write(odd)
        try:
            report = ingest_mod.import_matches(self.store, self.kit, csv_path)
            self.assertEqual(report.imported, 1)
            rows = self.store.matches(self.kit)
            self.assertEqual(rows[0]["total_cm"], 300)
        finally:
            os.unlink(csv_path)

    def test_percentage_conversion_is_flagged(self):
        pct = "Name,Shared DNA %\nPat Perry,25\n"
        csv_path = _write(pct)
        try:
            report = ingest_mod.import_matches(self.store, self.kit, csv_path)
            self.assertEqual(report.imported, 1)
            self.assertTrue(
                any("percentage" in w for w in report.warnings),
                "converting a percentage to cM must be flagged, since the "
                "percentage counts both chromosome copies",
            )
        finally:
            os.unlink(csv_path)

    def test_segments_join_by_name(self):
        m, s = _write(MATCHES_CSV), _write(SEGMENTS_CSV)
        try:
            ingest_mod.import_matches(self.store, self.kit, m, "Ancestry")
            report = ingest_mod.import_segments(self.store, self.kit, s, "Ancestry")
            self.assertEqual(report.imported, 4)
            self.assertEqual(report.unmatched_names, [])
            bob = [r for r in self.store.matches(self.kit) if r["name"] == "Bob Brown"][0]
            self.assertEqual(len(self.store.match_segments(bob["id"])), 1)
        finally:
            os.unlink(m)
            os.unlink(s)

    def test_name_normalisation_bridges_punctuation(self):
        self.assertEqual(
            ingest_mod.normalize_name("O'Brien,  Mary-Jane "),
            ingest_mod.normalize_name("OBrien Mary Jane"),
        )


class TestSideInference(unittest.TestCase):
    """The central inference: absence from a mother's list means paternal."""

    def setUp(self):
        self.store, self.path = _temp_store()
        self.child = self.store.create_kit("me", "self", "test", "37", "-")
        self.mother = self.store.create_kit("mum", "mother", "test", "37", "-")
        m, mm = _write(MATCHES_CSV), _write(MOTHER_CSV)
        try:
            ingest_mod.import_matches(self.store, self.child, m, "Ancestry")
            ingest_mod.import_matches(self.store, self.mother, mm, "Ancestry")
        finally:
            os.unlink(m)
            os.unlink(mm)

    def tearDown(self):
        helpers.cleanup(self.store, self.path)

    def _sides(self, **kw):
        report = sides_mod.assign_from_parent_list(
            self.store, self.child, self.mother, "mother",
            floor_multiple=0.0, min_paternal_cm=20.0, **kw
        )
        return {c.name: c.side for c in report.calls}, report

    def test_matches_on_the_mothers_list_are_maternal(self):
        sides, _ = self._sides()
        self.assertEqual(sides["Carol Clark"], "maternal")
        self.assertEqual(sides["Erin Ellis"], "maternal")

    def test_absence_from_the_mothers_list_means_paternal(self):
        sides, _ = self._sides()
        self.assertEqual(sides["Bob Brown"], "paternal")
        self.assertEqual(sides["Dan Davis"], "paternal")

    def test_the_mother_herself_is_not_assigned_a_side(self):
        sides, _ = self._sides()
        self.assertIsNone(
            sides["Mum Herself"],
            "someone sharing a whole genome and missing from the mother's own "
            "list is the mother; she has no side",
        )

    def test_the_mothers_other_child_is_a_half_sibling(self):
        # Alice shares her whole genome with the mother, so she is the
        # mother's child; sharing only 1750 cM with us makes her a half
        # sibling rather than a full one.
        sides, _ = self._sides()
        self.assertEqual(sides["Alice Adams"], "maternal")

    def test_cross_platform_absence_proves_nothing(self):
        other = self.store.create_kit("mum2", "mother2", "test", "37", "-")
        mm = _write(MOTHER_CSV)
        try:
            ingest_mod.import_matches(self.store, other, mm, "23andMe")
        finally:
            os.unlink(mm)
        report = sides_mod.assign_from_parent_list(
            self.store, self.child, other, "mother", floor_multiple=0.0
        )
        self.assertEqual(report.calls, [])
        self.assertTrue(any("platform" in w for w in report.warnings))

    def test_truncated_parent_list_suppresses_inference(self):
        # With the default safety margin the mother's list, which stops at
        # 92 cM, cannot justify calling a 95 cM match paternal.
        report = sides_mod.assign_from_parent_list(
            self.store, self.child, self.mother, "mother"
        )
        by_name = {c.name: c for c in report.calls}
        self.assertIsNone(by_name["Dan Davis"].side)
        self.assertEqual(by_name["Bob Brown"].side, "paternal")

    def test_full_sibling_is_related_on_both_sides(self):
        store, path = _temp_store()
        try:
            child = store.create_kit("me", "self", "test", "37", "-")
            mother = store.create_kit("mum", "mother", "test", "37", "-")
            store.upsert_match(child, "X", "sib", name="Sam Sibling", total_cm=2600)
            store.upsert_match(mother, "X", "sib", name="Sam Sibling", total_cm=3540)
            store.commit()
            report = sides_mod.assign_from_parent_list(store, child, mother, "mother")
            self.assertEqual(report.calls[0].side, "both")
        finally:
            helpers.cleanup(store, path)

    def test_half_sibling_is_maternal_only(self):
        store, path = _temp_store()
        try:
            child = store.create_kit("me", "self", "test", "37", "-")
            mother = store.create_kit("mum", "mother", "test", "37", "-")
            store.upsert_match(child, "X", "hs", name="Hal Half", total_cm=1750)
            store.upsert_match(mother, "X", "hs", name="Hal Half", total_cm=3540)
            store.commit()
            report = sides_mod.assign_from_parent_list(store, child, mother, "mother")
            self.assertEqual(report.calls[0].side, "maternal")
        finally:
            helpers.cleanup(store, path)


class TestClustering(unittest.TestCase):
    def setUp(self):
        self.store, self.path = _temp_store()
        self.kit = self.store.create_kit("me", "self", "test", "37", "-")

    def tearDown(self):
        helpers.cleanup(self.store, self.path)

    def test_two_disjoint_groups_become_two_clusters(self):
        group_a = [self.store.upsert_match(self.kit, "X", f"a{i}", name=f"A{i}",
                                           total_cm=100 + i) for i in range(5)]
        group_b = [self.store.upsert_match(self.kit, "X", f"b{i}", name=f"B{i}",
                                           total_cm=100 + i) for i in range(5)]
        for group in (group_a, group_b):
            for i, x in enumerate(group):
                for y in group[i + 1:]:
                    self.store.add_shared_match(self.kit, x, y)
        self.store.commit()
        result = cluster_mod.cluster_matches(self.store, self.kit, min_cm=40, max_cm=400)
        self.assertEqual(len(result.clusters), 2)
        for c in result.clusters:
            self.assertEqual(c.size, 5)
            self.assertEqual(c.cohesion, 1.0)

    def test_cluster_takes_the_side_of_its_members(self):
        ids = [self.store.upsert_match(self.kit, "X", f"p{i}", name=f"P{i}",
                                       total_cm=120, side="paternal") for i in range(4)]
        for i, x in enumerate(ids):
            for y in ids[i + 1:]:
                self.store.add_shared_match(self.kit, x, y)
        self.store.commit()
        result = cluster_mod.cluster_matches(self.store, self.kit, min_cm=40, max_cm=400)
        self.assertEqual(result.clusters[0].side, "paternal")

    def test_clustering_needs_shared_match_data(self):
        for i in range(5):
            self.store.upsert_match(self.kit, "X", f"z{i}", name=f"Z{i}", total_cm=150)
        self.store.commit()
        result = cluster_mod.cluster_matches(self.store, self.kit)
        self.assertEqual(result.clusters, [])
        self.assertTrue(any("shared-match" in w for w in result.warnings))


class TestTriangulation(unittest.TestCase):
    def setUp(self):
        self.store, self.path = _temp_store()
        self.kit = self.store.create_kit("me", "self", "test", "37", "-")
        self.gmap = GeneticMap.linear("37")

    def tearDown(self):
        helpers.cleanup(self.store, self.path)

    def _add(self, name, side, chrom, start, end, cm=250.0):
        mid = self.store.upsert_match(self.kit, "X", name, name=name,
                                      total_cm=cm, side=side)
        self.store.replace_match_segments(mid, [(chrom, start, end, None, 3000)])
        self.store.commit()
        return mid

    def test_overlapping_matches_form_a_group(self):
        self._add("A", "maternal", "1", 10_000_000, 50_000_000)
        self._add("B", "maternal", "1", 20_000_000, 60_000_000)
        groups = tri_mod.overlap_groups(self.store, self.kit, self.gmap)
        self.assertTrue(groups)
        self.assertEqual(groups[0].size, 2)
        self.assertEqual(groups[0].side(), "maternal")
        self.assertFalse(groups[0].conflict)

    def test_opposite_sides_overlapping_is_flagged_not_merged(self):
        self._add("A", "maternal", "1", 10_000_000, 50_000_000)
        self._add("B", "paternal", "1", 20_000_000, 60_000_000)
        groups = tri_mod.overlap_groups(self.store, self.kit, self.gmap)
        self.assertTrue(groups[0].conflict)
        self.assertEqual(groups[0].side(), "conflicting")

    def test_close_relatives_are_excluded(self):
        self._add("Parent", "maternal", "1", 1, 240_000_000, cm=3540)
        self._add("Cousin", "maternal", "1", 20_000_000, 60_000_000, cm=800)
        groups = tri_mod.overlap_groups(self.store, self.kit, self.gmap)
        self.assertFalse(
            any("Parent" in g.names for g in groups),
            "a parent overlaps every segment and locates nothing",
        )

    def test_coverage_counts_each_region_once(self):
        self._add("A", "maternal", "1", 10_000_000, 50_000_000)
        self._add("B", "maternal", "1", 10_000_000, 50_000_000)
        painted = tri_mod.paint(self.store, self.kit, self.gmap)
        cov = tri_mod.coverage(painted, self.gmap)
        single = self.gmap.length_cm("1", 10_000_000, 50_000_000)
        self.assertAlmostEqual(cov["maternal"], single, places=3)
        self.assertGreater(cov["unattributed"], 3000)


class TestGedcomAndKinship(unittest.TestCase):
    def setUp(self):
        self.store, self.path = _temp_store()
        self.ged = _write(GEDCOM, ".ged")
        gedcom_mod.import_gedcom(self.store, self.ged)
        self.idx = kinship.load_tree(self.store)

    def tearDown(self):
        helpers.cleanup(self.store, self.path)
        os.unlink(self.ged)

    def test_import_counts(self):
        self.assertEqual(len(self.idx.individuals), 6)
        self.assertEqual(len(self.idx.families), 3)

    def test_approximate_dates_yield_years(self):
        alfred = [i for i in self.idx.individuals.values()
                  if i.get("surname") == "Whitlock" and i.get("given") == "Alfred"][0]
        self.assertEqual(alfred["birth_year"], 1901)

    def test_between_dates_take_the_first_year(self):
        ruth = [i for i in self.idx.individuals.values()
                if i.get("given") == "Ruth"][0]
        self.assertEqual(ruth["birth_year"], 1932)

    def test_siblings(self):
        harold = self._xref("Harold")
        ruth = self._xref("Ruth")
        paths = kinship.relationship_between(self.idx, harold, ruth)
        self.assertEqual(paths[0].relationship.name, "full sibling")
        self.assertEqual(paths[0].ancestors, 2)

    def test_first_cousins(self):
        paths = kinship.relationship_between(self.idx, self._xref("Diane"),
                                             self._xref("Colin"))
        self.assertEqual(paths[0].relationship.name, "first cousin")

    def test_direct_line(self):
        paths = kinship.relationship_between(self.idx, self._xref("Alfred"),
                                             self._xref("Diane"))
        self.assertEqual(paths[0].relationship.name, "grandparent/grandchild")
        self.assertEqual(paths[0].up, 0)

    def test_unrelated_people_have_no_path(self):
        edith = self._xref("Edith")
        # Edith's husband is not her relative.
        self.assertEqual(
            kinship.relationship_between(self.idx, edith, self._xref("Alfred")), []
        )

    def test_ancestors_include_self_at_depth_zero(self):
        diane = self._xref("Diane")
        depths = kinship.ancestors(self.idx, diane)
        self.assertEqual(depths[diane], 0)
        self.assertEqual(depths[self._xref("Alfred")], 2)

    def test_tree_gaps_report_incompleteness(self):
        gaps = kinship.tree_gaps(self.idx, self._xref("Diane"), max_generations=3)
        self.assertEqual(gaps[1][1], 2)
        self.assertLess(gaps[1][0], 2, "Diane's mother is not in the file")

    def test_export_round_trips(self):
        fd, out = tempfile.mkstemp(suffix=".ged")
        os.close(fd)
        try:
            gedcom_mod.export_gedcom(self.store, out)
            parsed = gedcom_mod.parse_gedcom(out)
            self.assertEqual(len(parsed.individuals), 6)
            self.assertEqual(len(parsed.families), 3)
        finally:
            os.unlink(out)

    def _xref(self, given: str) -> str:
        for xref, ind in self.idx.individuals.items():
            if (ind.get("given") or "").startswith(given):
                return xref
        raise AssertionError(f"{given} not found")


class TestAdoptionIsNotGenetic(unittest.TestCase):
    def test_adoptive_links_are_excluded_from_genetic_paths(self):
        store, path = _temp_store()
        try:
            store.upsert_individual(xref="P", given="Pat", surname="Parent", sex="F")
            store.upsert_individual(xref="C", given="Chris", surname="Child", sex="M")
            store.upsert_family(xref="F", wife="P")
            store.add_child("F", "C", "adopted")
            store.commit()
            idx = kinship.load_tree(store)
            genetic = kinship.ancestors(idx, "C", genetic_only=True)
            documented = kinship.ancestors(idx, "C", genetic_only=False)
            self.assertNotIn("P", genetic)
            self.assertIn("P", documented)
        finally:
            helpers.cleanup(store, path)


if __name__ == "__main__":
    unittest.main(verbosity=2)

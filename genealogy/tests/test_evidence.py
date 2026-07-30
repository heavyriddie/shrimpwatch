"""Tests for the documentary evidence layer and the research planner."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from roots import evidence as ev
from roots.store import Store
from roots.tree import gedcom as gedcom_mod
from roots.tree import kinship
from tests import helpers

SOURCED_GEDCOM = """\
0 HEAD
1 GEDC
2 VERS 5.5.1
1 CHAR UTF-8
0 @S1@ SOUR
1 TITL 1911 Census of England and Wales
1 AUTH The National Archives
0 @S2@ SOUR
1 TITL GRO Marriage Certificate, Whitlock-Ramsay
0 @I1@ INDI
1 NAME Alfred /Whitlock/
1 SEX M
1 BIRT
2 DATE 1901
2 PLAC Leeds, Yorkshire, England
2 SOUR @S1@
1 SOUR @S2@
0 TRLR
"""


def _temp_store():
    fd, path = tempfile.mkstemp(suffix=".db")
    os.close(fd)
    os.unlink(path)
    return Store(path), path


def _write(text: str, suffix: str = ".ged") -> str:
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    return path


class TestJurisdiction(unittest.TestCase):
    def test_places_map_to_record_jurisdictions(self):
        self.assertEqual(ev.infer_country("Glasgow, Lanarkshire"), "scotland")
        self.assertEqual(ev.infer_country("Cork, Ireland"), "ireland")
        self.assertEqual(ev.infer_country("Belfast, Co. Antrim"), "northern ireland")
        self.assertEqual(ev.infer_country("Leeds, Yorkshire, England"), "england")
        self.assertEqual(ev.infer_country("Cardiff, Glamorgan"), "wales")

    def test_unknown_places_are_not_guessed(self):
        self.assertIsNone(ev.infer_country("Somewhere"))
        self.assertIsNone(ev.infer_country(None))


class TestResearchSuggestions(unittest.TestCase):
    def setUp(self):
        self.store, self.path = _temp_store()

    def tearDown(self):
        helpers.cleanup(self.store, self.path)

    def _person(self, xref, year, place, sex="M"):
        self.store.upsert_individual(
            xref=xref, given="Test", surname="Person", sex=sex,
            birth_year=year, birth_place=place,
        )
        self.store.commit()
        return kinship.load_tree(self.store)

    def test_english_wall_leads_with_the_marriage_certificate(self):
        idx = self._person("A", 1880, "Leeds, England")
        tasks = ev.suggest_for_person(idx, "A", depth=2)
        self.assertIn("marriage", tasks[0].record_set.lower())
        self.assertIn(
            "fathers of", tasks[0].rationale,
            "the reason a marriage certificate comes first is that it names "
            "both fathers, and the output should say so",
        )

    def test_scottish_wall_leads_with_the_death_register(self):
        idx = self._person("A", 1880, "Aberdeen, Scotland")
        tasks = ev.suggest_for_person(idx, "A", depth=2)
        self.assertIn("death", tasks[0].record_set.lower())
        self.assertEqual(tasks[0].repository, "ScotlandsPeople")

    def test_pre_civil_registration_goes_to_parish_registers(self):
        idx = self._person("A", 1800, "Devon, England")
        tasks = ev.suggest_for_person(idx, "A", depth=3)
        sets = " ".join(t.record_set.lower() for t in tasks)
        self.assertIn("parish", sets)
        self.assertNotIn(
            "gro birth index", sets,
            "civil registration did not exist in 1800, so suggesting the GRO "
            "index would send the user somewhere with no record to find",
        )

    def test_irish_research_offers_the_free_censuses(self):
        idx = self._person("A", 1870, "Galway, Ireland")
        tasks = ev.suggest_for_person(idx, "A", depth=2)
        free = [t for t in tasks if t.cost_band == "free"]
        self.assertTrue(free, "the surviving Irish censuses and BMD index are free")

    def test_closer_walls_outrank_distant_ones(self):
        idx = self._person("A", 1880, "Leeds, England")
        near = ev.suggest_for_person(idx, "A", depth=2)[0]
        far = ev.suggest_for_person(idx, "A", depth=8)[0]
        self.assertGreater(near.priority, far.priority)

    def test_paternal_walls_are_raised(self):
        idx = self._person("A", 1880, "Leeds, England")
        plain = ev.suggest_for_person(idx, "A", depth=3)[0]
        paternal = ev.suggest_for_person(idx, "A", depth=3, side="paternal")[0]
        self.assertGreater(paternal.priority, plain.priority)


class TestFrontier(unittest.TestCase):
    def setUp(self):
        self.store, self.path = _temp_store()
        self.store.upsert_individual(xref="ME", given="Me", surname="X", sex="F")
        self.store.upsert_individual(xref="MUM", given="Mum", surname="X", sex="F")
        self.store.upsert_individual(xref="GRAN", given="Gran", surname="Y", sex="F")
        self.store.upsert_family(xref="F1", wife="MUM")
        self.store.add_child("F1", "ME")
        self.store.upsert_family(xref="F2", wife="GRAN")
        self.store.add_child("F2", "MUM")
        self.store.commit()
        self.idx = kinship.load_tree(self.store)

    def tearDown(self):
        helpers.cleanup(self.store, self.path)

    def test_frontier_is_ancestors_without_recorded_parents(self):
        walls = dict(ev.frontier(self.idx, "ME"))
        self.assertIn("GRAN", walls)
        self.assertEqual(walls["GRAN"], 2)
        self.assertNotIn("MUM", walls, "her mother is recorded, so she is not a wall")

    def test_tasks_are_generated_and_persisted(self):
        tasks = ev.generate_tasks(self.store, root="ME", idx=self.idx)
        self.assertTrue(tasks)
        rows = self.store.db.execute("SELECT COUNT(*) FROM research_task").fetchone()
        self.assertEqual(rows[0], len(tasks))

    def test_no_root_and_no_kit_produces_nothing(self):
        self.assertEqual(ev.generate_tasks(self.store, root=None, idx=self.idx), [])


class TestSourcesAndEvidence(unittest.TestCase):
    def setUp(self):
        self.store, self.path = _temp_store()
        self.store.upsert_individual(xref="A", given="Alfred", surname="Whitlock")
        self.store.upsert_individual(xref="B", given="Edith", surname="Ramsay")
        self.store.commit()

    def tearDown(self):
        helpers.cleanup(self.store, self.path)

    def test_recording_a_source_and_a_claim(self):
        sid = ev.add_source(self.store, "GRO birth certificate", repository="GRO",
                            record_type="birth", cost=3.0)
        ev.add_evidence(self.store, sid, "A", "parentage",
                        "father recorded as Thomas Whitlock, labourer")
        rows = ev.evidence_for(self.store, "A")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["claim_type"], "parentage")
        self.assertEqual(rows[0]["repository"], "GRO")

    def test_contradicting_documents_are_kept_not_resolved(self):
        s1 = ev.add_source(self.store, "Census 1911")
        s2 = ev.add_source(self.store, "Birth certificate")
        ev.add_evidence(self.store, s1, "A", "birth", "born about 1901")
        ev.add_evidence(self.store, s2, "A", "birth", "born 1899", supports=False)
        conflicts = ev.contradictions(self.store)
        self.assertEqual(len(conflicts), 1)
        self.assertEqual(len(ev.evidence_for(self.store, "A")), 2,
                         "both documents stay on file; DNA adjudicates")

    def test_coverage_reports_the_unsourced_majority(self):
        sid = ev.add_source(self.store, "Census 1911", cost=0.0)
        ev.add_evidence(self.store, sid, "A", "residence", "Leeds")
        cov = ev.coverage(self.store)
        self.assertEqual(cov.individuals, 2)
        self.assertEqual(cov.with_evidence, 1)
        self.assertAlmostEqual(cov.fraction, 0.5)

    def test_costs_accumulate(self):
        ev.add_source(self.store, "cert one", cost=3.0)
        ev.add_source(self.store, "cert two", cost=12.5)
        self.assertAlmostEqual(ev.coverage(self.store).total_cost, 15.5)


class TestGedcomCitations(unittest.TestCase):
    def test_citations_import_from_a_gedcom(self):
        store, path = _temp_store()
        ged = _write(SOURCED_GEDCOM)
        try:
            gedcom_mod.import_gedcom(store, ged)
            stats = ev.import_gedcom_sources(store, ged)
            self.assertEqual(stats["sources"], 2)
            self.assertGreaterEqual(stats["citations"], 2)
            titles = {r["title"] for r in ev.sources(store)}
            self.assertIn("1911 Census of England and Wales", titles)
            kinds = {r["record_type"] for r in ev.sources(store)}
            self.assertIn("census", kinds)
            claims = {r["claim_type"] for r in ev.evidence_for(store, "I1")}
            self.assertIn("birth", claims,
                          "a SOUR under BIRT is evidence about the birth")
        finally:
            helpers.cleanup(store, path)
            os.unlink(ged)


if __name__ == "__main__":
    unittest.main(verbosity=2)

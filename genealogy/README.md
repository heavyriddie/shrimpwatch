# roots — genetic genealogy from your DNA and your mother's

An offline toolkit for building and investigating a deep family tree from
consumer DNA tests, match lists, and GEDCOM trees. It makes no network
requests. Your genotypes, your relatives' names, and everything derived from
them stay in one SQLite file on your machine.

It is built around one specific and unusually strong starting position:
**you have your own DNA and your mother's.** That combination is worth far
more than two separate tests, because it lets you split your genome — and
therefore every DNA match you will ever have — into a maternal half and a
paternal half, without your father ever being tested.

```
roots demo --out demo-data          # a synthetic family, to try it on
roots guide                         # the recommended order of operations
```

---

## What it does

**Reads raw DNA from any of the major vendors.** 23andMe, AncestryDNA,
MyHeritage, FamilyTreeDNA, LivingDNA — plain, gzipped or zipped. Chromosome
codes, no-call spellings and hemizygous conventions all differ between them
and are normalised on the way in.

**Merges a person's kits across vendors.** Two vendors type overlapping but
different SNP sets; combining a 23andMe and an Ancestry download of the same
person typically yields 25–35% more usable positions than either alone, which
sharpens every segment boundary downstream. Kits are checked for actually
being the same person before they are merged, and sites where the two vendors
disagree become no-calls rather than a coin flip.

**Finds shared DNA segments from genotypes directly.** No upload, no third
party. Give it two kits and it reports half-identical and fully-identical
regions and names the relationship for anything out to full siblings. On the
bundled synthetic family it recovers 3,543 cM of a true 3,544.6 cM
mother–child relationship and identifies the pair correctly.

**Phases you against your mother.** At every site where you are heterozygous
and she is homozygous, her contribution is fixed and the other allele must be
your father's. Genome-wide, this separates your two chromosome copies and
reconstructs half of an untested man. The paternal haplotype exports as an
uploadable pseudo-kit: every match it returns is a paternal-side relative, by
construction.

**Sorts your match list into maternal and paternal.** Anyone on your mother's
match list is maternal. Anyone absent from it who shares enough DNA that they
*would* have appeared is paternal. The inference refuses to run across
different testing platforms, and refuses to read anything into an absence
below the point where your mother's list is truncated — the two mistakes that
make this technique go wrong in practice.

**Clusters matches into ancestral lines** from shared-match data alone, so it
works on Ancestry results where no segment data exists. Average-linkage
agglomerative clustering on Jaccard distance, via the nearest-neighbour chain
algorithm, so a few thousand matches cluster in seconds.

**Triangulates** where segment data does exist, and is careful about the
distinction between matches that merely overlap the same coordinates and
matches that genuinely descend from the same ancestor. You have two copies of
every chromosome; two matches on opposite copies overlap perfectly and are
completely unrelated. Groups are labelled by which standard they actually
meet, and side conflicts are flagged rather than silently merged.

**Predicts relationships by simulation, not table lookup.** See below.

**Reads and writes GEDCOM 5.5.1**, tolerating the real-world messiness of
exports — CONC/CONT folding, `ABT 1832` and `BET 1830 AND 1840` dates, missing
FAMC backlinks, ANSEL declarations, xref collisions between files. Adoptive
and foster links are tracked but excluded from genetic reckoning, because an
adoptive parent transmits no DNA.

**Tracks documentary evidence.** Sources, citations, and what each record
actually asserts about a person — kept separate from the tree itself, so two
documents can disagree and both stay on file. It also generates a prioritised
list of records to look up next, derived from where your tree runs out and
where the DNA says the answers must be. The suggestions are
jurisdiction-aware: an England and Wales marriage certificate names both
fathers, so it advances two lines for one purchase and is usually the right
thing to buy first, whereas in Scotland a death certificate names both parents
of the deceased and is cheaper still.

**Fits a research plan to a budget.** Free lookups are never deferred,
subscription costs are treated as shared across every task they cover rather
than charged per record, and what remains is chosen by yield per pound. One
command, `roots auto`, runs every free and deterministic step end to end and
is safe to put on a cron job — it makes no network requests, so it costs
nothing to run forever.

**Investigates.** Checks documented relationships against measured sharing and
flags the ones that cannot both be true; projects your paternal matches onto
your untested father so they read as *his* relatives; identifies which clusters
connect to your tree and which are branches you have no documentation for; and
reports how much of your genome no match explains.

---

## Relationship prediction

Ask what a given amount of shared DNA means and you get a ranked list with
probabilities:

```
$ roots predict 212

relationship                    probability  avg cM  typical range
second cousin                   32.2%        200     80-345 cM
second cousin once removed      28.5%        101     24-209 cM
first cousin twice removed      19.0%        199     80-339 cM
half first cousin once removed   6.2%        206     70-367 cM
```

Those numbers are not looked up. The toolkit builds the actual pedigree,
gives each founder two uniquely labelled chromosome copies, simulates
recombination down every line as a Poisson process at one crossover per
Morgan, and measures what the two target people ended up sharing. Run a few
thousand times, that gives a real distribution — including the probability of
sharing nothing at all, which is the number people forget and which dominates
the interpretation of distant relationships.

Simulating the whole pedigree rather than using a closed form matters. The
common shortcut of treating two shared ancestors as independent contributors
predicts about 1525 cM for an aunt and her nephew; simulating the pedigree
gives 1742 cM, which is what aunt–nephew pairs actually show. The difference
is that both of the nephew's shared grandparents reach him through one
parent, so at any position his chromosome traces to exactly one of them.

Calibration against published averages, at a 7 cM threshold:

| relationship | simulated | published |
|---|---|---|
| parent/child | 3545 | 3485 (whole genome) |
| full siblings | 2621 | 2613 |
| half sib / grandparent / aunt-uncle | 1746 / 1766 / 1748 | ~1750 |
| first cousins | 854 | 866 |
| second cousins | 208 | 229 |

Close relationships land essentially on the published values. Cousin levels
come in slightly low and the gap widens with distance — which is expected,
since crowd-sourced averages are collected from pairs who *matched*, and pairs
who share nothing never get submitted. The simulator counts the zeros.

Three things feed the ranking: the simulated likelihood, a prior from how many
relatives of each type a person actually has (you have roughly 7 first cousins
and roughly 900 fourth cousins, which is why a 40 cM match is usually distant
rather than an unusually low close one), and a penalty for relationships
requiring an implausible age difference. Supply a real age gap with
`--age-gap` and it replaces the stand-in. `--flat-prior` shows raw likelihoods.

---

## Getting started

```bash
cd genealogy
python3 -m roots demo --out demo-data     # synthetic family with known truth
python3 -m roots guide                    # the full workflow, step by step
```

Install as a command with `pip install -e .`, which puts `roots` on your path.
Python 3.10+, standard library only, no dependencies.

The short version of the workflow:

```bash
roots init
roots import-dna you-23andme.txt  --person self   --label self-23
roots import-dna you-ancestry.txt --person self   --label self-anc
roots import-dna mum-23andme.txt  --person mother --label mum-23
roots merge-kits --person self
roots merge-kits --person mother
roots compare self-merged mother-merged                  # confirms she is your mother
roots phase --child self-merged --parent mother-merged \
            --export-other father-inferred.txt           # extracts your father
roots import-matches you-matches.csv --kit self-merged
roots import-matches mum-matches.csv --kit mother-merged
roots sides --child self-merged --parent mother-merged   # splits the match list
roots cluster --kit self-merged
roots import-tree yourtree.ged
roots import-citations yourtree.ged                      # what you already sourced
roots investigate --kit self-merged
roots research --kit self-merged --why                   # what to look up next
roots report --kit self-merged --out report.html
```

**[docs/DATA_SOURCES.md](docs/DATA_SOURCES.md)** covers the DNA side: what to
download from which testing site, which of them can triangulate at all, and
how to record relationships you already know.

**[docs/RECORD_SOURCES.md](docs/RECORD_SOURCES.md)** covers the documentary
side: births, marriages, deaths, censuses, parish registers and wills, what
each provider costs, and what order to spend money in.

**[docs/AUTOMATION.md](docs/AUTOMATION.md)** covers what can and cannot be
automated — which services have usable APIs, which forbid it, and what a
strict budget actually buys.

---

## Accuracy and honesty about limits

**The built-in genetic map is an approximation.** Recombination rates vary by
an order of magnitude along a chromosome, and the fallback map assumes a
uniform rate. It is fine for sanity checks and it is all the simulator needs
(which only uses total genetic lengths), but centimorgan figures computed for
specific segments will be off. Load a real recombination map with `--map` and
the warning goes away. Where a testing company has reported a segment's cM,
that value is used rather than recomputed.

**Small segments are mostly not real.** Below about 7 cM the majority of
reported segments in consumer data are not recent inheritance, and no amount
of triangulation makes them so. Thresholds default to 7 cM and 500 SNPs.

**Sharing at 1300–2300 cM is genuinely ambiguous.** Half siblings,
grandparents and aunts all sit near 1750 cM and total sharing cannot separate
them. Ages, the X chromosome, and which of them appears on your mother's list
can.

**A missing relative proves very little.** Around 10% of true third cousins
and a much larger share of fourth cousins share no detectable DNA at all.
`roots simulate` prints that probability for every relationship.

**Endogamy breaks the arithmetic.** If your ancestors come from a small or
intermarrying population, everyone is related through many lines at once,
totals run far above what a single relationship predicts, and clusters merge
into each other. The toolkit flags the symptoms — matches related through
multiple tree paths, excess sharing, low cluster cohesion — but cannot correct
for it.

---

## Privacy

Your DNA is not only yours. It discloses information about your mother, your
siblings, your children and relatives you have never met, none of whom
consented to anything you do with it. That is worth thinking about before
uploading a kit anywhere, and it applies with particular force to your
mother's kit, which is her data.

This toolkit makes no network requests at all. Generated HTML reports embed no
external resources and load offline forever. The project file and the reports
both contain identifiable genetic information, so treat them the way you would
treat medical records.

---

## Layout

```
roots/
  genome.py            reference geometry, genetic maps, build detection
  store.py             the SQLite project file
  dna/
    raw.py             vendor raw-file parsing and normalisation
    qc.py              quality control, concordance, cross-vendor merging
    relate.py          IBD segment detection from genotypes
    phase.py           parent-based phasing, untested-parent extraction
    simulate.py        pedigree meiosis simulator
    predict.py         Bayesian relationship ranking
  matches/
    ingest.py          match, segment and shared-match imports
    sides.py           maternal/paternal assignment
    cluster.py         Leeds-style clustering
    triangulate.py     shared-region grouping, chromosome painting
  tree/
    gedcom.py          GEDCOM 5.5.1 read/write, duplicate detection
    kinship.py         ancestors, descendants, relationship paths
  evidence.py          sources, citations, and the research planner
  budget.py            costed plans within a spending limit
  hypothesis.py        the investigation engine
  report.py            self-contained HTML reports
  demo.py              synthetic family generator with ground truth
  cli.py               command line interface
tests/                 100 tests, run with python3 -m unittest discover tests
docs/
  DATA_SOURCES.md      DNA testing sites: what to export, who can triangulate
  RECORD_SOURCES.md    record providers: what each costs and is good for
  AUTOMATION.md        what can be automated, and what a budget buys
```

Run the tests with:

```bash
python3 -m unittest discover -s tests -t .
```

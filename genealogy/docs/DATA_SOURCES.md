# Where the data comes from

You have two raw DNA files' worth of genotypes — yours and your mother's — and
that is all the raw DNA you will ever need. Everything else this toolkit
consumes is *derived* data that the testing sites compute for you and let you
export: who matches you, by how much, on which chromosomes, and which of your
matches also match each other.

This document is the practical answer to two questions: what do I download
from where, and which sites can actually triangulate.

---

## The short version

| What you need | Why | Where to get it |
|---|---|---|
| Your raw DNA | Everything | 23andMe, Ancestry (download) |
| **Your mother's raw DNA** | Splits every match maternal/paternal; lets you extract your father's DNA | Same |
| Match lists | The raw material | Every site |
| Segment data | Triangulation, chromosome mapping | 23andMe, MyHeritage, FTDNA, GEDmatch — **not Ancestry** |
| Shared-match lists | Clustering (works without segments) | Ancestry, 23andMe, MyHeritage, FTDNA |
| GEDCOM tree files | Names the ancestors the DNA points at | Ancestry, MyHeritage, FamilySearch, Geni, Geneanet |

The single highest-value action available to you: **upload both your kit and
your mother's kit to every site that accepts uploads.** It costs almost
nothing, and every site that holds both kits can then answer "does this person
match my mother too?" — which is the question that sorts your entire match
list into two piles.

---

## Site by site

### AncestryDNA — biggest database, no segment data

Roughly 25 million testers, far more than anywhere else, so this is where your
closest unknown relatives most likely are. The catch is that Ancestry has
never offered a chromosome browser: you get totals and shared-match lists, but
never "you share chr7:20–45 Mb". You therefore cannot triangulate on Ancestry
directly — you cluster instead, which is why `roots cluster` exists and why it
only needs in-common-with data.

- **Raw DNA download**: Settings → Download DNA Data. Ships a zip.
- **Match list / shared matches**: no official export. People use browser
  extensions to scrape it. Whatever tool you use, the importer here reads
  columns by meaning, so any sane CSV will load.
- **SideView**: Ancestry now labels many matches "Parent 1" / "Parent 2"
  automatically. It does not tell you *which* parent is which — but you know,
  because you have your mother's kit. If your export carries that column,
  `roots import-matches` picks it up as the `side` field.
- **ThruLines**: Ancestry's own guess at the common ancestor for a match. Very
  useful as a lead, entirely dependent on the accuracy of user-submitted
  trees. Treat as a hypothesis, not a finding.
- **Does not accept uploads** from other companies. Ancestry is a one-way door:
  data comes out, nothing goes in.

### 23andMe — segment data, and a real chromosome browser

- **Raw DNA download**: yes.
- **Relatives list export**: yes, a CSV of your DNA relatives.
- **Segment data**: the DNA Comparison tool gives you actual start/end
  positions. Export it and feed it to `roots import-segments`.
- **Shared matches**: visible per match ("Relatives in Common"), with the
  shared cM between *them*, which is genuinely triangulation-grade data.
- **Does not accept uploads.**
- Note 23andMe reports **percentages** rather than centimorgans in places. The
  percentage counts both chromosome copies, so a parent and a full sibling
  both show ~50% despite sharing very different amounts half-identically.
  Prefer a cM column; the importer warns when it has to convert.

### MyHeritage — free uploads, built-in triangulation

- **Accepts uploads** from 23andMe, Ancestry, FTDNA and others, free. Basic
  matching is free; the advanced DNA tools have historically required a
  one-off unlock fee or a subscription.
- **Chromosome browser**: yes, with a genuine triangulation view that shows
  where three or more people overlap on the same segment — not merely that
  they each overlap you.
- Upload both kits here. It is the cheapest route to real triangulation.

### FamilyTreeDNA — uploads, chromosome browser, and the Y

- **Accepts uploads** (autosomal) free, with a modest unlock fee for the full
  match list and chromosome browser.
- **Chromosome browser** and a **matrix** tool for pairwise comparison of your
  matches — that matrix is exactly the in-common-with data clustering wants.
- **Y-DNA testing.** This one deserves emphasis for your situation: if you are
  male, your Y chromosome came from your father, his father, and so on up an
  unbroken paternal line. A Y-37 or Big Y test can hand you a surname for a
  father you know nothing about, which no amount of autosomal analysis will do
  as directly. If you are female, a paternal-line male relative can test in
  your place — but identifying one is the very problem you are solving, so
  this is a later step.
- **mtDNA** traces the maternal line; less genealogically sharp because
  mitochondrial DNA mutates slowly, but it is definitive for confirming or
  excluding a maternal-line hypothesis.

### GEDmatch — the triangulation workhorse

Free tier plus a Tier 1 subscription (around $15/month, cancellable). Upload
raw data from any company. What it gives you that nothing else does:

- **One-to-many** and **one-to-one** comparison across every uploaded kit,
  regardless of which company tested them.
- **People who match both kits** — run it on you and your mother, and you have
  your maternal list; everyone else is paternal.
- **Triangulation** (Tier 1): finds groups that all match each other on the
  same segment.
- **Phasing utility**: give it a parent-child pair and it produces phased
  kits — precisely the maternal and paternal haplotype split that
  `roots phase` computes locally. Uploading a phased paternal kit and running
  one-to-many against it returns paternal matches only.
- **Lazarus** (Tier 1): reconstructs an untested person's genome from their
  relatives. Directly applicable to reconstructing your father.
- **Are your parents related?**: checks for runs of homozygosity, which matters
  because a positive result changes how you must interpret everything else.

Privacy note worth stating plainly: GEDmatch is an open-matching database and
has been used for law-enforcement identification. That is a real consideration
in deciding whether to upload, and it applies to your mother's kit as much as
your own — uploading her data is a decision about her, not only about you.

### DNA Painter — mapping and hypothesis testing

Free tier, paid subscription around $55/year. Not a matching database; a set
of tools.

- **Chromosome mapping**: assign segments to specific ancestors and watch your
  genome fill in. The same job `roots paint` does locally.
- **What Are The Odds (WATO)**: builds a hypothesis tree and scores where an
  unknown person fits given the cM they share with several known relatives.
  The closest thing on the market to what `roots investigate` does, and worth
  cross-checking your results against.
- **Shared cM Project tool**: the crowd-sourced relationship-probability table
  everyone cites. This toolkit computes its own distributions by simulation
  instead, and the two agree closely for close relationships — see the note in
  `roots/dna/simulate.py` about why they diverge at cousin level.

### Borland Genetics — reconstructing an untested parent

Free and paid tiers. Specialises in the exact manoeuvre your data supports:
tools that subtract a known parent from a child to produce the missing
parent's genome. `roots phase --export-paternal` produces the same artefact
locally, and you can upload the result there or to GEDmatch.

### Free trees, no DNA

**FamilySearch** (free, single shared world tree, excellent records),
**Geni**, **Geneanet**, **WikiTree** (free, and notable for taking DNA
evidence seriously in its structure). All export GEDCOM.

---

## What to actually download, in order

1. **Raw DNA** for you and your mother from 23andMe and Ancestry — all four
   files if you have them. Import each with `roots import-dna`, then
   `roots merge-kits` per person: the two vendors test overlapping but
   different SNP sets, so the union resolves segment boundaries better than
   either alone.
2. **Run `roots phase`** against your mother. This verifies she really is your
   mother from the data, then splits your genome into her half and your
   father's half, and exports the paternal half as an uploadable pseudo-kit.
3. **Match lists for both of you, from the same site.** This is the important
   constraint: absence from your mother's list only proves paternal descent if
   both lists come from the same database. `roots sides` refuses to infer
   across platforms for exactly this reason.
4. **Segment data** wherever it exists (23andMe, MyHeritage, FTDNA, GEDmatch)
   → `roots import-segments`.
5. **Shared-match lists** everywhere → `roots import-icw`. This is what
   clustering runs on, and it is your only analytical route on Ancestry.
6. **GEDCOM** from whichever site holds your tree → `roots import-tree`.
7. Link matches you have already identified to people in your tree with
   `roots link`, and record relationships you believe but have not proved with
   `roots claim` — the latter gets tested against the DNA rather than assumed.

---

## Recording relationships you already know

Two different situations, two different commands, and the distinction matters:

- `roots link --match 42 --xref I117` says "this match **is** this person in my
  tree". The relationship is then *derived* from the tree structure, and
  `roots investigate` checks whether the DNA agrees with it.
- `roots claim --match 42 --relationship 2c1r` says "I believe this match is my
  second cousin once removed" without placing them in the tree. The toolkit
  immediately tests that claim against the shared cM and tells you how well it
  fits.

Claims are treated as evidence to be tested, never as ground truth. That is
deliberate: an assumed relationship that is quietly wrong will corrupt every
inference built on top of it, and the whole point of having DNA is that you no
longer have to take a relationship on trust.

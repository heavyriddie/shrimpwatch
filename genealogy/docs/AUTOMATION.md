# What can be automated, and what a strict budget buys

The short answer: **the analysis is fully automatable and costs nothing. The
acquisition mostly is not, and that is where all the money goes.** The two
halves fail in opposite directions, so it is worth being precise about which
is which.

---

## Three tiers

### Tier 1 — already automated, £0, no network

Everything this toolkit does. Raw DNA import, cross-vendor merging,
segment detection, phasing against your mother, maternal/paternal
assignment, clustering, triangulation, relationship prediction, tree
reconciliation, the research plan, and the HTML report.

```bash
roots auto --budget 40 --report report.html
```

One command, every free and deterministic step, idempotent, safe to re-run.
It skips steps whose inputs are missing and says so rather than failing.
Put it on a cron job and it costs nothing forever, because it makes no
network requests at all — the entire system is one SQLite file on your disk.

The only meaningful cost is CPU: the relationship simulator takes a couple of
minutes on first run and is cached thereafter.

### Tier 2 — automatable if you build it, still £0

Two genealogy services have genuine public APIs:

| Service | API | Cost | What it could do for you |
|---|---|---|---|
| **FamilySearch** | Yes, OAuth, free developer key | Free | Search its record collections and read the shared tree programmatically. The largest free record set in existence. |
| **WikiTree** | Yes, public, no key for reads | Free | Check whether a match's surname and line already exist in a collaborative tree, with sources attached. |

These are the real automation opportunity, and both are free. A FamilySearch
integration could take the research plan this toolkit already produces and
attempt each lookup automatically — turning "find a marriage for Alfred
Whitlock around 1923" from a task into a result. That is a natural next
build, and it would not cost a penny to run.

It is not built yet. What exists is the plan it would consume.

### Tier 3 — not automatable at any price

| Service | Public API | Why not |
|---|---|---|
| Ancestry | No | Retired its public API years ago. Terms prohibit automated access, and it actively blocks it. |
| Findmypast | No | No public API. |
| 23andMe | No | Its genome API was withdrawn in 2018. |
| GEDmatch | No | None offered. |
| MyHeritage | Partner API only | DNA match data is not exposed through it. |

For all five, the workflow is: log in as a human, export, and hand the file
to `roots import-matches` / `import-segments` / `import-icw`. The importers
recognise columns by meaning precisely because the export formats are
whatever a browser extension happened to produce.

**This toolkit will not scrape these sites**, and that is a deliberate line
rather than a missing feature. Their terms forbid it, they detect and block
it, an account ban would cost you the match list you depend on — and the data
includes other living people who did not agree to be harvested. If someone
offers you a tool that automates Ancestry, understand what you are risking.

---

## What should stay manual even where it could be automated

Automation is not free of cost when it is wrong. Four decisions belong to a
human and the toolkit deliberately gates them:

**Buying a record.** Nothing spends money without you. `roots research`
proposes, you dispose.

**Deciding who a match is.** The toolkit produces probabilities and candidate
sets, not identifications. `roots claim` records your belief and immediately
tests it against the DNA rather than accepting it.

**Contacting a match.** These are living people. Automating outreach is how
you get blocked and how you upset relatives.

**Accepting a hypothesis into the tree.** `roots investigate` flags
contradictions; it does not resolve them. Merging duplicate individuals is
likewise proposed and never performed, because a bad automatic merge
propagates into every export you subsequently make.

---

## The budget model

Given a plan and a limit, `roots research --budget 40` chooses what to buy.
Three rules do the work:

**Free first, always.** Free tasks are never deferred and never counted
against the limit. The GRO birth index carries the mother's maiden surname
for nothing; the surviving Irish censuses are free; FamilySearch is free. Any
plan that spends before exhausting these is wrong at any budget.

**Subscription costs are shared, not per-record.** A census lookup and a
parish lookup and forty more cost *one month* between them. Charging each
separately makes subscription research look ruinous and pushes you toward
certificates you did not need. The planner costs two whole plans — one
taking a month's subscription, one not — and returns whichever yields more.
The practical corollary it prints: decide everything you need from a
subscription site, then take one month and do all of it.

**Then greedy by yield per pound.** A £3 GRO digital image that reveals a
mother's maiden name beats a £12.50 marriage certificate at four times the
price, even though the certificate is individually the more valuable
document. This is why the plans it produces tend to look cheap.

```
$ roots research --kit self-merged --budget 30

8 tasks selected (0 free, 8 paid), £21.99 of £30.00, 3 deferred

The plan includes one month of FindMyPast / Ancestry / FamilySearch at £9.99.
That fee covers every subscription task below, so do them all inside the same
month rather than spreading them out.
```

`roots budget --set 100` sets a cap, and `roots budget` reports actual spend
from the `cost` you record against each source — so the ledger reflects what
you really paid, not what was estimated.

Prices live in exactly one place, `roots/budget.py`, dated and sourced to
[RECORD_SOURCES.md](RECORD_SOURCES.md), overridable with `--cost-file`. If
that date is more than a year old, treat the plan as indicative.

---

## A realistic cadence

| When | What | Cost |
|---|---|---|
| Once | Import DNA, merge, phase, import tree | £0 |
| Weekly, cron | `roots auto` — re-clusters, re-sides, re-investigates as new matches arrive | £0 |
| Monthly | Re-export match lists from each site by hand, re-import | £0 |
| When a wall matters | `roots research --budget N`, buy the top items, record them with `roots source add --cost` | £3–£25 |
| Rarely | One month of a subscription, with a batched list of everything to pull | ~£10–£25 |

The honest bottom line on money: **the automation costs nothing and always
will.** A serious push on one English line back to 1837 runs about £105, of
which marriages are 60%; the same line is roughly £20 in Scotland and can be
free in Ireland. A £40 budget, spent in the order above, will typically buy
you two or three generations on the line you care about most.

The thing that no budget buys is the part DNA is uniquely good at: knowing
*which* line to spend on. That is what the free tier is for.

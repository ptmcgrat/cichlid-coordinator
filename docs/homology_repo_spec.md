# Homology decider repository: instruction sheet

*For the `cichlid-homology` repository (the name is proposed). October 2026. Written against cichlid-coordinator v0.11.0. Read first: `project_manifest.md` (why the resource exists), `species_final.md` (the species), `coordinator_spec.md` (layout, schemas, family contract, releases, the clade table), `trees_repo_spec.md` and `orthofinder_repo_spec.md` (what the two tree projects write), and the TreeFam project's `docs/OUTPUTS.md` (the per-gene tables it keeps).*

*The coordinator wrote this sheet from the project owner's goals and from what the other projects produce. Where it says "to settle in the pilot", the owner decides after the first results, and tells the coordinator if it changes anything here.*

## 1. Role

This repository makes **every orthology call** of the resource. It reads the families of the three methods (TreeFam, Compara, OrthoFinder) for a release and, for each **focal species** (the release's delivered species; release 1: *M. zebra* and *A. calliptera*), answers, gene by gene:

1. **Potential fusions:** which genes may be two genes merged.
2. **Very conserved proteins missing from the annotation:** a report of candidates only. Nothing is searched in the genome (the coordinator provides no genome sequence, and confirming them is not a goal).
3. **The closest human homolog, and the relationship:** one-to-one ortholog, one-to-two, and so on, or paralog.
4. **The closest zebrafish homolog, and the relationship:** the same.
5. **For genes with a one-to-one human ortholog:** the overall pattern in the teleosts. Are they mostly one-to-one?
6. **For genes with a one-to-two human ortholog:** do the copies look like **retained ohnologs** (from the teleost-specific genome duplication), or like **recent duplications** that belong at some level of the clade ladder?
7. **For paralogs:** at which level of the ladder did they arise?
8. **One decision per question across the three methods,** with every disagreement recorded.

It does **not** build trees or alignments (projects 1 and 2), download inputs, compute residues, variants or genotypes (project 4), or build the site (project 5).

**Rules (the same as for every repository):**

- Never build a data path by hand. Ask the coordinator (`get_layout()`, `cichlid-data where ...`).
- Never call rclone. Moving data is the coordinator's (`get_store()`, `cichlid-data pull | push | ensure`).
- Never download an input, and never copy the species list, the clade table or a species tree into this repository. Read them from the coordinator.
- Never write into an input folder.

Depend on the coordinator with a range pinned to one minor version, for example `cichlid-coordinator>=0.11,<0.12`, and use only its public API (`cichlid_data`) and its command line. If you need something it does not provide, write the coordinator a request that says what and why.

## 2. Repository layout (suggested)

```
cichlid-homology/
├── pyproject.toml           package cichlid_homology; depends on cichlid-coordinator>=0.11,<0.12
├── src/cichlid_homology/
│   ├── read.py              the families of each method, in one common in-memory form
│   ├── rooting.py           rooting a family tree (section 5)
│   ├── orthology.py         duplication and speciation nodes; closest homolog and relationship
│   ├── dating.py            the ladder level at which a duplication happened
│   ├── fusions.py, patterns.py, consensus.py, report.py, cli.py
├── tests/                   offline: small constructed trees
└── docs/
```

## 3. Inputs (all through the coordinator)

| Input | How to get it |
|---|---|
| Focal species, release, parts | `Catalog.load(get_layout()).releases["v1"]`: `deliver`, `parts`, `inputs` |
| Assignments, per method | `get_layout().results("v1", part).assign(prefix)`, for `part` in `TreeFam`, `Compara`, `OrthoFinder`. For TreeFam and Compara the `flags` column (`possible_fusion`, `low_family_coverage`, `close_runner_up`) and the fusion-partner columns are the fusion evidence; for OrthoFinder `best_family` is the orthogroup and the score columns are empty |
| Family trees and alignments | `results("v1", part).family(family).local`: `tree.nwk`, `combined.aln`, `done.json` (`tree_mode`: `place` or `denovo`; `rooted`) |
| Per-gene placement, TreeFam and Compara | the TreeFam project's `native/summary/gene_placement.tsv` (`gene_id, family, mode, status, lwr, n_placements`), `gene_status.tsv` and `skipped_families.tsv`. These are that project's own tables for now; **after your pilot, tell the coordinator which columns you read, and they become coordinator schemas** |
| The clade ladder | `Clades.load()`: `level(name)`, `common_level(names)`, `within(level)`, `lineage(name)`. Seven nested levels; `cyclostome` is **disputed** |
| The OrthoFinder species tree | `Catalog.load(layout).species_tree("v1", "OrthoFinder")` |
| A species tree for the other methods, and the master species tree | **To come from the coordinator** (section 5, rooting). Not yet available |
| Mappings of human and zebrafish identifiers across methods | **To come from the coordinator** (section 4). Not yet available |
| Reference gene names | the 2018 name tables, `Inputs/Annotations/Legacy2018` (a requirement of the Homology part) |
| Gene positions, for an optional synteny check | `layout.genes_table(prefix)` |

Run `cichlid-data check --release v1 --part Homology` before a run: it verifies the name tables and the clade table. The TreeFam and Compara results of release 1 are complete; **OrthoFinder's are not yet**, so start with the two that exist.

## 4. What each method can and cannot say

- **TreeFam and Compara, place mode.** Each focal gene is placed independently on a fixed reference tree. That gives its relationship to the reference genes (human, zebrafish and the rest) and the number of focal copies on a branch, but **it cannot order focal genes among themselves**: genes on one branch form a multifurcation. Only the de novo families (in release 1: 402 in TreeFam, 294 in Compara) are full trees of the focal genes.
- **OrthoFinder** builds a full tree of all 35 species, every cichlid included. It is the **authority for recent duplications among our own genes.**
- **Species coverage differs.** Teleost reference leaves: 8 in TreeFam, 37 in the pruned Compara trees (from the panel), 26 in OrthoFinder. The non-teleost bony fishes that anchor the teleost-specific duplication (spotted gar, reedfish, and bowfin if added): in Compara and OrthoFinder, **not in TreeFam**. Amphioxus: OrthoFinder only. Sea lamprey: TreeFam and Compara only.
- **Identifiers differ.** Human and zebrafish are Ensembl release-69-era proteins in TreeFam, current Ensembl proteins in Compara, and RefSeq genes in OrthoFinder. To compare methods, the closest human or zebrafish homolog must be mapped to common identifiers (symbols, Ensembl and NCBI gene IDs). The coordinator will provide the mappings as inputs (the TreeFam archive lists a file mapping its Ensembl IDs to RefSeq and UniProt; Compara needs an Ensembl-to-NCBI mapping for human and zebrafish). Until then, compare by gene symbol and say so.
- **Placement support varies.** In release 1, 4.7% of placed TreeFam genes and 6.4% of placed Compara genes have a best placement under 0.5 (the alternatives are in each family's `epa_result.jplace`). A gene in a de novo tree, or with no informative sites, is not placed.

## 5. Definitions (decided)

- **Relationship to a reference species R** (human or zebrafish): `1:1` (one R gene, one focal gene, no duplication since they separated), `1:many` (one R gene, several focal co-orthologs), `many:1`, `many:many`, **`paralog`** (the closest R gene is not an ortholog of the focal gene), or `none` (no R gene in the family).
- **Closest** means by position in the tree: the R genes in the smallest clade containing the focal gene and at least one R gene; branch length breaks ties.
- **Duplication or speciation:** by species overlap (a node is a duplication if its two sides share species; the threshold is to settle in the pilot). This needs no species tree.
- **Rooting:** **not midpoint rooting.** The TreeFam and Compara trees are unrooted on purpose. Root each tree by reconciling it with a species tree (the root that gives the fewest duplications), using the master species tree the coordinator will provide. OrthoFinder's trees may already be rooted (read `done.json`).
- **Dating a duplication:** the **common level** of the species below the duplication node, `Clades.common_level(species)`. Two cichlids meet at `cichlid`, a cichlid and a zebrafish at `teleost`, a cichlid and a human at `gnathostome`.
- **Retained ohnolog pair (teleost-specific):** both copies present across several teleost levels, with the non-teleost outgroup (gar, reedfish) holding one copy. **Recent duplication:** the second copy confined to one level or one species. Gene positions allow an optional synteny check for chromosome-level genomes.
- **Cyclostomes are disputed.** Do not date or conclude anything by the `cyclostome` level; flag it.
- **Fusions:** the union of the methods' evidence, kept as evidence and not collapsed: the TreeFam and Compara flags with the partner family and region, and, for OrthoFinder, a length outlier within the orthogroup. These are candidates, not proof (*Esur* shows how annotation artifacts inflate them).
- **Missing annotations:** families that are conserved across the other cichlids and teleosts but have no gene in the focal species. A report; no confirmation.
- **The teleost pattern of a one-to-one gene:** copies per teleost species and per level; the threshold for "mostly one-to-one" is to settle in the pilot. An absence may be an annotation gap, so use thresholds, not single species.

## 6. Combining the three methods

Disagreement is data, not noise.

1. **Same footing.** Each method fills the same per-gene call table (family, closest human and zebrafish gene with common identifiers, relationship, birth level, support), with a reason code for each call.
2. **Prefer the method built for the question.** Deep orthology to human and zebrafish: TreeFam and Compara, with OrthoFinder as a check. Recent duplications among our genes: OrthoFinder (or a de novo family). Teleost outgroup tests: Compara and OrthoFinder.
3. **Label the agreement** per question: `agree`, `majority`, `split`, or `single` (only one method could answer), and why: a different family, a different closest gene, a different relationship, or low support.
4. **Confidence** from the placement support, de novo status, family coverage and, where present, branch support.
5. **Independent tie-breakers:** alignment identity from `combined.aln`, and optionally conserved synteny.
6. **Keep every call.** The final call is stored with its status and the competing calls, and a **review list** collects the `split` cases.
7. **Calibrate** on a small benchmark of genes with known histories, to learn how often each method is right and to set the rules. The rules are code, with tests, and are versioned with the release.

## 7. Outputs

Under `get_layout().results("v1", "Homology")`:

| Path | Content |
|---|---|
| `run.json` | the run record: `release`, `method` (`Homology`), `method_version` (the version of the rules), `parameters` (every threshold, the rooting rule, the species-overlap setting), `started`, `finished`, `coordinator_version` (`cichlid_data.__version__`), `catalog_md5` (`Catalog.load(layout).md5()`), `proteomes_md5` (`Catalog.load(layout).proteomes_md5(species)`), `species`. Write it with `schemas.write_run` |
| the call tables | **their schema is designed with the coordinator after the pilot.** Until then write them under `native/calls/`. The existing `orthologs.tsv` schema has TreeFam-specific columns and will be replaced by one built for these eight questions |
| a review list | the `split` cases, under `native/review/` |

## 8. Pilot first

OrthoFinder is not ready, and the schema is not settled. Start with a **pilot** on TreeFam and Compara for one focal species (*M. zebra*), on a few hundred families: read the families, root them, classify the relationships, date the duplications, compare the two methods, and report what works and what does not. Settle the open thresholds (species overlap, "mostly one-to-one"), and propose the schema of the call tables. Only then run on all families and on OrthoFinder.

## 9. Validation and handing over

```bash
cichlid-data check --release v1 --part Homology             # the name tables and the clade table
R=$(cichlid-data where results v1 Homology)
cichlid-data validate "$R"/run.json
cichlid-data push 'results/v1/Homology/*' --dry-run         # then without --dry-run
```

The tests of this repository run offline on small constructed trees and must cover each case: a one-to-one ortholog, a one-to-two pair that is a retained teleost ohnolog, a recent duplication, a paralog, a fusion, a tie between methods, and an unrooted tree. There are no earlier results to compare with: the verification is the offline tests, the benchmark, `cichlid-data validate`, and a `run.json` that records the rules' version and every threshold.

## 10. What this repository must not do

- **No tree building, no alignment, no placement.** If a family needs rebuilding, that is project 1's or 2's.
- **No midpoint rooting** as the rule; no use of the disputed `cyclostome` level to date anything.
- **No genome searches,** no downloads, no hand-built paths, no copy of the species list, the clade table or a species tree.
- **No change to a family or an input file,** and nothing written into an input folder.
- **No change to a file format.** Ask the coordinator.

## 11. Acceptance criteria

1. The tests run offline and pass, and cover every case in section 9.
2. `cichlid-data check --release v1 --part Homology` passes.
3. `run.json` records the rules' version, every threshold and the rooting rule, and passes `cichlid-data validate`.
4. Every call has a reason code, and every disagreement between methods is recorded with the competing calls and listed for review.
5. The benchmark shows, for each question, how often each method agrees with the known answer.
6. The pilot's findings, and the schema of the call tables, are written down in `docs/` and agreed with the coordinator before the full run.
7. No other path, species list, clade table or cloud call appears anywhere in the code.

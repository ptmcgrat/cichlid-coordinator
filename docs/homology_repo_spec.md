# Homology decider repository: instruction sheet

*For the `cichlid-homology` repository (the name is proposed). October 2026. Written against cichlid-coordinator v0.14.1. Read first: `project_manifest.md` (why the resource exists), `species_final.md` (the species), `coordinator_spec.md` (layout, schemas, family contract, releases, the clade table), `trees_repo_spec.md` and `orthofinder_repo_spec.md` (what the two tree projects write), and the TreeFam project's `docs/OUTPUTS.md` (the per-gene tables it keeps).*

*The coordinator wrote this sheet from the project owner's goals and from what the other projects produce. Where it says "to settle in the pilot", the owner decides after the first results, and tells the coordinator if it changes anything here.*

## 1. Role

This repository makes **every orthology call** of the resource. It reads the families of the three methods (TreeFam, Compara, OrthoFinder) for a release and, for each **focal species** (the release's delivered species; release 1: *M. zebra* and *A. calliptera*), answers, gene by gene:

**Aim A, function inference.**
- **A1. The closest human gene, and the relationship:** one-to-one ortholog, one-to-two, and so on, or paralog.
- **A2. The closest zebrafish gene, and the relationship:** the same.

All three methods answer, **with equal standing**: none is preferred over another. A method whose family has no gene of the reference species **abstains**; that is not a disagreement. Every competing call is kept and listed for review.

**Aim B, the family in a teleost context.**
- **B1.** One copy kept (a duplicate lost early), both copies kept, or a later duplicate unrelated to the teleost genome duplication. Read it from the **copy number per teleost species** in the smallest clade holding the gene and human (or, when human is absent, another non-teleost anchor such as the gar), **not** from the position of single branches. Compara and OrthoFinder answer; TreeFam answers only as a fallback, is marked low resolution, and never overrules the other two.
- **B2.** For a later duplicate: how many cichlid species carry it, and where they meet in the master species tree. OrthoFinder is the authority; TreeFam and Compara give counts only, because place mode cannot order our own genes.

**Kept, at low priority:** potential fusions, and a report of very conserved proteins missing from the annotation (candidates only; nothing is searched in the genome, since the coordinator provides no genome sequence and confirming them is not a goal).

**Not done:** no comparison with the 2018 results (the *M. zebra* reference differs too much to link genes), and no dating of paralogs on the ladder as a goal of its own (the ladder is still used by B2). **Gene symbols** are not needed for the calls: they are for comparing methods and for checks.

**Mouse is not in scope for release 1:** human and zebrafish are the only reference species for the calls. Mouse stays a leaf in the trees of all three methods (as one of the vertebrate outgroups), but nothing is called against it and it needs no identifier mapping.

It does **not** build trees or alignments (projects 1 and 2), download inputs, compute residues, variants or genotypes (project 4), or build the site (project 5).

**Rules (the same as for every repository):**

- Never build a data path by hand. Ask the coordinator (`get_layout()`, `cichlid-data where ...`).
- Never call rclone. Moving data is the coordinator's (`get_store()`, `cichlid-data pull | push | ensure`).
- Never download an input, and never copy the species list, the clade table or a species tree into this repository. Read them from the coordinator.
- Never write into an input folder.

Depend on the coordinator with a range pinned to one minor version, for example `cichlid-coordinator>=0.14,<0.15`, and use only its public API (`cichlid_data`) and its command line. If you need something it does not provide, write the coordinator a request that says what and why.

## 2. Repository layout (suggested)

```
cichlid-homology/
├── pyproject.toml           package cichlid_homology; depends on cichlid-coordinator>=0.14,<0.15
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
| The species of every leaf | `LeafSpecies.load(get_layout(), part)` (`part` is `TreeFam`, `Compara` or `OrthoFinder`): `species(leaf)` gives a leaf's canonical species (the master tree's leaf name) or `None`, `of_leaves(leaves)` gives `{leaf: species}` and names any unknown leaf, `is_project_gene(leaf)`, `shadowed()`. Reference leaves come from the part's own mapping file (TreeFam's `protein_to_species.txt`, Compara's `prepared/protein_to_species.txt`), our genes by catalog prefix (`<prefix>_<id>`). The reference mapping is looked up first. `LeafSpecies.load(..., only=set_of_leaves)` keeps only those leaves, to bound memory |
| Per-gene placement, TreeFam and Compara | the TreeFam project's `native/summary/gene_placement.tsv` (`gene_id, family, mode, status, lwr, n_placements`), `gene_status.tsv` and `skipped_families.tsv`. These are that project's own tables for now; **after your pilot, tell the coordinator which columns you read, and they become coordinator schemas** |
| The clade ladder | `Clades.load()`: `level(name)`, `common_level(names)`, `within(level)`, `lineage(name)`. Seven nested levels; `cyclostome` is **disputed** |
| The OrthoFinder species tree | `Catalog.load(layout).species_tree("v1", "OrthoFinder")` |
| The master species tree, for rooting | `cichlid_data.master_tree.load_master_tree()`: one rooted tree over the 154 species that are leaves in any method's trees (the agreed OrthoFinder tree is its backbone). `master_tree.canonical(name)` gives a species' leaf name (TreeFam's mapping file and its tree spell a few species differently) |
| The names and official IDs of the reference proteins (human, zebrafish) | `ReferenceNames.load(get_layout(), part)` for `TreeFam` and `Compara`: `get(protein)` returns the gene, symbol, **official ID (HGNC or ZFIN)**, `symbol_source` and a `status` (`current`, `grch37`, `legacy_id`, `replaced`, `legacy_symbol`, `unmapped`); `reliable` is true for the first three. A frozen snapshot of the earlier pipeline's stage 6a, so the names are as of the day it was made. For OrthoFinder there is no table: its human and zebrafish genes are our own proteomes, with the symbol and NCBI gene ID in the genes tables; their HGNC and ZFIN IDs are a later addition |
| The 2018 name tables | `Inputs/Annotations/Legacy2018` is still a requirement of the Homology part (the reference names were built partly from it), but this repository does not read it |
| Gene positions, for an optional synteny check | `layout.genes_table(prefix)` |

Run `cichlid-data check --release v1 --part Homology` before a run: it verifies the name tables and the clade table. The TreeFam and Compara results of release 1 are complete; **OrthoFinder's are not yet**, so start with the two that exist.

## 4. What each method can and cannot say

- **TreeFam and Compara, place mode.** Each focal gene is placed independently on a fixed reference tree. That gives its relationship to the reference genes (human, zebrafish and the rest) and the number of focal copies on a branch, but **it cannot order focal genes among themselves**: genes on one branch form a multifurcation. Only the de novo families (in release 1: 402 in TreeFam, 294 in Compara) are full trees of the focal genes.
- **OrthoFinder** builds a full tree of all 35 species, every cichlid included. It is the **authority for recent duplications among our own genes.** Its family trees are expected to be unrooted IQ-TREE trees with branch lengths and **no support values**, one per orthogroup (the choice of tree is decided with the owner after its pilot).
- **Species coverage differs.** Teleost reference leaves: 8 in TreeFam, 37 in the pruned Compara trees (from the panel), 26 in OrthoFinder. The non-teleost bony fishes that anchor the teleost-specific duplication (spotted gar, reedfish, and bowfin if added): in Compara and OrthoFinder, **not in TreeFam**. Amphioxus: OrthoFinder only. Sea lamprey: TreeFam and Compara only.
- **Identifiers differ across methods:** Ensembl release-69-era proteins in TreeFam (many retired), current Ensembl proteins in Compara, RefSeq genes in OrthoFinder. The common key is the **official ID** (HGNC for human, ZFIN for zebrafish), from `ReferenceNames` for TreeFam and Compara. Compara's proteins are current, so nearly all resolve at once; TreeFam's older ones are resolved by several routes (in release 1: for zebrafish, 26,558 current, 11,124 by legacy ID and 4,417 unmapped of 42,146; for human, 79,001 current, 10,555 by GRCh37, 2,152 by legacy ID and 1,391 unmapped of 93,165), so **a call that rests on a `replaced`, `legacy_symbol` or `unmapped` protein is lower confidence, or the method abstains**. For OrthoFinder, compare by symbol and say so until its IDs are added.
- **Some reference species are also project species** (release 1):

| Method | Reference species that are also project species | Treatment in the trees |
|---|---|---|
| TreeFam | `gasterosteus_aculeatus` (Gacu), `oreochromis_niloticus` (Onil) | Both stay as Ensembl-annotated reference leaves **and** are placed as our RefSeq genes: shadowed |
| Compara | `betta_splendens` (Bspl), `gasterosteus_aculeatus` (Gacu) | Same: reference leaves (role `panel`) and placed genes: shadowed. The other five (*O. niloticus*, *N. brichardi*, *H. burtoni*, *M. zebra*, *A. calliptera*) have role `validation`: their Ensembl copies were removed from the trees when they were prepared, so they are not in the trees at all |
| OrthoFinder | none | Every species appears once, as our proteome |

  So in TreeFam and Compara a tree can hold Ensembl's stickleback (and tilapia, or betta) next to our placed genes of the same species.
- **De novo families differ by method** (your own survey of 100 per method that hold an *M. zebra* gene): Compara's de novo trees are mostly fish-only (median 10 species; human is present in 1% of them, zebrafish in 25%), TreeFam's are broader (median 28 species; human in 54%, zebrafish in 12%). So a gene in a Compara de novo family usually cannot get its human relationship from that method, and one in a TreeFam de novo family often cannot get its zebrafish one.
- **In a multifurcation** (genes placed on one branch) the order of the genes in the Newick text means nothing, and the branch lengths inside the group are not to be relied on. Use `placements.tsv`, `epa_result.jplace` and `combined.aln` for tie-breaks.
- **Placement support varies.** In release 1, 4.7% of placed TreeFam genes and 6.4% of placed Compara genes have a best placement under 0.5 (the alternatives are in each family's `epa_result.jplace`). A gene in a de novo tree is a leaf of that tree, not a placed gene. **A gene with no informative sites** has `mode` `place`, blank `lwr` and `n_placements`, **no leaf in `tree.nwk`** and no row in `placements.tsv`, but an all-gaps row in `combined.aln` (and an all-`-` row in `residue_map.tsv`): it has no position in the tree, so it gets no relationship from that method. **The rate varies by species:** in TreeFam 3.0% to 4.6% for ten species but *Esur* 9.9% and betta 12.9%; in Compara 3.8% to 6.7% for eleven species but *Esur* 17.3%; *M. zebra* is ordinary (3.4% and 5.4%). The shadowed species have the lowest rates where they are shadowed (stickleback in both, betta in Compara), which may be support inflated by placing a gene beside its own Ensembl copy (a hypothesis). So judge support against the species' and the method's own distribution, and take care with *Esur*, betta and stickleback.

## 5. Definitions (decided)

- **Relationship to a reference species R** (human or zebrafish): `1:1` (one R gene, one focal gene, no duplication since they separated), `1:many` (one R gene, several focal co-orthologs), `many:1`, `many:many`, **`paralog`** (the closest R gene is not an ortholog of the focal gene), or `none` (no R gene in the family). **`none` means no R gene in this method's family, not that no homolog exists:** in a fish-only de novo family it says nothing about human. Compare the methods before concluding that a gene has no human or zebrafish homolog.
- **Closest** means by position in the tree: the R genes in the smallest clade containing the focal gene and at least one R gene; branch length breaks ties.
- **Duplication or speciation:** by species overlap (a node is a duplication if its two sides share species; the threshold is to settle in the pilot). This needs no species tree.
- **Shadowed species are ignored for duplication calls.** A reference leaf whose species is also a project species (the table in section 4; `LeafSpecies.shadowed()`) is dropped before species overlap is computed, and our placed genes stand for the species. Ensembl's copy and our RefSeq gene are two annotations of the same genome, so counting them as shared species would call a duplication at every such pair. This is the rule the coordinator already applied to the five cichlids (role `validation`). Record it in `run.json`; the owner can overrule it.
- **Rooting:** **not midpoint rooting.** The TreeFam and Compara trees are unrooted on purpose. Root each tree by reconciling it with a species tree (the root that gives the fewest duplications), using the master species tree (`master_tree.load_master_tree()`). OrthoFinder's trees may already be rooted (read `done.json`).
- **Dating a duplication:** the **common level** of the species below the duplication node, `Clades.common_level(species)`. Two cichlids meet at `cichlid`, a cichlid and a zebrafish at `teleost`, a cichlid and a human at `gnathostome`.
- **Retained ohnolog pair (teleost-specific; Aim B1):** both copies present across several teleost levels, with the non-teleost outgroup (gar, reedfish) holding one copy. **Recent duplication:** the second copy confined to one level or one species. Gene positions allow an optional synteny check for chromosome-level genomes.
- **B2, where cichlid carriers meet:** the last common ancestor of the carrying species in the master species tree. The coordinator's `master_tree` has no public function for it: compute it from the tree (`newick.Node`: `iter`, `leaves`, `parent`). `Clades.common_level` gives only the ladder level, which inside the cichlids is always `cichlid`. *Pundamilia* sits in the haplochromine node with *Haplochromis*, *M. zebra* and *A. calliptera*, unresolved among them on purpose (published trees disagree), and *Amphilophus* is beside *Archocentrus* in the Neotropical pair. A duplicate carried by *M. zebra*, *A. calliptera* and *Pundamilia* therefore meets at the four-species haplochromine node, not at the whole cichlid clade. Several nodes are deliberately left as polytomies where published work disagrees (see the list in the coordinator's spec): a speciation that falls on one of them can read as a duplication, so treat a call that rests on such a node as lower confidence. Each placement has its source in `config/species_tree_placements.tsv`.
- **Cyclostomes are disputed.** Do not date or conclude anything by the `cyclostome` level; flag it.
- **Fusions:** the union of the methods' evidence, kept as evidence and not collapsed: the TreeFam and Compara flags with the partner family and region, and, for OrthoFinder, a length outlier within the orthogroup. These are candidates, not proof (*Esur* shows how annotation artifacts inflate them). `possible_fusion` is a lower-ranked family with coverage of at least 0.5 whose region on the gene overlaps the best family's by at most 20% of the shorter; the partner is the highest-scoring such family. The four `fusion_partner*` columns are filled if and only if `possible_fusion` is set (guaranteed by project 1's code and a test). `low_family_coverage` (best-family coverage below 0.5) means a different coverage in each method (TreeFam: of the HMM; Compara: of the best reference protein), so compare it only within a method.
- **Missing annotations:** families that are conserved across the other cichlids and teleosts but have no gene in the focal species. A report; no confirmation.
- **B1, the copy number:** copies per teleost species in the smallest clade holding the gene and human (or the gar, when human is absent), read as a pattern across species and levels, not from single branches. The thresholds are to settle in the pilot. An absence may be an annotation gap, so use thresholds, not single species.

## 6. Combining the three methods

Disagreement is data, not noise.

1. **Same footing.** Each method fills the same per-gene call table (family, closest human and zebrafish gene with common identifiers, relationship, birth level, support), with a reason code for each call.
2. **Equal standing for A1 and A2.** No method is preferred over another; a method whose family has no gene of the reference species abstains, which is not a disagreement. For B1, Compara and OrthoFinder answer and TreeFam only as a low-resolution fallback that never overrules them; for B2, OrthoFinder is the authority and the other two give counts only.
3. **Label the agreement** per question: `agree`, `majority`, `split`, or `single` (only one method could answer), and why: a different family, a different closest gene, a different relationship, or low support.
4. **Confidence** from the placement support, de novo status, family coverage and, where present, branch support.
5. **Independent tie-breakers:** alignment identity from `combined.aln`, and optionally conserved synteny.
6. **Keep every call.** The final call is stored with its status and the competing calls, and a **review list** collects the `split` cases.
7. **Calibrate later, if a benchmark can be identified.** A benchmark of genes with known histories is hard to generate; it is left to this project and implemented only if such genes can be found, with a source and citation for each. The rules are code, with tests, and are versioned with the release.

## 7. Outputs

Under `get_layout().results("v1", "Homology")`:

| Path | Content |
|---|---|
| `run.json` | the run record: `release`, `method` (`Homology`), `method_version` (the version of the rules), `parameters` (every threshold, the rooting rule, the species-overlap setting), `started`, `finished`, `coordinator_version` (`cichlid_data.__version__`), `catalog_md5` (`Catalog.load(layout).md5()`), `proteomes_md5` (`Catalog.load(layout).proteomes_md5(species)`), `species`. Write it with `schemas.write_run` |
| the call tables | **their schema is designed with the coordinator after the pilot.** Until then write them under `native/calls/`. The existing `orthologs.tsv` schema has TreeFam-specific columns and will be replaced by one built for these questions |
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

The tests of this repository run offline on small constructed trees and must cover each case: a one-to-one ortholog, a one-to-two pair that is a retained teleost ohnolog, a recent duplication, a paralog, a fusion, a tie between methods, and an unrooted tree. There are no earlier results to compare with: the verification is the offline tests, a benchmark if one can be identified, `cichlid-data validate`, and a `run.json` that records the rules' version and every threshold.

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
5. If a benchmark of genes with known histories could be identified, it shows, for each question, how often each method agrees with the known answer.
6. The pilot's findings, and the schema of the call tables, are written down in `docs/` and agreed with the coordinator before the full run.
7. No other path, species list, clade table or cloud call appears anywhere in the code.

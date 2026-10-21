# TreeFam and Compara repository: instruction sheet

*For the `cichlid-trees` repository. October 2026. Written against cichlid-coordinator v0.9.0 (releases: results are `Results/<release>/<part>/`). Read first: `project_manifest.md` (why the project exists), `species_final.md` (the species), and `coordinator_spec.md` (the layout, the schemas and the family contract, which this sheet relies on).*

*What the coordinator knows about this repository comes from its own documents and from the real results folders; it has not seen the repository's code. Where this sheet says "to confirm", the owner of this repository decides, and tells the coordinator if it changes anything here.*

## 1. Role

This repository builds the TreeFam 9 and Ensembl Compara 116 results. For each **placed** species (in release 1, the 12 project species) it assigns every protein to a family (TreeFam by HMM search, Compara by sequence search), aligns the proteins with the family's reference alignment, places them on the family's fixed reference tree, and writes what it finds in the project's standard formats. It owns:

1. **Assignment:** which family each protein belongs to, as `assign/<prefix>.assign.tsv`.
2. **Alignment and tree construction:** for each family, `combined.aln` and the final `tree.nwk` (stage 5), following the family contract.
3. **Its own working files** beside those (`queries.fa`, `ref.tree`, `epa_result.jplace`, `placements.tsv`, and so on), which are not part of the contract.

It does **not**:

- download or prepare inputs, build proteomes, or move data to and from the cloud (the coordinator);
- make orthology calls, combine methods, or name genes (the homology decider, project 3);
- compute residue maps, amino-acid changes or genotypes (the analysis project, project 4);
- build the site (project 5).

Only project 3 makes orthology calls, for all three methods with one rule. Calling and combining code written for earlier versions of this repository (the stages that call orthologs, reference gene names, and combined calls) moves to project 3.

**Rules (the same as for every repository):**

- Never build a data path by hand. Ask the coordinator (`get_layout()`, `cichlid-data where ...`).
- Never call rclone. Moving data to and from the cloud is the coordinator's (`get_store()`, `cichlid-data pull | push | ensure`).
- Never download an input yourself.
- Never copy the species list into this repository. Read it from the catalog.

Depend on the coordinator with a range pinned to one minor version, for example `cichlid-coordinator>=0.12,<0.13` (the minor you have tested), and use only its public API (`cichlid_data`) and its command line. Below 1.0.0 a minor version may change the layout API (v0.9.0 did: results are now `Results/<release>/<part>/`). When the coordinator moves to a new minor version, run your tests and `cichlid-data validate` against it, then move your pin to it; `cichlid_data.__version__` is the installed version.

## 2. What leaves this repository

Moving code out follows the coordinator's rule: it moves verbatim with its tests, changing only paths, transfers and formats; the outputs on real data are checked to be identical before it is accepted, and any difference is explained (section 6).

| Today in this repository | Goes to | Replaced by |
|---|---|---|
| Stage 1, genome and annotation download | the coordinator (done) | `cichlid-data download genomes` |
| Stage 2, proteomes | the coordinator (done) | `cichlid-data proteomes`; `layout.proteome(prefix)` |
| `newick.py` | the coordinator (done) | `cichlid_data.newick` |
| `compara_prepare.py` | the coordinator (done) | `cichlid-data prepare compara` |
| The download and summary scripts (annotations, expression, populations, h5ad and RNA-seq summaries) | the coordinator (done) | `cichlid-data download ...` and `summarize ...` |
| `export_seurat.R` | the coordinator (v0.9.1): it is the export step after `download expression --source daniocell`, and `summarize h5ad` reads its output | the packaged script `cichlid_data/downloaders/export_seurat.R` |
| Its own data-sync code | the coordinator | the store (`get_store()`) and `cichlid-data pull | push | ensure` |
| Ortholog calling, reference gene names (stage 6a) and combined calls (stage 7) | project 3 | per `cichlid-homology` |
| Residue tables (`residues.tsv.gz`, `residues.done.json`; stage 8) | project 4 | per `cichlid-analysis`; `residue_map.tsv` stays here (section 4) |
| The browser site (`web/`) | project 5 | per `cichlid-server` |
| Combined calls (stage 7) | project 3 | as above |
| `analyses` (variants; tissue, enrichment, gene sets), stage 9 | set aside, untouched | variant analyses are project 4's; the owner of the others is not decided; stage 9 is to be described by this project |

Stages are named as the coordinator's documents name them. **Stage 5 stays here:** it makes `tree.nwk`. In *place* mode it prunes the TreeFam reference tree to the alignment's rows, unroots it and resolves polytomies, evaluates branch lengths (raxml-ng), places the queries (EPA-ng) and grafts them (gappa). In *de novo* mode, for a family with no usable reference tree (286 of 8,930 TreeFam families in the previous run; in release 1, 402 of 9,090 TreeFam and 294 of 15,045 Compara families), it builds the tree from `combined.aln` (raxml-ng). *(As reported by this project.)* TreeFam's `family_audit.tsv` (stage 3a) also stays here: its statuses depend only on TreeFam 9's files. Compara's is made by the coordinator.

## 3. Inputs (all through the coordinator)

| Input | How to get it |
|---|---|
| The species | `Catalog.load(get_layout()).species(release="v1", part="TreeFam")` (and `part="Compara"`): the placed species, in the catalog's order |
| The proteomes | `get_layout().proteome(prefix)`: one protein per gene, headers `>gene_id`, already valid |
| TreeFam 9 | `get_layout().family_resource("TreeFam9")`: `hmm/`, `trees/<family>.aln` and `.nhx`, `protein_to_species.txt`, `source.json` |
| Compara 116 | `get_layout().family_resource("Compara116")`: `prepared/` (`trees/`, `search_db.fa`, `members.tsv`, `protein_to_species.txt`, `family_audit.tsv`, `prepare.json`); the Compara reference panel is `config/compara_panel.tsv`, owned by the coordinator |

Run `cichlid-data check --release v1 --part TreeFam --fetch` (and `--part Compara`) before every run; stop if it fails. The 2018 name tables (`Inputs/Annotations/Legacy2018`) are no longer part of this check: they belong to the Homology part.

**A fact about TreeFam 9's own files** (checked on family TF101001, and every one of the 8,836 trees parses, and the 8,836 alignments paired with them read, without error): a tree's leaf names are bare protein IDs, while the alignment's row names carry a `/start-end` suffix (`FBpp0177754/1-356`). Stripping `/\d+-\d+$` makes the two sets equal. Match them that way, and write the rows of your own `combined.aln` under exactly the tree's leaf names, as the existing results do.

**Per-gene fields for project 3.** Your `docs/OUTPUTS.md` documents `native/summary/gene_placement.tsv`, `gene_status.tsv` and `skipped_families.tsv`. They stay your own files for now. When project 3's sheet says which of them it reads, those become coordinator schemas, and you will be told.

**Inputs are read-only, and the DIAMOND database is yours.** `prepared/` holds `search_db.fa`, the unaligned sequences, which does not depend on any tool. A DIAMOND database depends on the DIAMOND version, so this repository builds its own, in a local cache (keyed by the DIAMOND version and the checksum of `search_db.fa`), and never writes into an input folder. The `search_db.dmnd` that sits in `prepared/` was written by `ensure_db` into the input folder and swept into the unit; ignore it. The coordinator's `prepare compara` does not make or recreate it.

**Compara's `family_audit.tsv` statuses** are `ok`, `tiny_tree`, `cichlid_only`, `no_panel_members`, `mismatch` and `tree_parse_error`, as the preparation script documents them; they were produced by the script when it lived in this repository. The Compara path here handles all six.

**Placed species and delivered species.** Release 1 places the 12 project species. Assignment tables are written for every placed species. Release 1 delivers calls for *M. zebra* and *A. calliptera* only; that is project 3's step, and it is settled in the release (`deliver`), not here.

## 4. Outputs

Everything is under `get_layout().results("v1", "TreeFam")` and `results("v1", "Compara")` (a `ResultsLayout`; `v1` is the release, which is the coordinator's to name, and the part is `TreeFam` or `Compara`). Use the API, not the names.

| Path | Content |
|---|---|
| `assign/<prefix>.assign.tsv` | one per placed species, the existing assignment schema (`best_family` is the TreeFam or Compara family; the score, coverage and fusion columns are yours) |
| `families/<family>/` | one folder per family, the family contract below |
| `native/<subfolder>/` | this part's own outputs: stage 3's hits, parameters and summaries, TreeFam's `family_audit.tsv` (stage 3a, `native/audit/`) and tables such as `ensembl_xref.tsv` (`native/xref/`), each kind in its own subfolder (`native/hits/`, `native/params/`, `native/summary/`). Each top-level folder is stored as its own tar; loose files go together into one extra tar that is repacked whenever any of them changes, so use subfolders |
| `run.json` | the run record: `release` (`v1`; it must match the folder it is in), `method` (`TreeFam` or `Compara`), `method_version`, `parameters` (search thresholds, aligner, placement program and their versions), `started`, `finished`, `coordinator_version` (`cichlid_data.__version__`), `catalog_md5` (`Catalog.load(layout).md5()`), `proteomes_md5` (`Catalog.load(layout).proteomes_md5(species)`, from coordinator v0.10: it ties the run to exactly the proteomes it used), `species` (the placed prefixes). Write it with `schemas.write_run`. |

### The family contract (summary; the full text is in `coordinator_spec.md`)

For every family, a folder named by the family ID (`TF101001`, `C116_012345`) holding:

- **`tree.nwk`:** the final tree. For a family with a usable reference tree, it is the reference tree, pruned to the alignment's rows and unrooted, with your genes placed on it; for a family without one, it is a tree built from `combined.aln`. Leaf names are unique: our gene IDs (`<prefix>_<id>`) for the placed proteins, and the reference proteins' own IDs for the rest.
- **`combined.aln`:** the aligned FASTA, with a row for every leaf of the tree, rows named like the leaves.
- **`done.json`:** a small JSON object, written **last**, at the end of stage 5: a folder without it is unfinished. It holds `"release": "v1"`, `"rooted": false` (stage 5 unroots the reference tree on purpose, and de novo trees are unrooted, so project 3 roots them by reconciliation with the master species tree; Compara's species-tree-aware root is lost, and keeping it would be a separate decision for project 3), `"tree_mode"` (`place` or `denovo`, as you reported; the coordinator will validate them), and the inputs hash that lets a rerun skip finished work. Name the stage records `stage4.json` and `stage5.json` (they are yours); the per-gene status list stays in `stage5.json` until project 3's sheet defines what it reads from a family.
- **`residue_map.tsv`:** (its format is now a coordinator schema, checked by `validate`; see `coordinator_spec.md`) keep writing it, from MAFFT `--mapout`, for the **delivered species of the release** (read them from the catalog: `Catalog.load(layout).releases["v1"]["deliver"]`), not from a hard-coded list. `combined.aln` alone cannot give an exact map here, because `--keeplength` drops query-only insertion columns. Send the coordinator the file's header and what each column means: it will become a schema, and project 4 will read it.
- **Skipped families have no folder.** A family with fewer than 4 sequences is skipped (in release 1: 10 in TreeFam and 116 in Compara); its genes keep their `best_family` in `assign/`. Record the skipped families, with their sequence counts, in `native/summary/` (say so in your documentation if you already do), so that project 3 can tell a skipped family from a missing one.
- **`stage5/` is scratch.** Keep each family's `stage5/` working folder locally as you wish: the `families` unit does not push anything inside it (the `stage5.json` record beside it is pushed).
- The other files you keep beside these (`queries.fa`, `ref.tree`, `epa_result.jplace`, `placements.tsv`, ...) are yours and are not checked. `residues.tsv.gz` and `residues.done.json` (stage 8) belong to project 4: leave them as they are and do not extend them.

## 5. Validation and handing over

```bash
cichlid-data check --release v1 --part TreeFam --fetch                 # before a run (and again with --part Compara)
R=$(cichlid-data where results v1 TreeFam)                             # never type the data root
cichlid-data validate "$R"/assign/*.assign.tsv                         # one per placed species
cichlid-data validate "$R"/families                                    # every family folder, leaf names checked against the catalog
cichlid-data validate "$R"/run.json
cichlid-data push 'results/v1/TreeFam/*' --dry-run                     # then without --dry-run
```

Hand over to project 3 only when all of these pass, for both methods.

## 6. Verification

There is **no comparison with the old results.** The earlier plan to reproduce `Rebuild2026` and `Compara2026` for the nine species was dropped, by decision of the project owner: the effort is not worth it, and a 12-species run aligns all of a family's queries together, so it could not match a 9-species run file by file anyway. In its place: this repository's own offline tests, `cichlid-data validate` on every output (section 5), and a `run.json` that records the version of every program used (MAFFT, DIAMOND, pyhmmer, raxml-ng, EPA-ng, gappa), so the results can be reproduced from now on. The old results folders may be kept or deleted at the owner's discretion; nothing here depends on them.

## 7. What this repository must not do

- **No orthology calls** and no combining of methods; those are project 3's.
- **No variant or genotype work, and no residue tables;** those are project 4's (you keep `residue_map.tsv`, section 4).
- **No hand-built paths,** no copy of the species list, no rclone, no downloads, no editing of the proteomes, the catalog or the TreeFam and Compara inputs.
- **No change to the TreeFam 9 or Compara reference sets** to make a run work. If they need changing, that is a decision recorded in the coordinator (`config/compara_panel.tsv`, and the input versions of a release).
- **No change to a file format.** If the contract is not enough, ask the coordinator; a schema change is a version change.

## 8. Acceptance criteria

1. The tests run offline and pass.
2. `cichlid-data check` passes for both methods, and every run is validated as in section 5.
3. `run.json` records the version and parameters of every program used, and the coordinator version and catalog checksum, and passes `cichlid-data validate`.
4. Every assignment file and every family folder passes `cichlid-data validate`, for the 12 placed species; every `done.json` is written last and says `release`, `rooted` and `tree_mode`.
5. The results are pushed, and `cichlid-data status 'results/v1/TreeFam/*'` (and `v1/Compara`) shows everything in sync.
6. No other path, species list, or cloud call, no calling code, and no residue code other than `residue_map.tsv`, appears anywhere in the repository; nothing is written into an input folder.

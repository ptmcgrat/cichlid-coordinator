# TreeFam and Compara repository: instruction sheet

*For the `cichlid-trees` repository. October 2026. Written against cichlid-coordinator v0.8.0; the release layout planned for v0.9.0 is noted where it changes a path. Read first: `project_manifest.md` (why the project exists), `species_final.md` (the species), and `coordinator_spec.md` (the layout, the schemas and the family contract, which this sheet relies on).*

*What the coordinator knows about this repository comes from its own documents and from the real results folders; it has not seen the repository's code. Where this sheet says "to confirm", the owner of this repository decides, and tells the coordinator if it changes anything here.*

## 1. Role

This repository builds the TreeFam 9 and Ensembl Compara 116 results. For each **placed** species (in release 1, the 12 project species) it assigns every protein to a family (TreeFam by HMM search, Compara by sequence search), aligns the proteins with the family's reference alignment, places them on the family's fixed reference tree, and writes what it finds in the project's standard formats. It owns:

1. **Assignment:** which family each protein belongs to, as `assign/<prefix>.assign.tsv`.
2. **Alignment and placement:** for each family, `combined.aln` and the final `tree.nwk`, following the family contract.
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

Depend on the coordinator with a range pinned to one minor version, for example `cichlid-coordinator>=0.8,<0.9`, and use only its public API (`cichlid_data`) and its command line. Below 1.0.0 a minor version may change the layout API; v0.9.0 is planned to (section 4).

## 2. What leaves this repository

Moving code out follows the coordinator's rule: it moves verbatim with its tests, changing only paths, transfers and formats; the outputs on real data are checked to be identical before it is accepted, and any difference is explained (section 6).

| Today in this repository | Goes to | Replaced by |
|---|---|---|
| Stage 1, genome and annotation download | the coordinator (done) | `cichlid-data download genomes` |
| Stage 2, proteomes | the coordinator (done) | `cichlid-data proteomes`; `layout.proteome(prefix)` |
| `newick.py` | the coordinator (done) | `cichlid_data.newick` |
| `compara_prepare.py` | the coordinator (done) | `cichlid-data prepare compara` |
| The download and summary scripts (annotations, expression, populations, h5ad and RNA-seq summaries) | the coordinator (done) | `cichlid-data download ...` and `summarize ...` |
| Its own data-sync code | the coordinator | the store (`get_store()`) and `cichlid-data pull | push | ensure` |
| Ortholog calling, reference gene names (stage 6a) and combined calls (stage 7) | project 3 | per `cichlid-homology` |
| Residue maps and residue tables | project 4 | per `cichlid-analysis` |
| The browser site | project 5 | per `cichlid-server` |

Stages are named as the coordinator's documents name them. **To confirm:** what stage 5 does (it leaves a `stage5/` folder and `stage5.json` in each family) and whether it stays here or moves.

## 3. Inputs (all through the coordinator)

| Input | How to get it |
|---|---|
| The species | `Catalog.load(get_layout()).species(method="TreeFam2026")` (and `"Compara2026"`): the placed species, in the catalog's order |
| The proteomes | `get_layout().proteome(prefix)`: one protein per gene, headers `>gene_id`, already valid |
| TreeFam 9 | `get_layout().family_resource("TreeFam9")`: `hmm/`, `trees/<family>.aln` and `.nhx`, `protein_to_species.txt`, `source.json` |
| Compara 116 | `get_layout().family_resource("Compara116")`: `prepared/` (`trees/`, `search_db.fa`, `members.tsv`, `protein_to_species.txt`, `family_audit.tsv`); the Compara reference panel is `config/compara_panel.tsv`, owned by the coordinator |

Run `cichlid-data check --method TreeFam2026 --fetch` (and `--method Compara2026`) before every run; stop if it fails. For now it also checks the 2018 name tables (`Inputs/Annotations/Legacy2018`), which belong to project 3, so that requirement moves there.

**A fact about TreeFam 9's own files** (checked on family TF101001, and every one of the 8,836 trees parses, and the 8,836 alignments paired with them read, without error): a tree's leaf names are bare protein IDs, while the alignment's row names carry a `/start-end` suffix (`FBpp0177754/1-356`). Stripping `/\d+-\d+$` makes the two sets equal. Match them that way, and write the rows of your own `combined.aln` under exactly the tree's leaf names, as the existing results do.

**Placed species and delivered species.** Release 1 places the 12 project species. Assignment tables are written for every placed species. Release 1 delivers calls for *M. zebra* and *A. calliptera* only; that is project 3's step, and it is settled in the release (`deliver`), not here.

## 4. Outputs

Everything is under `get_layout().results("TreeFam2026")` and `results("Compara2026")` (a `ResultsLayout`). Release 1 (coordinator v0.9.0) moves these to `Results/v1/TreeFam/` and `Results/v1/Compara/`; use the API, not the names, and the move will be a one-line change.

| Path | Content |
|---|---|
| `assign/<prefix>.assign.tsv` | one per placed species, the existing assignment schema (`best_family` is the TreeFam or Compara family; the score, coverage and fusion columns are yours) |
| `families/<family>/` | one folder per family, the family contract below |
| `run.json` | the run record: `method` (`TreeFam` or `Compara`), `method_version`, `parameters` (search thresholds, aligner, placement program and their versions), `started`, `finished`, `coordinator_version` (`cichlid_data.__version__`), `catalog_md5` (`Catalog.load(layout).md5()`), `species` (the placed prefixes). Write it with `schemas.write_run`. From release 1 it also has a `release` field |

### The family contract (summary; the full text is in `coordinator_spec.md`)

For every family, a folder named by the family ID (`TF101001`, `C116_012345`) holding:

- **`tree.nwk`:** the final tree: the reference tree with your genes placed on it. Leaf names are unique: our gene IDs (`<prefix>_<id>`) for the placed proteins, and the reference proteins' own IDs for the rest.
- **`combined.aln`:** the aligned FASTA, with a row for every leaf of the tree, rows named like the leaves.
- **`done.json`:** a JSON object with the inputs hash that lets a rerun skip finished work. Write it **last**: a folder without it is unfinished. The `"rooted"` field says whether the tree is already rooted; absent means unrooted, and project 3 midpoint-roots it. **To confirm with project 3:** TreeFam and Compara reference trees come rooted by their own species trees; if the placed tree keeps that root, write `"rooted": true` and `"rooting": "reference tree"`, but only if the root of `tree.nwk` really has two children (validation checks it).
- The files you keep beside these (`queries.fa`, `ref.tree`, `epa_result.jplace`, `placements.tsv`, ...) are yours and are not checked. The residue files (`residue_map.tsv`, `residues.tsv.gz`, `residues.done.json`) belong to project 4: leave them as they are, do not extend them, and do not rely on them.

Existing results already follow this: `cichlid-data validate` passes on every family folder of `Compara2026` and `Rebuild2026` and on every assignment file in them (18 files, 9 species each).

## 5. Validation and handing over

```bash
cichlid-data check --method TreeFam2026 --fetch                        # before a run (and again for Compara2026)
R=$(cichlid-data where results TreeFam2026)                            # never type the data root
cichlid-data validate "$R"/assign/*.assign.tsv                         # one per placed species
cichlid-data validate "$R"/families                                    # every family folder, leaf names checked against the catalog
cichlid-data validate "$R"/run.json
cichlid-data push 'results/TreeFam2026/*' --dry-run                    # then without --dry-run
```

Hand over to project 3 only when all of these pass, for both methods.

## 6. Replacing paths and transfers: verify on real data

When this repository switches from its own paths, transfers and formats to the coordinator's, the results it produces must be **identical** to what it produced before. The baselines are `Rebuild2026` (TreeFam) and `Compara2026` (Compara) in the data folder, kept for this purpose. They cover nine species (*Acal, Acen, Hbur, Mzeb, Nbri, Onil, Pmar, Pmic, Sdia*). *Etroplus*, betta and stickleback have no baseline: their results are new, so check them against the validators and by eye. Compare the assignment tables and family folders for the nine species file by file; explain every difference, in a note kept with the repository, before the change is accepted. Do not delete the baselines until this is done.

## 7. What this repository must not do

- **No orthology calls** and no combining of methods; those are project 3's.
- **No residue, variant or genotype work;** that is project 4's.
- **No hand-built paths,** no copy of the species list, no rclone, no downloads, no editing of the proteomes, the catalog or the TreeFam and Compara inputs.
- **No change to the TreeFam 9 or Compara reference sets** to make a run work. If they need changing, that is a decision recorded in the coordinator (`config/compara_panel.tsv`, and the input versions of a release).
- **No change to a file format.** If the contract is not enough, ask the coordinator; a schema change is a version change.

## 8. Acceptance criteria

1. The tests run offline and pass.
2. `cichlid-data check` passes for both methods, and every run is validated as in section 5.
3. `run.json` records the version and parameters of every program used, and the coordinator version and catalog checksum, and passes `cichlid-data validate`.
4. The nine baseline species reproduce the baseline results exactly, or every difference is explained (section 6).
5. Every assignment file and every family folder passes `cichlid-data validate`, for the 12 placed species.
6. The results are pushed, and `cichlid-data status 'results/TreeFam2026/*'` (and Compara2026) shows everything in sync.
7. No other path, species list, or cloud call, and no calling or residue code, appears anywhere in the repository.

# OrthoFinder repository: instruction sheet

*For the `cichlid-orthofinder` repository. October 2026. Written against cichlid-coordinator v0.9.0. Read first: `project_manifest.md` (why the project exists), `species_final.md` (the species), and `coordinator_spec.md` (the layout, the schemas and the family contract, which this sheet relies on).*

## 1. Role

This repository runs OrthoFinder on the project's 35 proteomes and turns its gene trees into the project's standard family format. It owns:

1. **The OrthoFinder run:** the settings in section 4, run reproducibly.
2. **The conversion:** OrthoFinder's orthogroups and gene trees into `Results/v1/OrthoFinder/assign/` and `families/`, following the family contract.
3. **OrthoFinder's own outputs,** kept as written in `Results/v1/OrthoFinder/native/`, so they are backed up and can be pulled by anyone.

It does **not** make orthology calls, name genes, find human or zebrafish orthologs, or analyze results. The homology decider (project 3, `cichlid-homology`) reads the families this repository writes and calls them, with the same rule it uses for TreeFam and Compara, so the three methods' calls are comparable: differences between them then reflect the trees, not the calling logic.

**Rules (the same as for every repository):**

- Never build a data path by hand. Ask the coordinator (`get_layout()`, `cichlid-data where ...`).
- Never call rclone. Moving data to and from the cloud is the coordinator's (`get_store()`, `cichlid-data pull | push | ensure`).
- Never download an input. The proteomes come from the catalog.
- Never copy the species list or the species tree into this repository. Read them from the catalog.

Depend on the coordinator with a range pinned to one minor version, for example `cichlid-coordinator>=0.9,<0.10`, and use only its public API (`cichlid_data`) and its command line. Below 1.0.0 a minor version may change the layout API (v0.9.0 did: results are now `Results/<release>/<part>/`).

## 2. Repository layout (suggested)

```
cichlid-orthofinder/
├── pyproject.toml           package cichlid_orthofinder; depends on cichlid-coordinator>=0.9,<0.10
├── src/cichlid_orthofinder/
│   ├── run.py               prepare the input folder, run OrthoFinder, record the settings
│   ├── convert.py           orthogroups and gene trees -> assign/ and families/
│   ├── native.py            what is kept under native/
│   └── cli.py               the cichlid-orthofinder command
├── tests/                   offline: small fake OrthoFinder output folders
└── docs/                    this sheet, notes on the runs
```

Moving code in from elsewhere follows the same rule as in the coordinator: verbatim, with its tests, changing only paths, transfers and formats.

## 3. Inputs (all through the coordinator)

| Input | How to get it |
|---|---|
| The species | `Catalog.load(get_layout()).species(release="v1", part="OrthoFinder")`: the 35 species, in the catalog's order |
| The proteomes | `get_layout().proteome(prefix)`: one protein per gene, headers `>gene_id`, already valid. Make sure they are local first: `cichlid-data check --release v1 --part OrthoFinder --fetch` (which also verifies the species tree) |
| The species tree | `Catalog.load(layout).species_tree("v1", "OrthoFinder")`, or `cichlid-data catalog tree --release v1 --part OrthoFinder`: rooted, fully resolved, leaves = the species prefixes |

Run `cichlid-data check --release v1 --part OrthoFinder --fetch` before every run; stop if it fails.

**Names matter.** OrthoFinder takes a species' name from its proteome file name, and `-s` needs the species tree's leaves to be those names. So the input folder OrthoFinder reads holds one file per species named `<prefix>.faa` (links or copies of the coordinator's proteomes; the file extension OrthoFinder accepts is *to confirm* against the installed version). Do not rename a species anywhere.

**Gene IDs.** Our gene IDs (`Mzeb_101467881`) contain an underscore, which is fine for the coordinator. OrthoFinder renames every sequence internally (`<species number>_<sequence number>`) and writes the mapping back to our names in its working directory (`SequenceIDs.txt` and `SpeciesIDs.txt`; names *to confirm* against the installed version). The conversion in section 5 maps every name back with it.

## 4. The OrthoFinder runs

Agreed settings:

| Setting | Value |
|---|---|
| Gene trees | multiple sequence alignments and a maximum-likelihood tree for every orthogroup: `-M msa`, with IQ-TREE (`-T iqtree`; option names *to confirm* against the installed OrthoFinder) |
| Species tree | our species tree, with `-s <species tree>` (the file from the catalog; do not edit it) |
| Second run | the same, **without** `-s`, so OrthoFinder infers its own species tree; kept for comparison only |
| Everything else | OrthoFinder's defaults, unless agreed otherwise. Thread counts are the machine's business and are not science |

Every setting actually used goes into `run.json` (`parameters`), together with the OrthoFinder version (`method_version`), the IQ-TREE and alignment-program versions, and the command line. A run that can't say how it was made is not accepted.

**Where the two runs go.** The first run, with our species tree, is the one that becomes `Results/v1/OrthoFinder/`: its families are converted and called. The second run, with OrthoFinder's own tree, is a comparison: keep its OrthoFinder output under `Results/OrthoFinder2026/native/<a name for it>/` (for example `native/own_species_tree/`), and do not convert it. *(This arrangement is an assumption to confirm: if calls are wanted from the second run it becomes a run of its own, with its own name and catalog entry.)* Each top-level folder under `native/` is stored as a separate tar, so adding the second run does not repack the first.

## 5. Outputs

Everything below is under `get_layout().results("v1", "OrthoFinder")` (a `ResultsLayout`; `v1` is the release, which is the coordinator's to name, and `OrthoFinder` the part):

| Path | Content |
|---|---|
| `assign/<prefix>.assign.tsv` | one per species, the existing assignment schema: `gene_id` and `protein_length` for every gene; `best_family` is the gene's orthogroup, **empty for a gene in none**; the score, coverage and fusion columns **empty** (they belong to placement methods) |
| `families/<family>/` | one folder per orthogroup, the family contract below |
| `native/` | OrthoFinder's own result folders, as it wrote them (see below) |
| `run.json` | the run record: `release` (`v1`; it must match the folder it is in), `method` (`OrthoFinder`), `method_version`, `parameters`, `started`, `finished`, `coordinator_version`, `catalog_md5`, `species` (the 35 prefixes). Write it with `schemas.write_run`, taking `coordinator_version` from `cichlid_data.__version__` and `catalog_md5` from `Catalog.load(layout).md5()` (never compute either by hand) |

### The family contract (summary; the full text is in `coordinator_spec.md`)

For every **orthogroup** (OrthoFinder's primary orthogroups, `OG0001234`; they have gene trees), a folder named by its family ID (letters, digits, `.`, `_`, `-`) holding:

- **`tree.nwk`**: OrthoFinder's gene tree with **every leaf renamed back to our gene ID** (`<prefix>_<id>`). All leaves are ours: every sequence comes from the 35 proteomes, so there are no reference-protein IDs. A leaf that still looks like `0_123` is refused by validation. Keep OrthoFinder's branch lengths and supports.
- **`combined.aln`**: required. The orthogroup's alignment (aligned FASTA), rows renamed to the same gene IDs, with a row for every leaf. (The name is the placement methods'; they write other files beside it, which OrthoFinder families do not have.)
- **`done.json`**: a JSON object with the inputs hash that lets a rerun skip finished families, `"release": "v1"`, and the rooting. OrthoFinder roots its gene trees using the species tree, so write `"rooted": true` and `"rooting": "orthofinder species tree"`, **but only if the root of that `tree.nwk` really has two children** (validation checks it). If a tree arrives unrooted, write `"rooted": false` and leave the rooting to the homology decider, which midpoint-roots it. Never write `"rooted": true` to mean "probably".

The families are the **orthogroups (`OG…`)**; the hierarchical orthogroups (`N0.HOG…`) are not converted and stay in `native/` with the rest of OrthoFinder's output. An orthogroup for which OrthoFinder wrote no gene tree gets no family folder; its genes still appear in `assign/` with their orthogroup, and in OrthoFinder's own tables. Report how many there are in the conversion's summary, so the gap is visible.

### What is kept under `native/`

OrthoFinder's result folder for each run, so that anyone can look at what it actually produced: the orthogroup tables, the gene trees, the alignments, the species tree, the duplication and ortholog tables, and the logs. Its bulk working files (sequence-search result files for every species pair, for example) may be tens of gigabytes: **decide with the coordinator what is kept before pushing the first run**, from the size of a real run. Never push `native/` before that.

## 6. Conversion (steps)

1. Read the mapping files for the run; check every species and sequence is accounted for (35 species; the count of sequences equals the count of proteins in the proteomes).
2. For every species, write `assign/<prefix>.assign.tsv` from the orthogroup table through `schemas.ASSIGN.write` (it refuses invalid rows). Every gene of the proteome appears, in the proteome's order.
3. For every orthogroup with a tree, write the family folder: renamed tree, renamed alignment, `done.json`. Write `done.json` **last**: a folder without it is unfinished, and a rerun redoes it.
4. Write `run.json` with `schemas.write_run`.
5. Validate (section 7). Then, and only then, push.

Conversion is deterministic: the same OrthoFinder output gives the same files, byte for byte.

## 7. Validation and handing over

```bash
cichlid-data check --release v1 --part OrthoFinder --fetch           # before a run: 35 proteomes and the species tree
R=$(cichlid-data where results v1 OrthoFinder)                        # never type the data root
cichlid-data validate "$R"/assign/*.assign.tsv                         # one per species
cichlid-data validate "$R"/families                                    # every family folder, leaf names checked against the catalog
cichlid-data validate "$R"/run.json
cichlid-data push 'results/v1/OrthoFinder/*' --dry-run                # then without --dry-run
```

Hand over to the homology decider only when all of these pass. The tests of this repository use small fake OrthoFinder output folders and run offline; they must cover the renaming (including a leaf that cannot be mapped), a gene in no orthogroup, rooted and unrooted trees, and the byte-for-byte repeatability of the conversion.

## 8. What this repository must not do

- **No orthology calls.** No human or zebrafish relationships, no "one-to-one", no names. OrthoFinder's own ortholog tables stay in `native/`; they are not project results.
- **No hand-built paths,** and no copy of the species list or the species tree.
- **No rclone,** no downloads, no editing of the proteomes or of the catalog.
- **No change to the species tree** to make a run work. If the tree needs changing, that is a decision recorded in the coordinator (`species_final.md` and `config/species_tree_orthofinder.nwk`, together).
- **No change to a file format.** If the contract is not enough, ask the coordinator; a schema change is a version change.

## 9. Acceptance criteria

1. The tests run offline and pass.
2. `cichlid-data check --release v1 --part OrthoFinder` passes, and the species tree used for `-s` is the catalog's, byte for byte.
3. `run.json` records the OrthoFinder version, every setting, and the coordinator version and catalog checksum, and passes `cichlid-data validate`.
4. Every `assign` file and every family folder passes `cichlid-data validate`, with no leaf left in OrthoFinder's internal naming.
5. Every gene of every proteome appears in its species' assignment file.
6. The first run's results and its native outputs are pushed, and `cichlid-data status 'results/v1/OrthoFinder/*'` shows everything in sync.
7. No other path, species list or cloud call appears anywhere in the code.

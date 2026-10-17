# Coordinator repository: instruction sheet

*For the Central Project Coordinator workspace and repository (`cichlid-coordinator`). October 2026. Current at v0.8.0.*

## 1. Role

The coordinator is the one package every other repository depends on. It owns:

1. **File locations:** one layout, identical relative paths on every machine and in the cloud.
2. **Transfers:** all movement between local machines and the cloud (Dropbox, via rclone).
3. **Inputs:** every download, the shared proteomes (one protein per gene), the species catalog, and provenance records.
4. **File formats** at the boundaries between repositories: readers and writers that check them.

It does **not** run tree methods, make orthology calls, analyze results or serve the browser. Those live in the method and analysis repositories, which use the coordinator through its public API only.

**Rule for every other repository:** never build a data path by hand, never call rclone directly, never download an input yourself. Ask the coordinator.

**The repositories that use it.** Two, and the boundary between them is the coordinator's (layout, schemas, catalog), so neither can drift from it:

- `cichlid-trees`: TreeFam and Compara, plus the shared core: alignment and placement for the placement methods, and orthology calling and residues for **any** method's results, so every method's calls come from one rule.
- `cichlid-orthofinder`: runs OrthoFinder on the catalog's 35 proteomes, converts its gene trees into the standard family format in `Results/OrthoFinder2026/`, and keeps OrthoFinder's own outputs alongside (`native/`). It never makes orthology calls; `cichlid-trees` calls its families. Its instruction sheet is `docs/orthofinder_repo_spec.md`.

## 2. Repository layout

```
cichlid-coordinator/
├── pyproject.toml           package cichlid_data; console script cichlid-data
├── config/
│   ├── species.tsv          project species (edited by hand)
│   ├── species_sets.toml    which species each method uses (and the species tree a method must use)
│   ├── species_tree_orthofinder.nwk   rooted, fully resolved tree of the 35 OrthoFinder species (for -s)
│   ├── compara_panel.tsv    Compara reference panel (role per Ensembl species)
│   └── sources.toml         every external resource: path, URL/accession, version, licence, citation
├── src/cichlid_data/
│   ├── __init__.py          public API
│   ├── config.py            Settings
│   ├── layout.py            Layout, ResultsLayout
│   ├── units.py             DataUnit registry
│   ├── store.py             Store, RcloneBackend, LocalBackend
│   ├── catalog.py           Catalog, Species
│   ├── provenance.py        SourceRecord, checksums
│   ├── schemas.py           boundary file formats: readers, writers, validators
│   ├── proteomes.py         one protein per gene
│   ├── downloaders/         ncbi, local, annotations, expression, h5ad, rnaseq, blumer, blumer_samples, compara
│   ├── newick.py            Newick parser, shared with the pipeline
│   ├── backup.py            push-what-you-wrote helper for the moved scripts
│   ├── check.py             cichlid-data check
│   ├── migrate.py           one-time move from the old layout
│   └── cli.py               the cichlid-data command
├── tests/                   all tests run offline against LocalBackend
└── docs/                    project manifest, species lists, this sheet, the OrthoFinder repository sheet
```

## 3. Data layout (the paths Layout must produce)

```
<root>/
├── Inputs/
│   ├── Proteins/
│   │   ├── catalog.json
│   │   ├── downloads/<prefix>/        genomic.gff, protein.faa, sequence_report.jsonl, source.json
│   │   └── proteomes/<prefix>.faa     + <prefix>.genes.tsv
│   ├── Families/
│   │   ├── TreeFam9/                  hmm/, trees/<family>.aln|.nhx, protein_to_species.txt, source.json
│   │   └── Compara116/                dumps/, prepared/ (trees/, search_db.fa, members.tsv,
│   │                                  protein_to_species.txt, family_audit.tsv), source.json
│   ├── Annotations/                   GO/, PANTHER/, Ensembl_names/, Legacy2018/   (each with source.json)
│   ├── Expression/                    HPA/, Zebrahub/, Daniocell/, Bgee/, Brawand/, raw/   (each with source.json)
│   └── Populations/Blumer2025/        VCFs, tables, reference/
├── Results/
│   ├── <Method><Year>/                TreeFam2026, Compara2026, OrthoFinder2026: identical structure
│   │   ├── assign/<prefix>.assign.tsv
│   │   ├── families/<family>/         alignment.fa, tree.nwk, done.json (the family contract, section 4;
│   │   │                              placement methods add their placements)
│   │   ├── names/
│   │   ├── native/                    the method's own outputs (OrthoFinder's result folders), as it wrote them
│   │   ├── orthology/<prefix>.orthologs.tsv
│   │   ├── residues/
│   │   └── run.json
│   └── Combined2026/                  orthology/, site/
├── Analyses/<Name>/                   Variants, Tissues, Duplications, Enrichment, ...
└── Archive/Project2018/               the 2018 resource's own files, kept as a record
```

## 4. Public API

Everything here is importable from `cichlid_data`. Anything else is internal and may change.

### Settings (`config.py`)

Read from `~/.config/cichlid-data/config.toml`, overridable by environment variables (`CICHLID_DATA_ROOT`, `CICHLID_DATA_REMOTE`, `CICHLID_DATA_RCLONE`):

```toml
local_root   = "~/Temp/TreeFamData"
cloud_remote = "ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/TreeFamData"
rclone       = "rclone"
```

A clear error if neither the file nor the variables exist, saying exactly what to create.

### Layout (`layout.py`)

`get_layout()` returns a `Layout`. Every method returns a path relative to the root, resolvable locally (`.local`) or in the cloud (`.cloud`):

| Method | Path |
|---|---|
| `proteome(prefix)`, `genes_table(prefix)` | `Inputs/Proteins/proteomes/<prefix>.faa`, `.genes.tsv` |
| `download(prefix)` | `Inputs/Proteins/downloads/<prefix>/` |
| `catalog()` | `Inputs/Proteins/catalog.json` |
| `family_resource(name)` | `Inputs/Families/<name>/` |
| `annotations(name)`, `expression(name)`, `population(name)` | `Inputs/Annotations/<name>/` and so on |
| `results(run)` | a `ResultsLayout` for `Results/<run>/` |
| `analysis(name)` | `Analyses/<name>/` |

`ResultsLayout` methods: `assign(prefix)`, `family(family)`, `orthology(prefix)`, `names()`, `native()`, `residues(prefix)`, `run_record()`, and `unit(part)` (the `DataUnit` for a part).

### DataUnit (`units.py`)

A named, transferable piece of the tree, with a policy:

| Policy | Meaning | Examples |
|---|---|---|
| `immutable` | Existing files never change. File-by-file units may gain new files (e.g. a `source.json`); tar-packed units can't change at all | downloads, TreeFam 9, local species sources |
| `regenerable` | Pipeline output; may be replaced by a newer run | proteomes, Results/*, Analyses/* |
| `on_demand` | Large; never pulled as part of a bulk pull, only by `ensure` | population VCFs, Daniocell raw data |

`units.all()`, `units.get(name)`, `units.for_path(path)`.

Every results run has the same units: `results/<run>/assign`, `families`, `names`, `native`, `orthology`, `residues`, `site` and `run`. `families` is sharded by family ID (the last three characters dropped, about 1,000 families per shard: `TF101001` in `TF101`, `OG0001234` in `OG0001`, `N0.HOG0001234` in `N0.HOG0001`). `native` holds a method's own outputs, **regenerable and tar-packed, one tar per top-level folder** (so a second run kept beside the first uploads only its own tar); it is backed up and pulled like everything else. Neither is ever packed by hand.

### Store (`store.py`)

`get_store()` returns a `Store` over a `Backend` (`RcloneBackend` in use, `LocalBackend` in tests):

| Method | Behavior |
|---|---|
| `ensure(path)` | Pull the file or folder if missing or different locally; no-op otherwise. Returns the local path |
| `pull(unit)` | Bring a unit local, checksums verified |
| `push(unit)` | Send a unit to the cloud. Refuses for immutable units that differ; takes a cloud lock first |
| `status(unit=None)` | Per unit: in sync, local newer, cloud newer, or missing |
| `dry_run=True` | Every method prints what it would do, without doing it |

Every transfer is logged (`<root>/.cichlid-data/transfers.log`: time, user, host, unit, direction, files, bytes).

**Locks:** `push` writes `<unit>/.lock` in the cloud (user, host, time) and removes it after; a lock older than 6 hours may be broken with `--force`.

### Catalog and Species (`catalog.py`)

Three configuration files, each with one job:

- `config/species.tsv`: what Stage 1 downloads (the download manifest).
- `config/species_sets.toml`: every organism the project uses (prefix, name, role, lineage, source; roles `project`, `relative`, `reference`, `outgroup`, `candidate`) and each method's species: `queries = "project"` for TreeFam and Compara, an explicit list for OrthoFinder.
- `config/compara_panel.tsv`: Compara's reference panel.

`Catalog.build(layout)` combines them with what is on disk into `catalog.json`; `Catalog.load(layout)` reads it.

`Species` fields: `prefix`, `species`, `common`, `role`, `lineage`, `source`, `accession`, `version`, `retrieved`, `in_manifest`, `status` (`proteome ready`, `downloaded`, `local source only`, `not downloaded`), `proteome`, `genes_table`, `proteome_md5`, `counts` (proteins, genes), `used_by` (method names), `notes`.

Queries: `catalog.species(method=..., role=..., status=...)`, `catalog.missing(method)` (organisms whose proteome is not ready), `catalog[prefix]`. `validate()` returns errors (bad prefixes, unknown roles, method members without an entry or still candidates) and warnings (project species not in `species.tsv`).

**Species tree.** A method may name a species tree in `species_sets.toml` (`species_tree = "species_tree_orthofinder.nwk"`, a file in `config/`). It must parse, be rooted and fully resolved (two children at the root and at every node), and have exactly the method's species as leaves. `catalog validate` and `catalog build` check it (an invalid tree stops the catalog from being written, so the tree and the species list cannot drift apart), and so does `cichlid-data check` for that method. `catalog.species_tree(method)` returns its path (`KeyError` if the method names none); its checksum is recorded with the other configuration files in `catalog.json`. The OrthoFinder tree is the one in `species_final.md`; a test keeps the two identical. Adding a species to OrthoFinder therefore means changing the list **and** the tree.

Commands: `cichlid-data catalog build | show | missing | validate | tree`, with `--method` and `--role` (`catalog tree --method M` prints the species tree's path, for scripts).

### SourceRecord (`provenance.py`)

Configuration files are TOML (read by Python's standard library, so the package has no dependencies).

`config/sources.toml` describes every external resource (path, name, type, version, URL or accession, licence, citation; unknowns marked "to confirm"). `cichlid-data provenance write` turns each entry into a `source.json` beside the data, with an inventory of its files (sizes, and MD5 checksums for files up to 1 GB; larger files by size). Species downloads are described by Stage 1's own `download.json` (accession, NCBI metadata, checksums, date), which `read_record()` reads as a `SourceRecord`. `provenance check` verifies every record (sizes; `--deep` also checksums), and reports input folders with no description.

The record:

```json
{"name": "...", "type": "ncbi_refseq | ensembl | treefam | authors | url",
 "accession": "...", "version": "...", "url": "...", "retrieved": "2026-10-15T12:00:00Z",
 "licence": "...", "citation": "...", "files": {"protein.faa": {"md5": "...", "bytes": 0}},
 "retrieved_by": "cichlid-data 1.2.0"}
```

### Schemas (`schemas.py`)

One module per boundary file, each with a `read()`, a `write()` and a `validate()` that checks columns, types and ID formats. Files are tab-separated with a header row.

**Gene IDs** everywhere: `<prefix>_<id>`, where `<prefix>` is letters and digits only and `<id>` uses letters, digits, `.` and `-` only (no underscores; the first underscore separates the prefix).

| File | Required columns (in order) |
|---|---|
| Proteome FASTA `<prefix>.faa` | One record per gene; header `>gene_id description`; standard amino acids, `X`; no internal `*` |
| Genes table `<prefix>.genes.tsv` | `gene_id, ncbi_gene_id, annotation_gene_key, symbol, description, biotype, seqid, start, end, strand, mitochondrial, n_transcripts, n_proteins, chosen_protein, chosen_transcript, protein_length, partial, low_quality, residues_replaced, internal_stops, product` |
| Assignment `<prefix>.assign.tsv` (`best_family` and `second_family` are family IDs; for OrthoFinder the orthogroup, with the score columns empty) | `gene_id, protein_length, best_family, score, evalue, protein_coverage, family_coverage, protein_region, n_families, second_family, second_score, fusion_partner, fusion_partner_score, fusion_partner_region, fusion_partner_family_coverage, flags` |
| Family folder `families/<family>/` | `tree.nwk` (final tree, leaf names = gene IDs or reference protein IDs), `alignment.fa`, `done.json` (inputs hash, rooting): see "The family contract" below |
| Orthology calls `<prefix>.orthologs.tsv` | Gene columns `gene_id, ncbi_gene_id, ncbi_symbol, ncbi_description, seqid, start, end, strand, mitochondrial, protein, protein_length, treefam_family, tree_mode, placement_status, placement_lwr`; then, for each reference (`human`, `mouse`, `zebrafish`, in that order), the block `<ref>_relationship, <ref>_ortholog_type, <ref>_symbol, <ref>_official_id, <ref>_ensembl_gene, <ref>_treefam_protein, <ref>_description, <ref>_name_status, <ref>_ref_in_subtree, <ref>_own_species_in_subtree, <ref>_median_cichlid_copies, <ref>_distance, <ref>_call_support`; then `human_uniprot, ncbi_symbol_agrees, confidence, caveats` |
| Run record `run.json` | `method, method_version, parameters, started, finished, coordinator_version, catalog_md5, species` |

Allowed values: `relationship` in {`Ortholog`, `Paralog`, `NoHomolog`}; `ortholog_type` in {`OneToOne`, `OneToTwo`, `OneToMany`, empty}; `confidence` in {`high`, `medium`, `low`}.

Booleans are written `True`/`False`; empty strings mean missing. Duplicate gene IDs within a species get a suffix (`Mzeb_101467881.2`), which the ID rules allow. The proteome check also cross-checks the genes table beside it: every gene with a chosen protein has exactly one protein, and vice versa. The orthology schema accepts any set of references, in blocks, so a method calling orthologs to other species stays valid. Stage 7's combined calls are a comparison output and are not yet specified (the comparison is being redesigned).

The schemas were checked against files written by the pipeline's own Stage 2, 3 and 6 code: all valid.

`cichlid-data validate PATH...` checks files (the kind is inferred from the name, or given with `--kind`: `genes`, `assign`, `orthologs`, `proteome`, `run`, `family`, `families`). A folder named `families` is checked family by family; given a catalog it also checks that leaf names starting with a species prefix are well-formed gene IDs.

#### The family contract

Every method writes `Results/<run>/` in the same way, whatever built the trees, so the calling code can be the same for all of them. For any run:

- **`assign/<prefix>.assign.tsv`**: the existing schema. For OrthoFinder, `best_family` is the orthogroup (empty for a gene in none) and the score columns may be empty.
- **`families/<family>/`**, one folder per family:
  - **Family ID**: letters, digits, `.`, `_` and `-`, starting with a letter or digit (`TF101001`, `C116_012345`, `OG0001234`, `N0.HOG0001234`). The ID is the folder name.
  - **`tree.nwk`**: the final tree. Leaf names are unique and are **our gene IDs** (`<prefix>_<id>`) or **reference-protein IDs**. OrthoFinder renames sequences internally (`0_123`); its trees must be mapped back before they are written, and a leaf that looks like an internal name is refused. A leaf that starts with a catalog prefix must be a well-formed gene ID.
  - **`alignment.fa`**: the aligned sequences (FASTA, rows of equal length), with a row for every leaf of the tree. Optional for validation; rows without a leaf are allowed. *(File name to be confirmed against the TreeFam2026 output before v0.8.0 is tagged.)*
  - **`done.json`**: a JSON object. It holds the inputs hash that lets a rerun skip finished work, and says whether the tree is already rooted: `"rooted": true` with a `"rooting"` text (for example `"rooting": "orthofinder species tree"`), or `"rooted": false`. `rooted` must be a true or false value; when true, `rooting` is required and the root of `tree.nwk` must have two children. **Consumers use an existing root and midpoint-root only otherwise**; a family that says nothing, or says something malformed, is treated as unrooted. `schemas.family_rooting(folder)` returns `(rooted, how)` so no consumer reimplements this rule.
- **`run.json`**: the run record, as above.

`schemas.validate_family(folder, prefixes=None)` and `schemas.validate_families(folder, prefixes=None)` check all of this; `cichlid-data validate Results/<run>/families` is the command.

A schema change is a version change of the package (Section 8).

### Proteomes (`proteomes.py`)

From each species' download, one protein per protein-coding gene: the pipeline's Stage 2, with its parsing and selection code moved over unchanged (the longest complete isoform, else the longest; NCBI's low-quality models last; mitochondrial genes flagged; residues outside the standard alphabet replaced and counted). Genes on alternate loci and patches (sequence roles `alt-scaffold`, `fix-patch`, `novel-patch`; e.g. human GRCh38) are left out: they are second copies that would look like recent duplications, and are counted in the QC table instead (a change from Stage 2). Roles are used rather than assembly-unit names, which differ between assemblies. If more than half of a species' genes would be skipped, the build stops instead. Second change from Stage 2: each protein is linked to every gene whose CDS names it (Stage 2 kept only the last such gene, so a protein annotated on a chromosome and on a patch was lost from the chromosome gene, about 870 human genes); a protein still shared by two loci is kept once, at the first, and the others counted (`shared_protein_loci`). Outputs are written through the schemas, and a species whose proteome fails validation leaves nothing behind. Species whose proteome is newer than their download are skipped unless `--force`; a missing download is fetched from the cloud first. Afterwards the catalog is rebuilt, and `--push` uploads the proteomes. QC rows for each species are kept in `stage2_summary.tsv`, one row per species.

`cichlid-data proteomes [--species P,...] [--select complete-longest|longest] [--force] [--push]`

### Downloaders (`downloaders/`)

The pipeline's download and preparation scripts, moved over with their logic unchanged: only their default folders (now from the layout) and their backups (now pushes through the store) changed. Each keeps its own options.

| Module | Command | Sources |
|---|---|---|
| `ncbi.py` | `download genomes` | Stage 1: genomes and annotations from NCBI Datasets (`source=ncbi`; verified, resumable, checked against NCBI's metadata) or local files (`source=local`; paths may be relative to the data root). Species already in the cloud are pulled rather than fetched again; the catalog is rebuilt afterwards |
| `local.py` | `prepare-local` | Converts a non-NCBI annotation into the download files (`prepare_local_species.py`) |
| `annotations.py` | `download annotations` | GO (go-basic.obo, human and ZFIN GAFs) and PANTHER classifications |
| `expression.py` | `download expression --source hpa\|bgee\|daniocell` | Human Protein Atlas, Bgee, Daniocell raw files |
| `h5ad.py` | `summarize h5ad` | Single-cell atlases (Zebrahub, Daniocell) to per-cell-type summaries; needs the `[expression]` extra (NumPy, SciPy) |
| `rnaseq.py` | `summarize rnaseq` | Tissue RNA-seq (nf-core outputs, e.g. Brawand) to per-tissue summaries |
| `blumer.py` | `download populations` | Blumer et al. 2025 VCFs and tables (about 200 GB; resumable, checksum-verified) |
| `blumer_samples.py` | `populations samples` | Sample metadata for the VCFs (Table S1, optionally the lab sheet; `[samples]` extra for gspread) |
| `compara.py` | `prepare compara` | Prunes the Compara dumps to `config/compara_panel.tsv`; the dumps, panel and output folder default to their layout paths |

`newick.py` (the pipeline's Newick parser) is a shared module of the coordinator, so the pipeline can import it from here.

### Check (`check.py`)

`cichlid-data check [--method M] [--fetch] [--deep]` confirms, for each method, that everything it needs is present and sound: the proteomes of its species (valid against the schemas, with the genes-table cross-check), its resources (TreeFam 9 or the prepared Compara trees, and the legacy name tables), each verified against its source record (sizes; `--deep` also checksums), and, for a method that names one, its species tree. `--fetch` pulls anything missing from the cloud first. It exits with an error if anything is missing or invalid, so run scripts can call it before starting.

## 5. Command line

```
cichlid-data init [--root DIR] [--remote REMOTE] [--force]   once per machine
cichlid-data where <thing> [args] [--cloud]          print a path, e.g. "where proteome Mzeb", "where results OrthoFinder2026 native"
cichlid-data env                                     shell exports of the main paths
cichlid-data units | status [pattern]
cichlid-data pull | push <pattern>... [--force] [--dry-run]
cichlid-data ensure <relative path>                  fetch one file or shard if missing
cichlid-data retire <unit key> [--dry-run]           delete a unit's or shard's cloud copy
cichlid-data migrate [--apply] [--cleanup]           the one-time move to the new layout
cichlid-data provenance write [--overwrite] [--md5-limit-gb N] | check [--deep] | show <folder>
cichlid-data catalog build | show | missing | validate | tree [--method M] [--role R]
cichlid-data validate <file or folder>... [--kind K] [--prefix P]
cichlid-data download genomes [--species P,...] [--check-only] [--force] [--keep-zip] [--no-fetch] [--push]
cichlid-data prepare-local --prefix P ...
cichlid-data proteomes [--species P,...] [--select complete-longest|longest] [--force] [--no-fetch] [--push]
cichlid-data download annotations | expression | populations ...
cichlid-data populations samples ...
cichlid-data prepare compara ...
cichlid-data summarize h5ad | rnaseq ...
cichlid-data check [--method M] [--fetch] [--deep]
```

The moved scripts' commands take the scripts' own options (`--help` lists them).

## 6. Migration from the current layout

`migrate.py` moves the project's **inputs** to the new layout, on the local machine and on the cloud by server-side moves (no re-upload), printing every move in a dry run first. **Results are not moved:** the TreeFam and Compara runs will be redone with the new species set (adding species changes nearly every family, so little would be reused), and the comparisons and analyses will be redone with new methods. The old results stay where they are, for comparison.

| Current | New |
|---|---|
| `Rebuild2026/genomes/<p>/` | `Inputs/Proteins/downloads/<p>/` (repacked: one shard per species) |
| `Rebuild2026/proteomes/` | `Inputs/Proteins/proteomes/` |
| `Data/hmm_lib/`, `Data/OriginalTreeFiles/`, `Data/protein_to_species.txt` | `Inputs/Families/TreeFam9/hmm/`, `trees/`, `protein_to_species.txt` |
| `Data/Compara116/` | `Inputs/Families/Compara116/` (`dumps/`, `prepared/`) |
| `Data/Annotations/` | `Inputs/Annotations/` |
| `Data/Expression/` | `Inputs/Expression/` |
| `Data/Blumer2025/` (VCFs, tables, reference) | `Inputs/Populations/Blumer2025/` |
| `Data/Esur/` | `Inputs/Proteins/local_sources/Esur/` |
| `Data/Proteins/`, `Data/Species/`, `Data/AllCichlidProteins.*`, `Data/NewTreeFiles/` (the 2018 resource's own files) | `Archive/Project2018/` (a record; not used by the new pipeline) |

**Left in place by design:** `Rebuild2026/`, `Compara2026/` and `Combined2026/` (phylogenies, calls, browser site, analyses), and `Data/Blumer2025/coding/` and `inversions/` (variant analyses).

Folders emptied by the moves are removed, so `Data/` ends up holding only the variant analyses. The dry run ends with two lists: what is left in place by design, and anything unexpected in `Data/`, so nothing is left behind silently. The proteomes' cloud index entry and sync state are carried over, so they are in sync immediately.

## 7. Acceptance criteria

1. Every public method and CLI command has tests, all runnable offline against `LocalBackend`.
2. The existing tests of the moved code (Stages 1–2, local species preparation, sync, downloaders) pass in the new repository.
3. `cichlid-data check` passes on `mzebra` after migration, and the pipeline (updated to the coordinator's paths) runs from the migrated inputs.
4. No other repository contains a hard-coded data path, an rclone call, or a download.
5. Every file in `Inputs/` has a `source.json`, and `catalog.json` lists every species with its source and `used_by`.
6. Every boundary file written by any repository passes `cichlid-data validate`.
7. Every family folder written by any method (TreeFam, Compara, OrthoFinder) passes `cichlid-data validate Results/<run>/families`, and `cichlid-data check --method OrthoFinder2026` passes, including its species tree.

## 8. Versioning

Semantic versions (`MAJOR.MINOR.PATCH`), tagged releases. Other repositories pin a range (`>=1.2,<2`).

- **MAJOR:** a layout change or a breaking schema change (a column removed or renamed).
- **MINOR:** new layout entries, units, downloaders, or optional columns.
- **PATCH:** fixes.

Every results folder's `run.json` records the coordinator version and catalog checksum used.

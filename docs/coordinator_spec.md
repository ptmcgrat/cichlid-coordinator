# Coordinator repository: instruction sheet

*For the Central Project Coordinator workspace and repository (`cichlid-coordinator`). October 2026.*

## 1. Role

The coordinator is the one package every other repository depends on. It owns:

1. **File locations:** one layout, identical relative paths on every machine and in the cloud.
2. **Transfers:** all movement between local machines and the cloud (Dropbox, via rclone).
3. **Inputs:** every download, the shared proteomes (one protein per gene), the species catalog, and provenance records.
4. **File formats** at the boundaries between repositories: readers and writers that check them.

It does **not** run tree methods, make orthology calls, analyze results or serve the browser. Those live in the method and analysis repositories, which use the coordinator through its public API only.

**Rule for every other repository:** never build a data path by hand, never call rclone directly, never download an input yourself. Ask the coordinator.

## 2. Repository layout

```
cichlid-coordinator/
├── pyproject.toml           package cichlid_data; console script cichlid-data
├── config/
│   ├── species.tsv          project species (edited by hand)
│   ├── species_sets.toml    which species each method uses
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
│   ├── downloaders/         ncbi.py, local.py (done); families.py, annotations.py, expression.py, populations.py
│   ├── migrate.py           one-time move from the old layout
│   └── cli.py               the cichlid-data command
├── tests/                   all tests run offline against LocalBackend
└── docs/                    project manifest, species lists, this sheet
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
│   │   ├── families/<family>/         alignment, tree.nwk, placements, done.json
│   │   ├── names/
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

`ResultsLayout` methods: `assign(prefix)`, `family(family)`, `orthology(prefix)`, `names()`, `residues(prefix)`, `run_record()`, and `unit(part)` (the `DataUnit` for a part).

### DataUnit (`units.py`)

A named, transferable piece of the tree, with a policy:

| Policy | Meaning | Examples |
|---|---|---|
| `immutable` | Existing files never change. File-by-file units may gain new files (e.g. a `source.json`); tar-packed units can't change at all | downloads, TreeFam 9, local species sources |
| `regenerable` | Pipeline output; may be replaced by a newer run | proteomes, Results/*, Analyses/* |
| `on_demand` | Large; never pulled as part of a bulk pull, only by `ensure` | population VCFs, Daniocell raw data |

`units.all()`, `units.get(name)`, `units.for_path(path)`.

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

Commands: `cichlid-data catalog build | show | missing | validate`, with `--method` and `--role`.

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
| Assignment `<prefix>.assign.tsv` | `gene_id, protein_length, best_family, score, evalue, protein_coverage, family_coverage, protein_region, n_families, second_family, second_score, fusion_partner, fusion_partner_score, fusion_partner_region, fusion_partner_family_coverage, flags` |
| Family folder `families/<family>/` | `tree.nwk` (final tree, leaf names = gene IDs or reference protein IDs), the alignment, `done.json` (inputs hash) |
| Orthology calls `<prefix>.orthologs.tsv` | Gene columns `gene_id, ncbi_gene_id, ncbi_symbol, ncbi_description, seqid, start, end, strand, mitochondrial, protein, protein_length, treefam_family, tree_mode, placement_status, placement_lwr`; then, for each reference (`human`, `mouse`, `zebrafish`, in that order), the block `<ref>_relationship, <ref>_ortholog_type, <ref>_symbol, <ref>_official_id, <ref>_ensembl_gene, <ref>_treefam_protein, <ref>_description, <ref>_name_status, <ref>_ref_in_subtree, <ref>_own_species_in_subtree, <ref>_median_cichlid_copies, <ref>_distance, <ref>_call_support`; then `human_uniprot, ncbi_symbol_agrees, confidence, caveats` |
| Run record `run.json` | `method, method_version, parameters, started, finished, coordinator_version, catalog_md5, species` |

Allowed values: `relationship` in {`Ortholog`, `Paralog`, `NoHomolog`}; `ortholog_type` in {`OneToOne`, `OneToTwo`, `OneToMany`, empty}; `confidence` in {`high`, `medium`, `low`}.

Booleans are written `True`/`False`; empty strings mean missing. Duplicate gene IDs within a species get a suffix (`Mzeb_101467881.2`), which the ID rules allow. The proteome check also cross-checks the genes table beside it: every gene with a chosen protein has exactly one protein, and vice versa. The orthology schema accepts any set of references, in blocks, so a method calling orthologs to other species stays valid. Stage 7's combined calls are a comparison output and are not yet specified (the comparison is being redesigned).

The schemas were checked against files written by the pipeline's own Stage 2, 3 and 6 code: all valid.

`cichlid-data validate PATH...` checks files (the kind is inferred from the name, or given with `--kind`).

A schema change is a version change of the package (Section 8).

### Proteomes (`proteomes.py`)

From each species' download, one protein per protein-coding gene: the pipeline's Stage 2, with its parsing and selection code moved over unchanged (the longest complete isoform, else the longest; NCBI's low-quality models last; mitochondrial genes flagged; residues outside the standard alphabet replaced and counted). Outputs are written through the schemas, and a species whose proteome fails validation leaves nothing behind. Species whose proteome is newer than their download are skipped unless `--force`; a missing download is fetched from the cloud first. Afterwards the catalog is rebuilt, and `--push` uploads the proteomes. QC rows for each species are kept in `stage2_summary.tsv`, one row per species.

`cichlid-data proteomes [--species P,...] [--select complete-longest|longest] [--force] [--push]`

### Downloaders (`downloaders/`)

Each subclass of `Downloader` implements `fetch()` and `verify()`; the base class writes the `SourceRecord`, updates the catalog where relevant, and pushes. Required sources:

| Module | Sources |
|---|---|
| `ncbi.py` | Done (v0.6.0): the pipeline's Stage 1, moved over unchanged. Genomes and annotations for every `species.tsv` row: from NCBI Datasets (`source=ncbi`; verified against NCBI's checksums, resumable, checked against NCBI's metadata first) or copied from local files (`source=local`; paths may be relative to the data root). A species missing locally but in the cloud is pulled rather than fetched again. `cichlid-data download genomes [--species ...] [--check-only] [--force] [--push]`, then the catalog is rebuilt |
| `local.py` | Done (v0.6.0): `prepare_local_species.py`, moved over unchanged. Converts a non-NCBI annotation into the download files (`cichlid-data prepare-local --prefix P --cds ... [--proteins ...]`; output defaults to `Inputs/Proteins/local_sources/<P>/prepared/`) |
| `families.py` | TreeFam 9 (HMMs, trees); Ensembl Compara dumps and pruning to the panel |
| `annotations.py` | GO (go-basic.obo, GAFs), PANTHER classification, Ensembl names, 2018 legacy files |
| `expression.py` | HPA, Bgee, Zebrahub, Daniocell, Brawand (nf-core outputs), with their summaries |
| `populations.py` | Blumer et al. 2025 VCFs and supplementary tables |

## 5. Command line

```
cichlid-data where <thing> [args]          print a path, e.g. "where proteome Mzeb"
cichlid-data env                            shell exports of the main paths
cichlid-data status [unit]
cichlid-data pull <unit> | push <unit>      --dry-run on every transferring command
cichlid-data download <source> [--species ...]
cichlid-data proteomes [--species ...]
cichlid-data catalog build | validate
cichlid-data validate <file>                check any boundary file against its schema
cichlid-data check [--method M]             every input a method needs is present and verified
cichlid-data migrate --dry-run | --apply
```

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

## 8. Versioning

Semantic versions (`MAJOR.MINOR.PATCH`), tagged releases. Other repositories pin a range (`>=1.2,<2`).

- **MAJOR:** a layout change or a breaking schema change (a column removed or renamed).
- **MINOR:** new layout entries, units, downloaders, or optional columns.
- **PATCH:** fixes.

Every results folder's `run.json` records the coordinator version and catalog checksum used.

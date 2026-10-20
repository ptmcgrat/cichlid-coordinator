# Coordinator repository: instruction sheet

*For the Central Project Coordinator workspace and repository (`cichlid-coordinator`). October 2026. Current at v0.10.0.*

## 1. Role

The coordinator is the one package every other repository depends on. It owns:

1. **File locations:** one layout, identical relative paths on every machine and in the cloud.
2. **Transfers:** all movement between local machines and the cloud (Dropbox, via rclone).
3. **Inputs:** every download, the shared proteomes (one protein per gene), the species catalog, and provenance records.
4. **File formats** at the boundaries between repositories: readers and writers that check them.

It does **not** run tree methods, make orthology calls, analyze results or serve the browser. Those live in the method and analysis repositories, which use the coordinator through its public API only.

**Rule for every other repository:** never build a data path by hand, never call rclone directly, never download an input yourself. Ask the coordinator.

**The repositories that use it.** Five projects, each with its own repository. The boundary between them is the coordinator's (layout, schemas, catalog, transfers), so none can drift from it. Names marked * are proposed, to be confirmed.

| # | Repository | Reads (through the coordinator) | Writes | Instruction sheet |
|---|---|---|---|---|
| 1 | `cichlid-trees`: TreeFam and Compara | proteomes, the catalog, TreeFam 9 and Compara 116 | per run: family folders (tree, `combined.aln`, `done.json`), assignment tables for each placed species, `run.json` | `trees_repo_spec.md` |
| 2 | `cichlid-orthofinder` | proteomes, the catalog, the OrthoFinder species tree | per run: family folders (orthogroups), assignment tables, `native/`, `run.json` | `orthofinder_repo_spec.md` |
| 3 | `cichlid-homology`*: the homology decider | the families of 1 and 2 (all three methods), a species-to-clade table, reference gene names | orthology calls per method and combined (human, mouse, zebrafish), retained ohnologs, cichlid-specific duplications and paralogs | to be written |
| 4 | `cichlid-analysis`*: additional data | population variation (VCFs, inversions, genotypes), alignments, the calls of 3 | residue maps and amino-acid changes, inversion status, genotypes | to be written |
| 5 | `cichlid-server`*: the web server | all of the above, expression, annotations | the public site | to be written |

The flow is 1 and 2 (independent of each other), then 3, then 4, then 5: the server is built from the analysis outputs. Two rules follow from it:

- **Only project 3 makes orthology calls**, and it makes them for all three methods with one rule, so differences between methods reflect the trees, not the calling logic. Projects 1 and 2 stop at families and assignments.
- **The residue step belongs to project 4.** It maps protein positions onto alignment columns. Whether `combined.aln` alone is enough depends on the part: the placement methods trim query-only insertion columns (MAFFT `--keeplength`), so `combined.aln` cannot give them an exact map. A part like that supplies `residue_map.tsv` for the delivered species of the release (an optional by-product, its format to become a schema); project 4 derives the map itself where `combined.aln` holds every residue.

**Releases (v0.9.0).** One release number is shared by all projects (`v1`, `v2`, ...). A release is a named state of the whole resource: the species, the input versions (TreeFam 9, Compara 116, the OrthoFinder version), the species tree and the delivered species. A new species, a better annotation or a new Compara release makes a new release. A release is defined in `config/releases.toml`, and its status is one of `planned`, `building`, `frozen`, `published` and `superseded`; the owner of the project changes it by editing that file. Results live under `Results/<release>/<Part>/` (parts `TreeFam`, `Compara`, `OrthoFinder`, `Homology`, `Analysis`). `cichlid-data release freeze v1` archives the proteomes, catalog and species trees of a release in `Releases/v1/`, which never changes afterwards, so the release can be reproduced after a species' annotation improves. Release 1 delivers calls for *M. zebra* and *A. calliptera* (`deliver`), while TreeFam and Compara place the 12 project species. The *A. calliptera* gene models are the newest annotation, and the variant data are on fAstCal1.2 (GCF_900246225.1), so project 4 lifts protein coordinates over; both accessions are recorded (the second in the release's `inputs`).

## 2. Repository layout

```
cichlid-coordinator/
├── pyproject.toml           package cichlid_data; console script cichlid-data
├── config/
│   ├── species.tsv          project species (edited by hand)
│   ├── species_sets.toml    every organism the project uses, and its role
│   ├── releases.toml        the releases: status, delivered species, input versions, and each part's species and species tree
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
└── docs/                    project manifest, species lists, this sheet, the sheets for the TreeFam/Compara and OrthoFinder repositories
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
│   ├── <release>/<Part>/              v1/TreeFam, v1/Compara, v1/OrthoFinder, v1/Homology, v1/Analysis:
│   │   │                              identical structure
│   │   ├── assign/<prefix>.assign.tsv
│   │   ├── families/<family>/         combined.aln, tree.nwk, done.json (the family contract, section 4;
│   │   │                              placement methods add queries.fa, ref.tree, placements.tsv, ...)
│   │   ├── names/
│   │   ├── native/                    the part's own outputs (OrthoFinder's result folders), as it wrote them
│   │   ├── orthology/<prefix>.orthologs.tsv
│   │   ├── residues/
│   │   └── run.json
├── Releases/<release>/                the frozen archive of a release (`release freeze`): release.json, catalog.json,
│                                      proteomes/, genes/, config/ (species trees, reference panel)
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
| `results(release, part)` | a `ResultsLayout` for `Results/<release>/<part>/` |
| `release(release)` | a `ReleaseLayout` for `Releases/<release>/` |
| `runs()`, `releases()` | the (release, part) pairs and the releases present on this machine |
| `analysis(name)` | `Analyses/<name>/` (the earlier analyses; new ones go in `results(release, "Analysis")`) |

A release is `v` and a number; a part is one of `TreeFam`, `Compara`, `OrthoFinder`, `Homology`, `Analysis` (anything else is refused). `ResultsLayout` methods: `assign(prefix)`, `family(family)`, `orthology(prefix)`, `names()`, `native()`, `residues(prefix)`, `run_record()`, `section(name)`, and `unit(section)` (the `DataUnit` for a section: `assign`, `families`, `names`, `native`, `orthology`, `residues`, `site`). `ReleaseLayout` methods: `manifest()`, `catalog()`, `proteome(prefix)`, `genes_table(prefix)`, `species_tree(name)`, `root()`.

### DataUnit (`units.py`)

A named, transferable piece of the tree, with a policy:

| Policy | Meaning | Examples |
|---|---|---|
| `immutable` | Existing files never change. File-by-file units may gain new files (e.g. a `source.json`); tar-packed units can't change at all | downloads, TreeFam 9, local species sources |
| `regenerable` | Pipeline output; may be replaced by a newer run | proteomes, Results/*, Analyses/* |
| `on_demand` | Large; never pulled as part of a bulk pull, only by `ensure` | population VCFs, Daniocell raw data |

`units.all()`, `units.get(name)`, `units.for_path(path)`.

Every part of a release has the same units: `results/<release>/<part>/assign`, `families`, `names`, `native`, `orthology`, `residues`, `site` and `run`. Each frozen release is one immutable unit, `releases/<release>`. `families` is sharded by family ID (the last three characters dropped, about 1,000 families per shard: `TF101001` in `TF101`, `C116_012345` in `C116_012`, `OG0001234` in `OG0001`). `native` holds a part's own outputs, **regenerable and tar-packed, one tar per top-level folder** (loose files go together into one extra tar, `_top`, repacked whenever any of them changes, so put files in subfolders) (so a second run kept beside the first uploads only its own tar); it is backed up and pulled like everything else. Neither is ever packed by hand.

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
- `config/species_sets.toml`: every organism the project uses (prefix, name, role, lineage, source; roles `project`, `relative`, `reference`, `outgroup`, `candidate`).
- `config/releases.toml`: each release (`status`, `description`, `deliver`, `inputs`) and its parts, each with its species: `queries = "project"` for TreeFam and Compara, an explicit list for OrthoFinder, none for Homology and Analysis. A part may name a `species_tree` or `reference_panel` (a file in `config/`).
- `config/compara_panel.tsv`: Compara's reference panel.

`Catalog.build(layout)` combines them with what is on disk into `catalog.json`; `Catalog.load(layout)` reads it.

`Species` fields: `prefix`, `species`, `common`, `role`, `lineage`, `source`, `accession`, `version`, `retrieved`, `in_manifest`, `status` (`proteome ready`, `downloaded`, `local source only`, `not downloaded`), `proteome`, `genes_table`, `proteome_md5`, `counts` (proteins, genes), `used_by` (method names), `notes`.

`catalog.md5()` is the MD5 of the `catalog.json` file the catalog was loaded from or written to (a catalog that was only built has no file yet and raises `FileNotFoundError`). **It identifies that one file, not what is in it:** the file records when it was built, so the value changes on every rebuild, even if nothing else did, and it can never match a frozen release's archive. To tie a run to the proteomes it used, use `catalog.proteomes_md5(prefixes)`: the MD5 of one line `<prefix> <proteome md5>` per species, sorted by prefix. It stays the same however often the catalog is rebuilt, and a frozen release can recompute it from its manifest (`proteomes_digest` in `catalog.py` is the recipe). It raises `KeyError` for a species without a proteome.

Queries: `catalog.species(release=..., part=..., role=..., status=...)` (a release alone gives the union of its parts), `catalog.missing(release, part)` (organisms whose proteome is not ready), `catalog.release_name(release=None)` (a release, checked; or the only one that is not superseded), `catalog.part(release, part)`, `catalog.releases[release]`, `catalog[prefix]`. `validate(release=None)` returns errors (bad prefixes, unknown roles; for a release: an unknown status or part, no species, members without an entry or still candidates, a delivered species that is not a project species or not in every tree part, an invalid species tree) and warnings (project species not in `species.tsv`). A `catalog.json` from v0.8 (which has methods, not releases) is refused with a message to run `catalog build`.

**Species tree.** A part may name a species tree in `releases.toml` (`species_tree = "species_tree_orthofinder.nwk"`, a file in `config/`). It must parse, be rooted and fully resolved (two children at the root and at every node), and have exactly the part's species as leaves. `catalog validate` and `catalog build` check it (an invalid tree stops the catalog from being written, so the tree and the species list cannot drift apart), and so does `cichlid-data check` for that part. `catalog.species_tree(release, part)` returns its path (`KeyError` if the part names none); its checksum is recorded with the other configuration files in `catalog.json`. The OrthoFinder tree is the one in `species_final.md`; a test keeps the two identical. Adding a species to OrthoFinder therefore means changing the list **and** the tree.

Commands: `cichlid-data catalog build | show | missing | validate | tree`, with `--release`, `--part` and `--role` (`catalog tree --part P` prints the species tree's path, for scripts; `--release` may be left out while only one release is live).

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
| Family folder `families/<family>/` | `tree.nwk` (final tree, leaf names = gene IDs or reference protein IDs), `combined.aln`, `done.json` (inputs hash, rooting): see "The family contract" below |
| Orthology calls `<prefix>.orthologs.tsv` | Gene columns `gene_id, ncbi_gene_id, ncbi_symbol, ncbi_description, seqid, start, end, strand, mitochondrial, protein, protein_length, treefam_family, tree_mode, placement_status, placement_lwr`; then, for each reference (`human`, `mouse`, `zebrafish`, in that order), the block `<ref>_relationship, <ref>_ortholog_type, <ref>_symbol, <ref>_official_id, <ref>_ensembl_gene, <ref>_treefam_protein, <ref>_description, <ref>_name_status, <ref>_ref_in_subtree, <ref>_own_species_in_subtree, <ref>_median_cichlid_copies, <ref>_distance, <ref>_call_support`; then `human_uniprot, ncbi_symbol_agrees, confidence, caveats` |
| Run record `run.json` | `release, method, method_version, parameters, started, finished, coordinator_version, catalog_md5, species`, and optionally `proteomes_md5`. `coordinator_version` is `cichlid_data.__version__`; `catalog_md5` is `Catalog.md5()` (the MD5 of that `catalog.json` file); `proteomes_md5` is `Catalog.proteomes_md5(species)`, which is what to use to show that a run used exactly the frozen release's proteomes (v0.10; write it in every new run record) |
| `residue_map.tsv` (in a family folder) | Tab-separated, header `gene_id, protein_length, aln_length, residues_not_in_alignment, columns`; one row per protein, sorted by `gene_id`, no duplicates. `columns` is a comma-separated list with one entry per residue in sequence order: the 1-based column of `combined.aln`, or `-` for a residue the alignment left out (an insertion MAFFT `--keeplength` dropped). Columns increase along the sequence; `aln_length` is the same on every row and is the length of `combined.aln`; `residues_not_in_alignment` counts the `-` entries. Optional; checked when present (`validate_residue_map`). Readers need `csv.field_size_limit(sys.maxsize)` |

Allowed values: `relationship` in {`Ortholog`, `Paralog`, `NoHomolog`}; `ortholog_type` in {`OneToOne`, `OneToTwo`, `OneToMany`, empty}; `confidence` in {`high`, `medium`, `low`}.

Booleans are written `True`/`False`; empty strings mean missing. Duplicate gene IDs within a species get a suffix (`Mzeb_101467881.2`), which the ID rules allow. The proteome check also cross-checks the genes table beside it: every gene with a chosen protein has exactly one protein, and vice versa. The orthology schema accepts any set of references, in blocks, so a method calling orthologs to other species stays valid. Stage 7's combined calls are a comparison output and are not yet specified (the comparison is being redesigned).

The schemas were checked against files written by the pipeline's own Stage 2, 3 and 6 code: all valid.

`cichlid-data validate PATH...` checks files (the kind is inferred from the name, or given with `--kind`: `genes`, `assign`, `orthologs`, `proteome`, `run`, `family`, `families`). A folder named `families` is checked family by family; given a catalog it also checks that leaf names starting with a species prefix are well-formed gene IDs.

#### The family contract

Every part writes `Results/<release>/<Part>/` in the same way, whatever built the trees, so the homology decider can read all of them alike. For any part of a release:

- **`assign/<prefix>.assign.tsv`**: the existing schema. For OrthoFinder, `best_family` is the orthogroup (empty for a gene in none) and the score columns may be empty.
- **`families/<family>/`**, one folder per family:
  - **Family ID**: letters, digits, `.`, `_` and `-`, starting with a letter or digit (`TF101001`, `C116_012345`, `OG0001234`). The ID is the folder name. OrthoFinder's families are its primary orthogroups (`OG…`), which have gene trees; its hierarchical orthogroups (`N0.HOG…`) are kept in `native/`. The ID rule would accept them, and sharding handles them, should that change.
  - **`tree.nwk`**: the final tree. Leaf names are unique and are **our gene IDs** (`<prefix>_<id>`) or **reference-protein IDs**. OrthoFinder renames sequences internally (`0_123`); its trees must be mapped back before they are written, and a leaf that looks like an internal name is refused. With the catalog's prefixes, a leaf whose text before its first `_` is a catalog prefix must be a well-formed gene ID; any other name (a bare prefix, or a reference ID that merely starts with the same letters) is a reference ID and is left alone.
  - **`combined.aln`**: required. The aligned sequences (FASTA, rows of equal length), with a row for every leaf of the tree; rows without a leaf are allowed. The placement methods write it beside their own files (`queries.fa`, `ref.tree`, `epa_result.jplace`, `placements.tsv`, `residue_map.tsv`, ...), which are not part of the contract and are not checked.
  - **`done.json`**: a JSON object. It holds the inputs hash that lets a rerun skip finished work, and says whether the tree is already rooted: `"rooted": true` with a `"rooting"` text (for example `"rooting": "orthofinder species tree"`), or `"rooted": false`. A part that builds some trees by placing genes on a fixed reference tree and others from scratch should also say which, in `"tree_mode"`: `place` or `denovo` (validated since v0.10; OrthoFinder writes `denovo`). `rooted` must be a true or false value; when true, `rooting` is required and the root of `tree.nwk` must have two children. **Consumers use an existing root and midpoint-root only otherwise**; a family that says nothing, or says something malformed, is treated as unrooted. `schemas.family_rooting(folder)` returns `(rooted, how)` so no consumer reimplements this rule.
- **A family may be missing.** A part skips a family it cannot build (the placement parts skip families with fewer than 4 sequences; OrthoFinder writes no gene tree for some orthogroups). Such a family has **no folder**, while the genes assigned to it keep their `best_family` in `assign/`. A consumer must not read that as an error; a part says how many families it skipped in its own summary (`native/summary/`).
- **A part's scratch is not part of the contract and is not pushed.** The `families` unit leaves out everything in a family's `stage5/` folder (the placement run's working files), RAxML checkpoint and start-tree files, and the residue tables; the `stage5.json` record beside the folder is kept.
- **`run.json`**: the run record, as above. Its `release` (and a family's `done.json` `release`, if present) must be the release of the folder it sits in.

`schemas.validate_family(folder, prefixes=None)` and `schemas.validate_families(folder, prefixes=None)` check all of this; `cichlid-data validate Results/<release>/<Part>/families` is the command.

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
| `export_seurat.R` | run with `Rscript` | The step after `download expression --source daniocell`: writes a Seurat object's matrix and cell metadata as plain files (`counts_x.bin`, `counts_i.bin`, `counts_p.bin`, `dims.txt`, `genes.txt`, `cells.txt`, `meta.tsv`, `source.txt`) for `summarize h5ad`. Moved verbatim from the pipeline (v0.9.1); the downloader prints the command with its path. Needs R and Seurat (see below) |
| `h5ad.py` | `summarize h5ad` | Single-cell atlases (Zebrahub, Daniocell) to per-cell-type summaries; needs the `[expression]` extra (NumPy, SciPy) |
| `rnaseq.py` | `summarize rnaseq` | Tissue RNA-seq (nf-core outputs, e.g. Brawand) to per-tissue summaries |
| `blumer.py` | `download populations` | Blumer et al. 2025 VCFs and tables (about 200 GB; resumable, checksum-verified) |
| `blumer_samples.py` | `populations samples` | Sample metadata for the VCFs (Table S1, optionally the lab sheet; `[samples]` extra for gspread) |
| `compara.py` | `prepare compara` | Prunes the Compara dumps to `config/compara_panel.tsv`; the dumps, panel and output folder default to their layout paths The output `prepared/` holds `trees/`, `members.tsv`, `protein_to_species.txt`, `search_db.fa` (unaligned sequences, tool-independent), `family_audit.tsv` and `prepare.json`. It holds **no DIAMOND database**: a database depends on the DIAMOND version, so the part that searches builds its own (and never writes into an input folder). The audit statuses are `ok`, `tiny_tree`, `cichlid_only`, `no_panel_members`, `mismatch` and `tree_parse_error`; a consumer handles all of them. |

`newick.py` (the pipeline's Newick parser) is a shared module of the coordinator, so the pipeline can import it from here.

**System requirements outside Python.** `export_seurat.R` needs R and the Seurat package, which are not pip packages. On mzebra they are in a small conda environment, `seurat` (R 4.5.3, Seurat 5.5.1); the rest of the coordinator runs in the `cichlid-trees` environment. The export is a one-off per Daniocell release: `conda activate seurat`, run the printed `Rscript` command, then go back to the other environment for `summarize h5ad`.

### Check (`check.py`)

`cichlid-data check [--release R] [--part P] [--fetch] [--deep]` confirms, for each part of a release (default: every release that is not superseded), that everything it needs is present and sound: the proteomes of its species (valid against the schemas, with the genes-table cross-check), its resources (TreeFam 9 or the prepared Compara trees for those parts; the legacy name tables for Homology), each verified against its source record (sizes; `--deep` also checksums), and, for a part that names one, its species tree. For a frozen or published release it also checks its archive (present and intact) and that the current proteomes and configuration have not drifted from it. `--fetch` pulls anything missing from the cloud first. It exits with an error if anything is missing or invalid, so run scripts can call it before starting.

### Releases (`releases.py`)

`freeze(layout, catalog, release)` writes the archive `Releases/<release>/` from a freshly built catalog. It refuses unless the release is valid, its status in `releases.toml` is `frozen` or `published` (freezing is the owner's decision), every species it uses has a proteome, and it has not been frozen already (a frozen release never changes). The archive holds `release.json` (the release as defined, its status, the coordinator version, the species with their proteome checksums, and the size and MD5 of every file), `catalog.json`, `proteomes/<prefix>.faa`, `genes/<prefix>.genes.tsv` and `config/` (the species trees and reference panel). It is written to a temporary folder and renamed, so it appears whole or not at all. Its unit is immutable: after `cichlid-data push releases/<release>` an edit is refused. Runs of a frozen release can read the proteomes from the archive (`layout.release(r).proteome(p)`).

`run_records(layout, release)` compares each `Results/<release>/<part>/run.json` with the archive: a run whose `proteomes_md5` differs from what the release froze (for the run's species) is a problem, a run that records a species the release does not have is a problem, and a run without `proteomes_md5` is only noted. `release validate` prints these. The manifest records `proteomes_md5` for all the release's species.

`verify_archive(layout, release, deep=True)` returns what differs in the archive from when it was frozen (missing, resized or changed files, and files that were added). `drift(layout, release)` returns how the current proteomes and configuration differ from a frozen release: not an error in the data, but a sign that a newer release is due. `load_manifest(layout, release)` reads `release.json`.

## 5. Command line

```
cichlid-data init [--root DIR] [--remote REMOTE] [--force]   once per machine
cichlid-data where <thing> [args] [--cloud]          print a path, e.g. "where proteome Mzeb", "where results v1 OrthoFinder native", "where release v1"
cichlid-data env                                     shell exports of the main paths
cichlid-data units | status [pattern]
cichlid-data pull | push <pattern>... [--force] [--dry-run]
cichlid-data ensure <relative path>                  fetch one file or shard if missing
cichlid-data retire <unit key> [--dry-run]           delete a unit's or shard's cloud copy
cichlid-data migrate [--apply] [--cleanup]           the one-time move to the new layout
cichlid-data provenance write [--overwrite] [--md5-limit-gb N] | check [--deep] | show <folder>
cichlid-data catalog build | show | missing | validate | tree [--release R] [--part P] [--role X]
cichlid-data release show | validate | freeze [R] [--fast]       releases and their frozen archives
cichlid-data validate <file or folder>... [--kind K] [--prefix P]
cichlid-data download genomes [--species P,...] [--check-only] [--force] [--keep-zip] [--no-fetch] [--push]
cichlid-data prepare-local --prefix P ...
cichlid-data proteomes [--species P,...] [--select complete-longest|longest] [--force] [--no-fetch] [--push]
cichlid-data download annotations | expression | populations ...
cichlid-data populations samples ...
cichlid-data prepare compara ...
cichlid-data summarize h5ad | rnaseq ...
cichlid-data check [--release R] [--part P] [--fetch] [--deep]
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
7. Every family folder written by any part (TreeFam, Compara, OrthoFinder) passes `cichlid-data validate Results/<release>/<Part>/families`, and `cichlid-data check --release v1 --part OrthoFinder` passes, including its species tree.
8. A frozen release's archive is intact (`cichlid-data release validate`), and `cichlid-data check` shows no drift from it.

## 8. Versioning

Semantic versions (`MAJOR.MINOR.PATCH`), tagged releases. From 1.0.0, other repositories pin the major version (`>=1.2,<2`). Before 1.0.0 a MINOR version may change the layout API, so they pin one minor (`>=0.8,<0.9`).

- **MAJOR:** a layout change or a breaking schema change (a column removed or renamed).
- **MINOR:** new layout entries, units, downloaders, or optional columns.
- **PATCH:** fixes.

Every results folder's `run.json` records the coordinator version and catalog checksum used.

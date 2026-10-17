# Coordinator handoff

*This file is the coordinator project's record of its own state: what it owns, where things are, how routine jobs are done, what was decided and why, what is open. A new session, or a new person, starts here. The other projects do not edit it; their contract with the coordinator is `coordinator_spec.md`. Handoff of the Central Project Coordinator role, October 2026; updated at coordinator v0.8.0 (the family contract). Read with `project_manifest.md` (why the project exists) and `coordinator_spec.md` (what the coordinator must provide).*

## 0. Where things stand (v0.8.0)

- **All 35 OrthoFinder proteomes are ready.** Betta (*Bspl*) and stickleback (*Gacu*) are in `species.tsv` and `species_sets.toml`, downloaded, their proteomes built and pushed; `catalog missing --method OrthoFinder2026` reports none missing, and `check` reports TreeFam2026 16 of 16, Compara2026 18 of 18, OrthoFinder2026 ready.
- **The organization is decided: five projects, each its own repository** (updated from two): 1 `cichlid-trees` (TreeFam and Compara: assignment, alignment, placement), 2 `cichlid-orthofinder` (builds families with OrthoFinder; `native/` keeps its own outputs), 3 the homology decider (all orthology calls, for all three methods, with one rule; ohnologs, cichlid-specific duplications, paralogs), 4 the analysis project (residue maps, amino-acid changes, inversion status, genotypes), 5 the server (built from the analysis outputs). The boundary is the coordinator's: layout, schemas, catalog, transfers. Names of 3 to 5 (`cichlid-homology`, `cichlid-analysis`, `cichlid-server`) are proposed, to be confirmed. See `coordinator_spec.md` section 1.
- **v0.8.0 is the coordinator's side of that boundary:** the family contract (including whether a tree is already rooted), the OrthoFinder species tree as configuration, and the `native/` unit. The instruction sheets are `orthofinder_repo_spec.md` and `trees_repo_spec.md`; the sheets for projects 3 to 5 are to be written.

## 1. What the coordinator owns

- **The data tree:** one layout (`Inputs/`, `Results/`, `Analyses/`, `Archive/`), identical on every machine and in the cloud.
- **All transfers** between machines and Dropbox. No other repository calls rclone.
- **All inputs:** genome and annotation downloads, the 35 proteomes, TreeFam 9, the prepared Compara 116 trees, GO, PANTHER, the legacy name tables, expression and population data.
- **The species decisions:** `config/species.tsv` (what is downloaded), `config/species_sets.toml` (roles, and which species each method uses), `config/compara_panel.tsv` (Compara's reference panel), `config/species_tree_orthofinder.nwk` (the species tree OrthoFinder must use).
- **Provenance:** a source record for every input (`config/sources.toml`, `source.json`, Stage 1's `download.json`).
- **The file formats** that cross repositories (`cichlid_data.schemas`), including the family contract every method's results follow, and the shared Newick parser.
- **The project documents:** the manifest, the species lists and figure, this handoff, the specification.

Other repositories use the coordinator only through its public API (`cichlid_data`) and command line, and pin a version range.

## 2. Where everything is

| Thing | Location |
|---|---|
| Repository | `github.com/ptmcgrat/cichlid-coordinator` (branch `main`, tags `v0.x.y`) |
| Laptop clone | `~/GaTech Dropbox/CoS/BioSci/BioSci-McGrath/Applications/cichlid-coordinator` |
| Compute server | `bio-mcgrath-dropbox@mzebra.ebb.gatech.edu`, clone `~/cichlid-coordinator`, conda env `cichlid-trees` |
| Data root (mzebra) | `~/Temp/TreeFamData` |
| Cloud copy | `ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/TreeFamData` (cloud index in `.cichlid-data/index.json`) |
| Machine settings | `~/.config/cichlid-data/config.toml` (written by `cichlid-data init`) |
| TreeFam and Compara repository | `github.com/ptmcgrat/cichlid-trees` (not yet using the coordinator; see section 7; instruction sheet `trees_repo_spec.md`) |
| OrthoFinder repository | `cichlid-orthofinder` (the student's; instruction sheet `orthofinder_repo_spec.md`) |
| Other projects | the homology decider, analysis and server (repository names proposed, section 0); no repository yet |

Old result folders are still in the data root, untouched: `Rebuild2026/`, `Compara2026/`, `Combined2026/` (phylogenies, calls, browser site, analyses of the previous runs), and `Data/Blumer2025/coding` and `inversions` (variant analyses). They are kept for comparison with the new runs.

## 3. How to work on the coordinator

**Code changes.** Every change comes with tests; the suite runs offline (fake cloud, fake rclone, fake web):

```bash
PYTHONPATH=src:tests python -m unittest discover tests        # 154 tests at v0.8.0 (103 at v0.7.1)
```

**Moving code in from elsewhere** (as with Stages 1 and 2): move the logic verbatim, port its tests with their fixtures and expectations unchanged, and change only paths, transfers and formats. Then verify on real data that outputs are identical (as was done for the proteomes: byte-for-byte), and explain any difference before accepting it.

**Versions.** Semantic: MAJOR for layout or breaking schema changes, MINOR for new features, PATCH for fixes. Bump `pyproject.toml` and `src/cichlid_data/__init__.py` together, commit, tag `vX.Y.Z`, push with `--tags`. Tag 1.0.0 once the pipeline runs on the coordinator.

**Installing on a machine.** Always through the environment's own Python (on mzebra, a stray Python 3.6 pip in `~/.local` otherwise wins):

```bash
conda activate cichlid-trees
python -m pip install -e ~/cichlid-coordinator          # once; editable, so git pull updates it
python -m pip install -e "~/cichlid-coordinator[expression]"   # adds NumPy and SciPy for summarize h5ad
```

**Delivering a change from a Claude conversation.** Claude produces a zip of the repository; on the laptop: unzip, `rsync -av --exclude config/species.tsv --exclude .git` over the clone, run the tests, commit, tag, push; then `git pull` on mzebra. `species.tsv` is excluded so hand edits are never overwritten; edit it directly in the repository and commit.

## 4. Procedures

**New machine.** Clone the repository, install as above, `cichlid-data init` (edit `~/.config/cichlid-data/config.toml` if the data root differs), then `cichlid-data pull 'proteins/*' 'families/*'` and whatever else the work needs.

**Add a species.**
1. Look up its RefSeq reference annotation (lookup script in the conversation record; it queries NCBI's Datasets API for `filters.reference_only` and `filters.has_annotation`).
2. Add a row to `config/species.tsv` and an `[[organism]]` entry (role, lineage) to `config/species_sets.toml`; add the prefix to a method's list if it should be used. Prefixes: letters and digits, unique ignoring case. **For OrthoFinder, also add the species to `config/species_tree_orthofinder.nwk` (and the tree in `species_final.md`)**: `catalog build` refuses to write while the tree and the species list disagree.
3. `cichlid-data download genomes --species P --check-only`, then without `--check-only` and with `--push`.
4. `cichlid-data proteomes --species P --push`; review its row in `Inputs/Proteins/proteomes/stage2_summary.tsv` (proteins vs NCBI's count, genes without a protein, partial models, median length 350–460 aa).
5. `cichlid-data catalog build`; `cichlid-data check`.
6. Update `species_final.md` and, if the lineage is new, the figure (`make_species_phylogeny.py`).

**A species not annotated by NCBI.** `cichlid-data prepare-local --prefix P --cds ... --genome ... [--functional ...] [--proteins ...]`; a `source=local` row with paths relative to the data root (as for `Esur`). Ask the authors for their protein FASTA: translating their coding sequences is unreliable (it was for *Etroplus*).

**Refresh an input** (a new NCBI annotation, a new GO release). Downloads are immutable: refreshing a species means `download genomes --species P --force`, then `push --force` for its shard, then rebuilding its proteome. Record the reason in the commit message. For other inputs, rerun the downloader, then `provenance write --overwrite` for that resource, then push.

**Before any method runs.** `cichlid-data check --method <Method> --fetch`. For OrthoFinder this includes its species tree; the path for `-s` comes from `cichlid-data catalog tree --method OrthoFinder2026`.

**Check a method's families.** `cichlid-data validate Results/<run>/families` checks every family folder against the family contract (the catalog supplies the species prefixes). Run it on the converter's output before anything calls it.

**Transfer conflicts.** `status` shows `local changed`, `cloud changed` or `both changed`. Pull before pushing; `--force` only when you are sure which copy is right. A stale lock (a push that died) can be overridden with `--force` after six hours. A shard that should no longer exist: `cichlid-data retire <key>`.

**Provenance.** After adding or refreshing any input: `cichlid-data provenance write` then `provenance check`. Fill in the "to confirm" licences and citations in `config/sources.toml` (a good student task), then `provenance write --overwrite`.

## 5. Decisions made, and why

| Decision | Reason |
|---|---|
| One coordinator package that every repository depends on | One implementation of paths, transfers and formats; no drift between repositories |
| Layout `Inputs/ Results/ Analyses/ Archive/`, identical relative paths locally and in the cloud | Any path can be resolved on any machine |
| Units with policies: immutable (existing files never change; file-by-file units may gain files, tar units can't change at all), regenerable, on-demand | Downloads must never be silently replaced; huge data is never bulk-pulled |
| Tars hold paths relative to the unit's parent; families sharded by ID prefix, downloads by species | Units can move without repacking; reruns upload only changed shards |
| Migration moved inputs only; old results left in place; the 2018 resource's files archived | TreeFam and Compara are rerun with the new species anyway; old results stay for comparison |
| Stage 1's `download.json` is the species' source record | It already holds accession, metadata, checksums and date; adding a file to an immutable shard is refused |
| Configuration in TOML, not YAML | Read by the standard library; the package needs only `requests` (plus optional extras) |
| Proteomes: genes on alternate loci and patches are left out, identified by sequence role | They are second copies (human GRCh38: 3,384) that would look like recent duplications; unit names differ between assemblies (mouse's primary unit is "C57BL/6J"), roles don't |
| Proteomes: proteins linked per gene; a protein shared by two loci kept once | Fixes a Stage 2 bug that left ~870 human genes without proteins; avoids identical false duplicates (human X/Y pseudoautosomal genes) |
| Line endings `\n` in tables written by the coordinator | Stage 2 wrote `\r\n` by default; content unchanged |
| Download reports (`download_summary.tsv`, `check_summary.tsv`) are their own regenerable units | They change every run; inside the immutable downloads unit they blocked pushes |
| OrthoFinder species set: 35 (see `species_final.md`); lamprey and *Ciona* dropped, eel added; bowfin and *Amphilophus* candidates only (no RefSeq annotation) | Serves the two aims: human/zebrafish relationships and teleost duplicate retention |
| *(Superseded by the five-project row below.)* Two method repositories: `cichlid-trees` (TreeFam, Compara, and a shared core that calls any method's families) and `cichlid-orthofinder` (runs OrthoFinder, converts its trees, never calls) | Every method's calls come from one rule, so differences between methods reflect the trees, not the calling logic; the student's repository has a small, checkable job |
| The family contract is the boundary: family folder = `tree.nwk` (leaves are our gene IDs or reference-protein IDs), `combined.aln` (required), `done.json`; `assign` for OrthoFinder has the orthogroup as `best_family` and empty scores | The core reads any method's results the same way; OrthoFinder's internal sequence names (`0_123`) can't leak through because validation refuses them |
| Rooting is recorded in `done.json` (`"rooted": true`, `"rooting": "..."`); absent or malformed means unrooted; the consumer (project 3) uses an existing root and midpoint-roots only otherwise | OrthoFinder roots its gene trees with the species tree, which is better than a midpoint; placement-method trees stay as they are today. Failing safe (treat as unrooted) never trusts a malformed claim |
| Family IDs: letters, digits, `.`, `_`, `-`, not starting with `_`; shards drop the last three characters | Covers `TF101001`, `C116_012345`, `OG0001234` (and `N0.HOG0001234`, should hierarchical orthogroups ever become families); a leading `_` is refused because `_top` is the shard for loose files. Fixed in v0.8.0: a folder named like `N0.HOG0001234` that was not yet on the machine was taken for a loose file (shard `_top`), so `ensure` would have looked in the wrong shard |
| The alignment of a family is `combined.aln`, and it is required | It is what Stage 4 of the placement methods already writes beside `queries.fa`, `ref.tree`, `epa_result.jplace`, `placements.tsv` and `residue_map.tsv`; one name for every method means the shared core reads any method's families the same way. (An earlier draft said `alignment.fa`; checked against the real folders and changed.) |
| OrthoFinder's families are its primary orthogroups (`OG…`); the hierarchical orthogroups (`N0.HOG…`) stay in `native/` | Orthogroups have gene trees and the hierarchical ones do not; the homology decider needs trees |
| `Catalog.md5()` supplies `catalog_md5` for every `run.json`; `coordinator_version` is `cichlid_data.__version__` | Repositories never compute either by hand, so records cannot disagree about what a checksum is |
| A leaf is checked against the catalog only if the text before its first `_` is a catalog prefix | `<prefix>_` is the gene-ID form; a reference ID that merely starts with the same letters, or is the bare prefix, is not ours and is not refused |
| Five projects (trees, OrthoFinder, homology decider, analysis, server), in that order of flow; only project 3 makes orthology calls | Replaces the earlier plan of a shared calling core inside `cichlid-trees`. One rule for all three methods; projects 1 and 2 stop at families and assignments |
| The residue step (protein positions on alignment columns) belongs to project 4 | It needs only `combined.aln` and the proteome, so one tool serves every method. The residue files in existing family folders are left as they are |
| One release number shared by all projects (`v1`, `v2`, ...); a new species, a better annotation or a new Compara release makes a new major release; statuses `planned`, `building`, `frozen`, `published`, `superseded`, changed by the owner | Project 3 combines the three methods at one state, so they need one version to describe it. Inputs change under a result otherwise. **Decided, not built (v0.9.0)** |
| Release layout `Results/<release>/<Part>/` (parts `TreeFam`, `Compara`, `OrthoFinder`, `Homology`, `Analysis`), replacing `Results/<Method><Year>/`; `release freeze` archives the proteomes, catalog and species tree | No project writes into the layout yet, so changing it is cheap now and costly later. Proteomes are regenerable, so without a freeze a better annotation would overwrite what a release was built from |
| Release 1 delivers calls for *M. zebra* and *A. calliptera* (`deliver`); TreeFam and Compara place the 12 project species | Delivered is separate from placed; more species can be requested later |
| The *A. calliptera* gene models are the newest annotation; the variant data are on fAstCal1.2 (GCF_900246225.1); project 4 lifts protein coordinates over | Quality of the newest annotation; both accessions are recorded in the release. The accession in the real `species.tsv` for Acal is to be checked |
| Below 1.0.0, other repositories pin one minor version (`>=0.8,<0.9`) | A minor version may change the layout API (v0.9.0 will) |
| The old results (`Rebuild2026`, `Compara2026`, `Combined2026`) are kept, not deleted, until the `cichlid-trees` update has been verified against them | They are the baseline for "outputs identical on real data" and were the real-data test set for the validators |
| The OrthoFinder species tree is configuration (`config/species_tree_orthofinder.nwk`), named by `species_tree` in `species_sets.toml`, validated against the catalog's species, with its checksum in `catalog.json` | One source for what OrthoFinder is given and what the paper reports; the tree cannot silently disagree with the species list |
| `Results/<run>/native/` for a method's own outputs: regenerable, tar-packed, **one tar per top-level folder** | Backed up and pullable like everything else. One tar per folder rather than one for the whole unit, so a second OrthoFinder run kept beside the first (`native/<name>/`) uploads only its own tar |
| The second OrthoFinder run (OrthoFinder's own species tree, for comparison) lives inside the main run's `native/`, not as a run of its own *(assumption, to confirm)* | It is a comparison, with no conversion to families and no calls; it needs no new run name, catalog method or check. Promote it to its own run if calls are ever wanted from it |
| `catalog build` and `check` validate the species tree | An edit to the species list without the tree stops at the earliest point |

## 6. Open items

- **"To confirm" fields** in `config/sources.toml` (TreeFam and PANTHER terms, Zebrahub and Daniocell licences, exact citations for Blumer, Katneni, Brawand).
- **Old cloud copies from the migration**, removed by hand with rclone, each after its replacement was checked (`migrate --cleanup` cannot reach them: it deletes an old copy only when the same run repacked it, and the trees were never repacked by a run that completed, while the genomes were repacked in a run without `--cleanup`; a later run plans nothing). Removed on 2026-10-17: `Rebuild2026/genomes.tar` (its 9 species, 45 files, are in the new per-species shards) and the loose folder `Data/OriginalTreeFiles/` (24,572 original TreeFam 9 files plus 17,678 hand-edited copies, `.fixed.aln`, `.fixed.nex` and one `.fixed.nhx`, made in 2018 so the Python library of the time could read the TreeFam files; all 8,836 originals and their alignments are read by the coordinator's parser). `Data/OriginalTreeFiles.tar` (replaced by `Inputs/Families/TreeFam9/trees.tar`, 24,573 entries equal to the local folder) was removed the same day.
- **Species caveats to remember:** *Etroplus* (authors' models; median protein 351 aa, shorter than all others); *A. centrarchus* (possible retained second haplotypes inflating its species-specific duplications); *H. burtoni* (old scaffold-level annotation, 9% partial models); elephant shark (scaffold-level). Bowfin and *Amphilophus* have no RefSeq annotation; *Pholidichthys* (sister family of cichlids) is unannotated: an NCBI annotation request would be worthwhile.
- **Small proteome remainders:** 6 duplicate-ID genes each in human and zebrafish (one GeneID at two loci with different proteins), kept as `….2`.
- **Cloud index repaired (October 2026):** the first migration run stopped at the `NewTreeFiles` error before repacking the TreeFam trees, and its index carry-over covered only the proteomes. TreeFam9 `hmm` and `protein_to_species` and Compara116 `dumps`, moved server-side, were verified present and identical in the cloud (`rclone check --size-only`, and a dry-run `copy` that transferred nothing), then pushed, which wrote only the index. `TreeFam9/trees` (3.55 GiB) and `Compara116/prepared` (3.7 GiB) were rebuilt and uploaded as tars (the old tar's paths, `OriginalTreeFiles/…`, differ from the unit's, `trees/…`). Every `families/*` unit is in sync. A future server-side move leaves the index behind in the same way: verify, then push (an `adopt` command is possible if a tar is ever too large to re-upload). `expression/raw` (33 files, 23.994 GiB) was also "local only"; it was verified identical in the cloud, then pushed (index only, nothing transferred). No unit is "local only" now.
- **TreeFam 9 input, as found:** all 8,836 trees parse and all 15,736 alignments read with equal row lengths. In a family's alignment the row names carry a `/start-end` suffix that the tree leaf names do not (checked on TF101001: equal once stripped). In the old *results* folders, every tree leaf has a row in the family's `combined.aln` (validated on all of Compara2026 and Rebuild2026). `Inputs/Families/Compara116/prepared.json` (in the cloud) is a leftover record from an older version of the Compara preparation; nothing reads it.
- **Size of `native/`:** OrthoFinder's output for 35 proteomes may be tens of GB. If so, the OrthoFinder project keeps the standard results and the files that matter (orthogroups, gene trees, alignments, species tree, logs) and leaves out the bulk working files, or `native` gets excludes. To decide once the first run exists.
- **v0.9.0 (decided, not built):** `config/releases.toml` (status, `deliver`, species per part, species tree, input versions); layout `Results/<release>/<Part>/` (`Layout.results(release, part)`); `Releases/<release>/` archive written by `release freeze`; `release show | validate | freeze`; `check --release`; a `release` field in `run.json` and `done.json`; the per-method species lists move from `species_sets.toml` into releases; the 2018 name tables (`Inputs/Annotations/Legacy2018`) leave the TreeFam and Compara requirements of `check` and become project 3's. `release export` (public CSVs with release and annotation accession) follows in v0.10.
- **Schemas to define with the other projects:** the orthology and combined-call tables (project 3; today's `orthologs.tsv` has TreeFam-specific columns), the outputs of project 4, and what project 5 reads. A species-to-clade table (every species of every method mapped to cichlid, Ovalentaria, percomorph, teleost, pre-genome-duplication and so on) is configuration the coordinator owns; `compara_panel.tsv` has a coarse `group` for 203 species, TreeFam 9's 109 species are not mapped, and lineages should come from NCBI Taxonomy.
- **Licences** for public release (`sources.toml`: TreeFam, PANTHER, Zebrahub, Daniocell) are a gate for project 5; owner: the project's PI.
- **Version 1.0.0** once the pipeline runs on the coordinator.

## 7. What comes next

Five projects plus the coordinator, each in its own project. The coordinator's part is the boundary between them: layout, schemas, catalog, transfers, releases.

1. **`cichlid-trees` update:** takes its paths, transfers, proteomes, inputs and formats from the coordinator; removes what now lives here (Stages 1 and 2, `newick.py`, `compara_prepare.py`, the download and summary scripts, its own sync code); hands calling to project 3 and residues to project 4; writes a `run.json` per run; reruns TreeFam and Compara with the 12 project species; verifies the nine baseline species against `Rebuild2026` and `Compara2026`. From `trees_repo_spec.md`.
2. **`cichlid-orthofinder`:** from `orthofinder_repo_spec.md`.
3. **Homology decider:** reads the families of 1 and 2, makes the calls for all three methods with one rule, combines them. Boundary sheet and schemas to be written with its owner.
4. **Analysis:** residue maps, amino-acid changes (lifting protein coordinates from fAstCal1.2 to the newest *A. calliptera* annotation), inversion status, genotypes. Boundary sheet to be written.
5. **Server:** built from the analysis outputs. Gated on licences. Boundary sheet to be written.
6. **Coordinator (this project):** tag v0.8.0; build v0.9.0 (releases); write the boundary sheets for 3 to 5 and the species-to-clade table; answer the other projects' requests (a missing API, a format that needs a field, a version to pin); keep `docs/` and `check` current; tag 1.0.0 once projects 1 and 2 run on the coordinator.

A request from another project is a request to the coordinator, not an edit to its code: the other projects do not change schemas, paths or transfers on their own side.

## 8. Documents for the coordinator workspace

| Document | Status |
|---|---|
| `HANDOFF.md` (this file) | Current |
| `project_manifest.md` | Current |
| `coordinator_spec.md` | Current (v0.8.0) |
| `species_final.md` | Current: the agreed species list (and the species tree; `config/species_tree_orthofinder.nwk` is the working copy) |
| `orthofinder_repo_spec.md` | Current (v0.8.0): the instruction sheet for the OrthoFinder repository |
| `trees_repo_spec.md` | Current (v0.8.0): the instruction sheet for the TreeFam and Compara repository |
| Sheets for the homology decider, analysis and server | To be written |
| `species_reference.md`, `species_phylogeny.svg/.png`, `make_species_phylogeny.py` | Current: the three tiers and the figure |
| `orthofinder_species.md` | Background: the draft that led to `species_final.md` (superseded where they differ) |
| `reorganization_proposal.md` | Historical: the proposal the layout was built from (implemented; see the spec) |

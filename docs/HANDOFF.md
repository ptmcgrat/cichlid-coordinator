# Coordinator handoff

*Handoff of the Central Project Coordinator role, October 2026; updated at coordinator v0.8.0 (the family contract). Read with `project_manifest.md` (why the project exists) and `coordinator_spec.md` (what the coordinator must provide).*

## 0. Where things stand (v0.8.0)

- **All 35 OrthoFinder proteomes are ready.** Betta (*Bspl*) and stickleback (*Gacu*) are in `species.tsv` and `species_sets.toml`, downloaded, their proteomes built and pushed; `catalog missing --method OrthoFinder2026` reports none missing, and `check` reports TreeFam2026 16 of 16, Compara2026 18 of 18, OrthoFinder2026 ready.
- **The repository organization is decided: two repositories.** `cichlid-trees` (TreeFam and Compara, plus the shared core: alignment and placement for the placement methods, and calling and residues for any method's results, so all calls come from one rule) and `cichlid-orthofinder` (runs OrthoFinder on the catalog's 35 proteomes, converts its gene trees into the standard family format, keeps OrthoFinder's own outputs in `native/`, never makes calls). The boundary is the coordinator's: layout, schemas, catalog.
- **v0.8.0 is the coordinator's side of that boundary:** the family contract (including whether a tree is already rooted), the OrthoFinder species tree as configuration, and the `native/` unit. The student's instruction sheet is `orthofinder_repo_spec.md`.

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
| Pipeline repository | `github.com/ptmcgrat/cichlid-trees`: TreeFam, Compara and the shared core (not yet using the coordinator; see section 7) |
| OrthoFinder repository | `cichlid-orthofinder` (the student's; instruction sheet `orthofinder_repo_spec.md`) |

Old result folders are still in the data root, untouched: `Rebuild2026/`, `Compara2026/`, `Combined2026/` (phylogenies, calls, browser site, analyses of the previous runs), and `Data/Blumer2025/coding` and `inversions` (variant analyses). They are kept for comparison with the new runs.

## 3. How to work on the coordinator

**Code changes.** Every change comes with tests; the suite runs offline (fake cloud, fake rclone, fake web):

```bash
PYTHONPATH=src:tests python -m unittest discover tests        # 147 tests at v0.8.0 (103 at v0.7.1)
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
| Two method repositories: `cichlid-trees` (TreeFam, Compara, and the shared core that calls any method's families) and `cichlid-orthofinder` (runs OrthoFinder, converts its trees, never calls) | Every method's calls come from one rule, so differences between methods reflect the trees, not the calling logic; the student's repository has a small, checkable job |
| The family contract is the boundary: family folder = `tree.nwk` (leaves are our gene IDs or reference-protein IDs), `alignment.fa`, `done.json`; `assign` for OrthoFinder has the orthogroup as `best_family` and empty scores | The core reads any method's results the same way; OrthoFinder's internal sequence names (`0_123`) can't leak through because validation refuses them |
| Rooting is recorded in `done.json` (`"rooted": true`, `"rooting": "..."`); absent or malformed means unrooted; the core uses an existing root and midpoint-roots only otherwise | OrthoFinder roots its gene trees with the species tree, which is better than a midpoint; placement-method trees stay as they are today. Failing safe (treat as unrooted) never trusts a malformed claim |
| Family IDs: letters, digits, `.`, `_`, `-`, not starting with `_`; shards drop the last three characters | Covers `TF101001`, `C116_012345`, `OG0001234` and `N0.HOG0001234`; a leading `_` is refused because `_top` is the shard for loose files. Fixed in v0.8.0: a folder named like `N0.HOG0001234` that was not yet on the machine was taken for a loose file (shard `_top`), so `ensure` would have looked in the wrong shard |
| The OrthoFinder species tree is configuration (`config/species_tree_orthofinder.nwk`), named by `species_tree` in `species_sets.toml`, validated against the catalog's species, with its checksum in `catalog.json` | One source for what OrthoFinder is given and what the paper reports; the tree cannot silently disagree with the species list |
| `Results/<run>/native/` for a method's own outputs: regenerable, tar-packed, **one tar per top-level folder** | Backed up and pullable like everything else. One tar per folder rather than one for the whole unit, so a second OrthoFinder run kept beside the first (`native/<name>/`) uploads only its own tar |
| The second OrthoFinder run (OrthoFinder's own species tree, for comparison) lives inside the main run's `native/`, not as a run of its own *(assumption, to confirm)* | It is a comparison, with no conversion to families and no calls; it needs no new run name, catalog method or check. Promote it to its own run if calls are ever wanted from it |
| `catalog build` and `check` validate the species tree | An edit to the species list without the tree stops at the earliest point |

## 6. Open items

- **"To confirm" fields** in `config/sources.toml` (TreeFam and PANTHER terms, Zebrahub and Daniocell licences, exact citations for Blumer, Katneni, Brawand).
- **Old cloud copies from the migration** (`Rebuild2026/genomes.tar`, `Data/OriginalTreeFiles.tar`) are removed by `cichlid-data migrate --apply --cleanup` once the new copies have been checked, if not done already.
- **Species caveats to remember:** *Etroplus* (authors' models; median protein 351 aa, shorter than all others); *A. centrarchus* (possible retained second haplotypes inflating its species-specific duplications); *H. burtoni* (old scaffold-level annotation, 9% partial models); elephant shark (scaffold-level). Bowfin and *Amphilophus* have no RefSeq annotation; *Pholidichthys* (sister family of cichlids) is unannotated: an NCBI annotation request would be worthwhile.
- **Small proteome remainders:** 6 duplicate-ID genes each in human and zebrafish (one GeneID at two loci with different proteins), kept as `….2`.
- **To confirm before tagging v0.8.0:** the alignment's file name in existing TreeFam2026 and Compara2026 family folders (the contract says `alignment.fa`), and that every existing family and assignment file passes the stricter validation (commands in the delivery notes).
- **Size of `native/`:** OrthoFinder's output for 35 proteomes may be tens of GB. If so, the student keeps the standard results and the files that matter (orthogroups, gene trees, alignments, species tree, logs) and leaves out the bulk working files, or `native` gets excludes. To decide once the first run exists.
- **`run.json` helper:** `catalog_md5` and `coordinator_version` are filled in by each repository by hand today; a helper that builds the record (a candidate for v0.9) would stop them diverging.
- **Which orthogroup level becomes `families/`** (root-level hierarchical orthogroups `N0.HOG…` or OrthoFinder's primary orthogroups `OG…`) is open between the student and `cichlid-trees`; the contract allows either.
- **Stage 7 (combined calls) and the browser** have no schemas yet: the comparison is being redesigned.
- **Version 1.0.0** once the pipeline runs on the coordinator.

## 7. What comes next

1. **The student's repository (`cichlid-orthofinder`)**, from `orthofinder_repo_spec.md`: run OrthoFinder on the 35 proteomes with our species tree (and once with its own), convert the gene trees to the family contract, keep the native outputs; validate with `cichlid-data validate` and `check`.
2. **The pipeline update (`cichlid-trees`):** it takes its paths, transfers, proteomes, inputs and formats from the coordinator; Stages 1 and 2, `newick.py`, `compara_prepare.py` and the download and summary scripts are removed from it (they live here now); the calling and residues code becomes the shared core that reads any method's families (using `schemas.family_rooting`); each run writes a `run.json`; TreeFam and Compara are rerun with the 12 project species.
3. **Tag 1.0.0** once the pipeline runs on the coordinator.
4. **The comparison and browser** are to be redesigned afterwards.

## 8. Documents for the coordinator workspace

| Document | Status |
|---|---|
| `HANDOFF.md` (this file) | Current |
| `project_manifest.md` | Current |
| `coordinator_spec.md` | Current (v0.8.0) |
| `species_final.md` | Current: the agreed species list (and the species tree; `config/species_tree_orthofinder.nwk` is the working copy) |
| `orthofinder_repo_spec.md` | Current (v0.8.0): the instruction sheet for the OrthoFinder repository |
| `species_reference.md`, `species_phylogeny.svg/.png`, `make_species_phylogeny.py` | Current: the three tiers and the figure |
| `orthofinder_species.md` | Background: the draft that led to `species_final.md` (superseded where they differ) |
| `reorganization_proposal.md` | Historical: the proposal the layout was built from (implemented; see the spec) |

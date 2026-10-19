# cichlid-coordinator

File locations, cloud transfers and inputs for the cichlid protein phylogeny project. Every other
project repository (the TreeFam/Compara pipeline, OrthoFinder, analyses, the gene browser) depends on
this package for where data lives and how it moves. Read `docs/HANDOFF.md` first (current state,
procedures, decisions), then `docs/project_manifest.md` and `docs/coordinator_spec.md` (the full
specification).

## Install

```bash
pip install -e .                    # from a clone, for development
# or, in another repository:   pip install "cichlid-coordinator @ git+https://github.com/<org>/cichlid-coordinator@v0.1.0"
cichlid-data init                   # once per machine: writes ~/.config/cichlid-data/config.toml
```

## Use

```python
from cichlid_data import get_layout, get_store
L, S = get_layout(), get_store()
S.ensure(L.proteome("Mzeb"))                       # fetch if missing locally
R = L.results("v1", "Compara")
print(R.orthology("Mzeb").local)                   # Results/v1/Compara/orthology/Mzeb.orthologs.tsv
S.push(R.unit("orthology").name)
```

```bash
cichlid-data where proteome Mzeb
eval "$(cichlid-data env)"
cichlid-data status 'results/*'
cichlid-data pull 'families/TreeFam9/*'
cichlid-data push results/v1/Compara/orthology
cichlid-data migrate              # dry run of moving the inputs to the new layout; --apply to do it
```

## Status (v0.9.2)

Feature-complete against `docs/coordinator_spec.md`: configuration, layout, data units, transfers
(fingerprints, conflict detection, immutable units, locks, transfer log), the layout migration
(applied on mzebra), provenance records, the species catalog, file-format schemas, genome downloads
and local-species preparation, proteomes (Stage 2, with two fixes), all the other downloaders
(annotations, expression, populations, the Compara preparation), and `cichlid-data check`.

New in 0.9.0: releases. One release number (`v1`, `v2`, ...) is shared by every project and defined in
`config/releases.toml` (status, delivered species, input versions, each part's species and species
tree). Results live under `Results/<release>/<part>/` (`Layout.results(release, part)`), and
`cichlid-data release freeze` archives a release's proteomes, catalog and species trees in
`Releases/<release>/`, which never changes; `check` reports drift from it. `run.json` and `done.json`
carry the release. 0.9.2: a family's `stage5/` scratch folder is no longer pushed. 0.9.1 moves `export_seurat.R` (the Daniocell export step) into the package. Earlier (0.8.0): the family contract (`schemas.validate_family`), the OrthoFinder
species tree as configuration, `native/`, and `Catalog.md5()`. Instruction sheets for the repositories
that use the coordinator are in `docs/` (`orthofinder_repo_spec.md`, `trees_repo_spec.md`).

Optional extras: `pip install -e ".[expression]"` (NumPy, SciPy) for `summarize h5ad`;
`".[samples]"` (gspread) for `populations samples --sheet`.

Next: the pipeline uses the coordinator for its paths and transfers.

## Tests

```bash
PYTHONPATH=src:tests python -m unittest discover tests
```



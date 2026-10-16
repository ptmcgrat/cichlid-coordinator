# cichlid-coordinator

File locations, cloud transfers and inputs for the cichlid protein phylogeny project. Every other
project repository (the TreeFam/Compara pipeline, OrthoFinder, analyses, the gene browser) depends on
this package for where data lives and how it moves. Read `docs/project_manifest.md` first, then
`docs/coordinator_spec.md` (the full specification).

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
R = L.results("Compara2026")
print(R.orthology("Mzeb").local)                   # Results/Compara2026/orthology/Mzeb.orthologs.tsv
S.push(R.unit("orthology").name)
```

```bash
cichlid-data where proteome Mzeb
eval "$(cichlid-data env)"
cichlid-data status 'results/*'
cichlid-data pull 'families/TreeFam9/*'
cichlid-data push results/Compara2026/orthology
cichlid-data migrate              # dry run of moving the inputs to the new layout; --apply to do it
```

## Status (v0.6.0)

Done: configuration, layout, data units, transfers (rclone and a local test backend; fingerprints,
conflict detection, immutable units, locks, transfer log), the layout migration (applied on mzebra),
provenance records (`config/sources.toml`, `cichlid-data provenance write|check|show`), the species
catalog (`config/species_sets.toml`, `cichlid-data catalog build|show|missing|validate`), the
file-format schemas (`cichlid_data.schemas`, `cichlid-data validate`), genome downloads (the
pipeline's Stage 1, `cichlid-data download genomes`) and local-species preparation
(`cichlid-data prepare-local`), proteomes (the pipeline's Stage 2, `cichlid-data proteomes`), and
the command line. Depends on `requests`.

Next: the remaining downloaders (annotations, the Compara panel preparation, expression,
population data), then `cichlid-data check`.

## Tests

```bash
PYTHONPATH=src:tests python -m unittest discover tests
```

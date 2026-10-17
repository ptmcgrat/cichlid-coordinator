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

## Status (v0.7.1)

Feature-complete against `docs/coordinator_spec.md`: configuration, layout, data units, transfers
(fingerprints, conflict detection, immutable units, locks, transfer log), the layout migration
(applied on mzebra), provenance records, the species catalog, file-format schemas, genome downloads
and local-species preparation, proteomes (Stage 2, with two fixes), all the other downloaders
(annotations, expression, populations, the Compara preparation), and `cichlid-data check`.

Optional extras: `pip install -e ".[expression]"` (NumPy, SciPy) for `summarize h5ad`;
`".[samples]"` (gspread) for `populations samples --sheet`.

Next: the pipeline uses the coordinator for its paths and transfers.

## Tests

```bash
PYTHONPATH=src:tests python -m unittest discover tests
```

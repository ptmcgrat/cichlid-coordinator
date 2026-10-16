# Proposed data reorganization (draft, October 2026; on hold, to return to)

## Layout

```
TreeFamData/
├── Inputs/
│   ├── Proteins/                      every organism's proteins, one place for all methods
│   │   ├── catalog.json               ALL organisms: source, versions, files, checksums, counts, roles, used_by
│   │   ├── downloads/<prefix>/        as downloaded: genomic.gff, protein.faa, sequence_report.jsonl, source.json
│   │   └── proteomes/<prefix>.faa     one protein per gene (Stage 2), plus <prefix>.genes.tsv
│   ├── Families/                      methods with precomputed families
│   │   ├── TreeFam9/                  hmm/, trees/ (<family>.aln, .nhx), protein_to_species.txt, source.json
│   │   └── Compara116/                dumps/, prepared/ (trees/, search_db.fa, members.tsv,
│   │                                  protein_to_species.txt, family_audit.tsv), panel.tsv, source.json
│   ├── Annotations/                   GO/, PANTHER/, Ensembl_names/ (Stage 6a cache), Legacy2018/; each with source.json
│   ├── Expression/                    HPA/, Zebrahub/, Daniocell/, Bgee/, Brawand/, raw/; each with source.json
│   └── Populations/                   Blumer2025/ (VCFs, Table S1, Table S7, reference/)
├── Results/
│   ├── TreeFam2026/, Compara2026/, OrthoFinder2026/   same internal layout per method:
│   │       assign/<prefix>.assign.tsv, families/<family>/, orthology/<prefix>.orthologs.tsv,
│   │       residues/, names/, run.json (method, version, parameters, inputs used)
│   └── Combined2026/                  orthology/ (Stage 7 merge), site/ (browser)
└── Analyses/                          Variants/ (from Data/Blumer2025 outputs), Tissues/, Duplications/,
                                       Enrichment/, OrthologyClasses/
```

## Key ideas

1. `Inputs/Proteins/catalog.json`: one record per organism (project and reference), generated from species.tsv
   plus the downloads; `used_by` lists the methods that include it.
2. Reference species need their own proteomes for OrthoFinder (human, zebrafish, medaka, gar, ...), downloaded
   into Inputs/Proteins with role "reference".
3. Every method's results share one internal layout; the Compara run's `treefam/` folder becomes `assign/`;
   shared genome data leave the TreeFam run's folder.
4. Every downloaded resource gets a source.json (what, where from, when, version, licence, citation).
5. Analyses move out of results and input folders (Blumer-derived outputs to Analyses/Variants).

## How

1. One path configuration file (layout.py) used by every stage, sync.py and run_pipeline.sh.
2. Migration script with a dry run: moves on mzebra and server-side moves on Dropbox (rclone), writes
   catalog.json and source.json files from the existing download records.
3. A check stage confirming every expected file is in place; then a pipeline run should skip all work.
4. Updated coordinator documentation.

Open choices: Populations/ under Inputs; browser site in Combined2026/site/.

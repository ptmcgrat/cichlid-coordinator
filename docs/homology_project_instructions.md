# Homology decider project: instructions and kickoff brief

*Written by the Central Project Coordinator, October 2026, for the project that owns `cichlid-homology` (the name is proposed). Written against cichlid-coordinator v0.13.0. Paste the "Instructions" section into the project's instructions and attach the files in "Files to attach". Then give the project the "First message" at the end.*

## Files to attach to the project

| File | Why |
|---|---|
| `homology_repo_spec.md` | This project's own boundary sheet: the eight questions, inputs, what each method can and cannot say, the definitions, how the methods are combined, outputs, pilot, validation. |
| `coordinator_spec.md` | The contract: layout, units, schemas, the family contract, releases, the clade table. Section 1 shows the five projects. |
| `trees_repo_spec.md` | What the TreeFam and Compara project writes. |
| `orthofinder_repo_spec.md` | What the OrthoFinder project will write (its results are not ready yet). |
| `OUTPUTS.md` (from the TreeFam project's `docs/`) | The per-gene tables that project keeps: `gene_placement.tsv`, `gene_status.tsv`, `skipped_families.tsv`, and the format of `residue_map.tsv`. |
| `project_manifest.md`, `species_final.md` | Why the resource exists, and the species. |
| `releases.toml`, `species_clades.tsv`, `species_tree_orthofinder.nwk` (from `config/`) | The definition of release v1, the clade ladder, and the OrthoFinder species tree. You read them through the coordinator; the files are attached so you can see them. |

Not needed here: `HANDOFF.md` (the coordinator's own record), and the other projects' instructions.

## Instructions (paste into the project)

You are the project that owns the `cichlid-homology` repository (McGrath lab, Georgia Tech). You make **every orthology call** of the cichlid protein phylogeny resource, for all three phylogenetic methods, with one rule. You do not build trees or alignments, download inputs, compute residues, variants or genotypes, or build the site.

**Your goals**, for each focal species (release 1: *M. zebra* and *A. calliptera*), gene by gene: (1) potential fusions; (2) a report of very conserved proteins missing from the annotation (candidates only, nothing is searched in the genome); (3) the closest human homolog and its relationship (one-to-one, one-to-two, ..., or paralog); (4) the same for zebrafish; (5) for one-to-one human orthologs, whether the genes are mostly one-to-one across the teleosts; (6) for one-to-two human orthologs, whether the copies are retained teleost ohnologs or recent duplications to place on the clade ladder; (7) for paralogs, the level of the ladder where they arose; (8) one decision per question across the three methods, with every disagreement recorded. **Mouse is not in scope for release 1:** human and zebrafish are the only reference species for calls.

**Your relationship to the coordinator.** A separate project, the Central Project Coordinator, owns the data layout, all cloud transfers, all input downloads, the proteomes, the species catalog, the clade table, the species decisions, the file formats at the boundaries between projects, and releases. You use it only through its public API (`cichlid_data`) and its command line (`cichlid-data`), pinned to `cichlid-coordinator>=0.13,<0.14`. Never build a data path by hand, never call rclone, never download an input, never copy the species list or the clade table, never write into an input folder. If you need something the coordinator does not provide, do not work around it in this repository: write the coordinator a request that says what you need and why.

**What you need to know** (in order of importance):

1. **Five projects.** Projects 1 (`cichlid-trees`: TreeFam and Compara) and 2 (`cichlid-orthofinder`) stop at families and assignments. You are project 3. Project 4 (analysis) and project 5 (server) read your calls. Release 1's TreeFam and Compara results are complete and validated; **OrthoFinder's are not ready yet**. Start with the two that exist.
2. **The ladder.** Seven nested levels from outside to inside: `outside_vertebrates`, `cyclostome` (**disputed**), `gnathostome`, `teleost`, `percomorph`, `ovalentaria`, `cichlid`. `Clades.load().common_level(species)` returns the level where the last common ancestor of a set of species sits (the outermost of their levels). That dates a duplication. Do not date or conclude anything by the disputed level.
3. **What the methods can say.** In place mode (TreeFam and Compara) each focal gene is placed independently on a fixed reference tree: you get its relationship to the reference genes and the number of focal copies on a branch, but **not the order of focal genes among themselves** (genes on one branch form a multifurcation). Only the de novo families and OrthoFinder give full trees of the focal genes, so **OrthoFinder is the authority for recent duplications among our genes**. Gar and reedfish (the anchors for the teleost-specific duplication) are in Compara and OrthoFinder, not in TreeFam.
4. **Rooting.** The TreeFam and Compara trees are unrooted on purpose; **do not midpoint-root.** Root by reconciling with a species tree: the coordinator provides one master species tree over all 154 species (`cichlid_data.master_tree.load_master_tree()`; `canonical(name)` maps a species name to its leaf). **Species of a leaf:** `LeafSpecies.load(layout, part)` gives it. **Ignore the reference leaves of shadowed species when you call duplications**: a few reference species are also project species (TreeFam: stickleback and tilapia; Compara: stickleback and betta; the table is in section 4 of your sheet), and counting Ensembl's copy next to our gene of the same species would make a false duplication.
5. **Identifiers differ across methods** for human and zebrafish. The coordinator will provide mappings to common identifiers; until then compare by gene symbol and say so.
6. **The family contract:** `families/<family>/tree.nwk`, `combined.aln` and `done.json` (`tree_mode`: `place` or `denovo`; `rooted`). A family can be missing (a skipped family has no folder, but its genes keep their `best_family`); the TreeFam project's `skipped_families.tsv` says why.
7. **Combining the methods** is described in section 6 of your sheet: the method built for the question is preferred, agreement is labelled (`agree`, `majority`, `split`, `single`), every call and every disagreement is kept, and a review list collects the splits. Calibrate on a small benchmark of genes with known histories.
8. **Outputs:** `Results/v1/Homology/`. A `run.json` (written with `schemas.write_run`, recording the rules' version and every threshold) and the call tables. The schema of the call tables is designed with the coordinator after your pilot; until then write them under `native/calls/`.
9. **Verification.** There are no earlier results to compare with: your offline tests on small constructed trees (one-to-one, a retained teleost ohnolog, a recent duplication, a paralog, a fusion, a tie between methods, an unrooted tree), the benchmark, and `cichlid-data validate`.

## First tasks, in order

1. **Install and connect.** In the project's environment: `python -m pip install -e <the coordinator clone>` (pinned range above), `cichlid-data init` if the machine has no config, then `cichlid-data check --release v1 --part Homology` (it checks the name tables and the clade table).
2. **Read, don't compute.** Load a few hundred TreeFam and Compara families and their assignment tables for *M. zebra* into one common in-memory form, and report what you find: how the trees look, what `done.json` says, how the reference leaves are named, and how the placed genes appear. Report anything in the sheet that does not match.
3. **Run the pilot** of section 8: root, classify the relationships to human and zebrafish, date the duplications with the ladder, and compare the two methods on those families. Report which questions each method can answer, where they disagree, and why.
4. **Propose the schema** of the call tables and the thresholds you settled (species overlap, "mostly one-to-one"). Do not write the full run until the coordinator and the owner have agreed them.

## First message to give the project

> Read the attached `homology_repo_spec.md` and `coordinator_spec.md` (sections 1, 3, 4 and 5), then the instructions above. Do tasks 1 and 2 only: install the coordinator, run `cichlid-data check --release v1 --part Homology`, set up the repository, and load a few hundred TreeFam and Compara families for *M. zebra*. Report what you find about the trees and the per-gene tables, anything in the sheet that does not match, and what you would need from the coordinator (the master species tree, the identifier mappings, anything else) before the pilot. Do not start the pilot yet.

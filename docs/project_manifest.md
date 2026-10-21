# Cichlid protein phylogenies: project manifest

*Central Project Coordinator knowledge file. October 2026. Read this first.*

## Purpose

Most of what is known about any gene's function comes from human, mouse and zebrafish. Cichlids, and other non-traditional models such as betta and stickleback, are powerful systems for studying how behavior, development and form evolve, but their genetic tools are limited. **Orthology is the bridge:** it decides whether knowledge from model organisms can be applied to a gene in our species, and how confidently. And because many fish genes exist in two copies after the teleost genome duplication, knowing which duplicate copies were kept, lost, or newly made in cichlids points to where fish- and cichlid-specific biology may lie.

This resource rebuilds and extends the 2018 cichlid protein phylogenies (Ramakrishnan Varadarajan et al., *BMC Evolutionary Biology* 18:1).

## The two aims

### Aim 1: human and zebrafish orthology for every gene of the project species

For each gene: its best human and best zebrafish homolog, and the relationship:

- **one-to-one ortholog:** one copy in our species for one human gene. The most reliable case for transferring function.
- **one-to-two ortholog:** two copies for one human gene.
- **one-to-many ortholog,** **paralog** (the gene's own human ortholog is missing; the closest human gene is a relative from an older duplication), or **no homolog**.

Each call carries a confidence level, the methods that support it, and the tree evidence behind it.

### Aim 2: what happened to duplicated genes in teleost evolution

- **One-to-two pairs, by origin.** Whether the duplication is **shared across teleosts** (the teleost genome duplication, 3R) or **more specific,** and if specific, the clade where it arose: Ovalentaria, percomorphs, or within cichlids, down to a branch of the cichlid tree.
- **Retention and loss across teleosts.** For each gene, how many copies other teleosts kept; and where the second copy was lost.
- **Cichlid-specific duplications,** placed on the cichlid species tree, rooted with the earliest-branching cichlids and checked against the closest non-cichlid relatives.
- **Paralogs, by where the human ortholog went missing:** lost in mammals, lost in tetrapods, ray-finned fish only, or teleost-only.

### Further aims

Gene losses in the cichlid lineage (confirmed against the genome sequence); guidance on transferring human and zebrafish annotations given each relationship; consistent gene naming; links to population variation, expression and gene sets.

## Species

The final list, with sources and which method uses each, is in `species_final.md`. In brief:

1. **Project species,** which get orthology calls: ten cichlids, from *Etroplus suratensis* (the earliest-branching cichlid lineage) to the Lake Malawi radiation, plus **betta and stickleback** as non-traditional models.
2. **The cichlids' relatives** (Ovalentaria: medaka, killifishes, guppy, platyfish, clownfish; ideally the sister family *Pholidichthys*), to tell cichlid-specific duplications from older ones.
3. **Reference species:** teleosts across the lineages that split after the genome duplication (from eel and bonytongues to percomorphs), fish from before it (gar, reedfish, ideally bowfin), other vertebrates including human and mouse, and an outgroup.

Lineages with their own extra genome duplications (salmon, carp, sturgeons) are deliberately excluded.

## The three phylogenetic approaches

We use three independent approaches, and compare and combine them rather than trusting any one. All three use the same rule to read relationships off the trees (species overlap: a node with the same species on both sides is a duplication), so differences between them reflect the trees, not the calling logic.

| | TreeFam 9 | Ensembl Compara 116 | OrthoFinder |
|---|---|---|---|
| **What it is** | Curated gene families with HMMs and trees from 109 genomes (2013); the method of the 2018 resource | Ensembl's current gene trees, from about 200 genomes; we keep a 57-species panel | Families and trees built from scratch from our own species set |
| **How our genes enter** | Assigned to a family by HMM, then **placed** onto its fixed tree | Assigned by sequence search, then **placed** onto the fixed tree | **Built into** the trees with every other species |
| **Fish sampling** | Thin: 8 teleosts, 3 ovalentarians, no gar | Dense: about 40 teleosts, 9 ovalentarians, gar and reedfish | Chosen by us: 35 species, including *Etroplus*, eel and gar |
| **Strengths** | Broad, stable families: about 96% of genes assigned. Continuity with 2018. Deep outgroups for very old families | Current annotations. Dense fish sampling resolves teleost duplicates best, giving the best zebrafish one-to-one calls. Trees built with the species tree in mind | Every species treated alike, project species included. We choose the outgroups that define cichlid-specific events. Duplications mapped directly onto our species tree. Hierarchical families at every level of the tree. Independent of families designed around human or other model species |
| **Weaknesses** | Old data and annotations. Too few fish to resolve the teleost genome duplication well, so more duplications look cichlid-specific than really are. Its reference set can't be changed | Some families are split into separate fish and tetrapod trees, hiding about 2,300 human orthologs per species. Reference species limited to what Ensembl annotates. Trees change between releases | Most sensitive to annotation quality: missing or fragmented gene models in any species become apparent losses or duplications. Costly, scaling roughly with the square of the number of species. Needs a fully resolved species tree, so debated relationships must be chosen. Families are clustered automatically, without curation |
| **Best for** | A broad first call for nearly every gene; comparison with 2018 | Zebrafish relationships; teleost duplicates and their retention | Dating duplications on our own tree; an independent third opinion on every call |

**Placement versus building.** TreeFam and Compara keep their reference trees fixed and add our genes onto them. That is fast and stable, and our genes can't distort the trees, but they also can't improve them: if the reference tree is wrong for a family, our genes inherit the error. OrthoFinder builds every tree with our genes included, so our species shape their own trees; the price is that their annotation errors shape them too.

**Combining them.** Where TreeFam and Compara both call a human ortholog, they pick the same gene 97% of the time. The combined calls use each method where it is strongest: TreeFam fills Compara's split fish-only trees, Compara supplies the zebrafish calls. OrthoFinder will add a third independent call, so each relationship can carry "supported by 3 of 3".

## Principles

- **Show the evidence:** every call records the tree node that decided it; the browser shows the tree.
- **Say how sure we are:** confidence levels, agreement between methods, flags for weak gene models.
- **Track every source:** every genome, annotation and dataset has a recorded origin, version, date and licence.
- **Reproducible and resumable:** the pipeline reruns to the same results, and each stage skips work that is still valid.
- **Annotation quality first:** gene-model errors are the largest single source of error, so species are chosen and flagged with that in mind.

## Products

- **Gene browser:** gene pages with orthology calls, trees, alignments, variants, expression and gene sets.
- **Tables:** orthology calls per species; one-to-two pairs with their origin; cichlid-specific duplications by branch; paralogs by where the human ortholog went missing.
- **Gene sets:** inversion candidates, Y-linked genes, cichlid-specific duplications by branch.
- **Versioned releases** for citation, and documentation for the lab and the community.

## What we have learned so far

- **TreeFam and Compara agree** on the human ortholog 97% of the time where both make a call.
- **One-to-one orthologs are robust:** for cichlid one-to-one genes, other teleosts have exactly one copy 85% of the time.
- **One-to-two pairs fall into two clear groups:** about 60% are shared by nearly all teleosts (the teleost genome duplication), and the rest by few or none. Of those, most are pairs kept only in cichlids from older duplications, and a smaller set (about 190) arose within cichlids, mostly on deep cichlid branches.
- **Cichlid-specific duplications** on the lake-radiation branches reproduce Brawand et al. (2014) closely.

## Known limitations

- Gene-model quality varies between species and limits every call.
- The youngest duplications are the hardest to place, and assembly artifacts can mimic them.
- Some relationships in the species tree are still debated, which affects where duplications are placed.
- The cichlids' sister family (*Pholidichthys*) and the Madagascar cichlids (Ptychochrominae) have no annotated genomes yet.

# Final species list

*Central Project Coordinator knowledge file. Agreed October 2026. Supersedes the draft lists in `orthofinder_species.md` where they differ.*

Aims this list serves: (1) human and zebrafish orthology relationships for project species; (2) how duplicates from the teleost genome duplication (3R) were kept or lost across teleosts.

**How each method uses species.** TreeFam 9 and Ensembl Compara bring their own reference species in precomputed trees; we only choose which **project species** are placed onto them (and, for Compara, which reference species to keep: `compara_panel.tsv`). OrthoFinder builds families from scratch, so it uses exactly the proteomes we give it.

**Status:** *downloaded* = proteome in hand; *to download* = NCBI RefSeq, fetched by Stage 1; *in resource* = present inside the TreeFam or Compara data, no separate download. As of October 2026 all 35 OrthoFinder proteomes are built and in the cloud (`cichlid-data catalog missing --method OrthoFinder2026` reports none missing).

Columns: **TF** = TreeFam 9, **C** = Compara 116, **OF** = OrthoFinder. "query" = our proteins placed on the method's trees, with calls made; "ref" = a reference species in the method's trees.

## 1. Project species (orthology calls in every method)

| Prefix | Species | Source | Status | TF | C | OF |
|---|---|---|---|---|---|---|
| Esur | *Etroplus suratensis* | Authors' proteins; positions from figshare 26303968 | downloaded (pipeline not yet run) | query | query | yes |
| Acen | *Archocentrus centrarchus* | RefSeq GCF_007364275.1 | downloaded | query | query | yes |
| Pmar | *Pelmatolapia mariae* | RefSeq GCF_036321145.2 | downloaded | query | query | yes |
| Onil | *Oreochromis niloticus* | RefSeq GCF_001858045.2 | downloaded | query (TreeFam also has an older tilapia) | query | yes |
| Nbri | *Neolamprologus brichardi* | RefSeq GCF_056152115.1 | downloaded | query | query | yes |
| Pmic | *Perissodus microlepis* | RefSeq GCF_051167445.1 | downloaded | query | query | yes |
| Sdia | *Simochromis diagramma* | RefSeq GCF_900408965.1 | downloaded | query | query | yes |
| Hbur | *Haplochromis burtoni* | RefSeq GCF_018398535.1 | downloaded | query | query | yes |
| Mzeb | *Maylandia zebra* | RefSeq GCF_041146795.1 | downloaded | query | query | yes |
| Acal | *Astatotilapia calliptera* | RefSeq GCF_964374335.1 | downloaded | query | query | yes |
| Bspl | *Betta splendens* (betta) | RefSeq (accession in `config/species.tsv`) | downloaded (new project species) | query | query (Ensembl copy becomes validation) | yes |
| Gacu | *Gasterosteus aculeatus* (stickleback) | RefSeq (accession in `config/species.tsv`) | downloaded (new project species) | query (TreeFam also has it as ref) | query (Ensembl copy becomes validation) | yes |

## 2. Cichlids' relatives, for placing duplications

| Species | Source | TF | C | OF | Why |
|---|---|---|---|---|---|
| *Amphilophus citrinellus* | Ensembl | – | ref | candidate (if annotation is good) | Second Neotropical cichlid |
| *Pundamilia nyererei* | Ensembl | – | ref | – | Lake Victoria; fragmented assembly |
| Medaka (*Oryzias latipes*) | RefSeq (OF) | ref | ref | yes | Atherinomorpha |
| Turquoise killifish (*Nothobranchius furzeri*) | RefSeq (OF) | – | ref | yes | Cyprinodontiformes |
| Mummichog (*Fundulus heteroclitus*) | RefSeq (OF) | – | ref | yes | Cyprinodontiformes |
| Guppy (*Poecilia reticulata*) | RefSeq (OF) | – | ref | yes | Cyprinodontiformes |
| Platyfish (*Xiphophorus maculatus*) | RefSeq (OF) | ref | ref | yes | Cyprinodontiformes |
| Clownfish (*Amphiprion ocellaris*) | RefSeq (OF) | – | ref | yes | Damselfishes |
| *Acanthochromis polyacanthus*, *Cyprinodon variegatus*, *Kryptolebias marmoratus* | Ensembl | – | ref | – | Extra ovalentarians in the Compara trees |
| *Pholidichthys leucotaenia* | Unannotated (GCA_020510985.1) | – | – | when annotated | Sister family of cichlids |

## 3. Reference species

### Teleosts, across the lineages that diverged after the teleost genome duplication

| Species | Group | TF | C | OF |
|---|---|---|---|---|
| Fugu (*Takifugu rubripes*) | Percomorphs | ref | ref | yes |
| 13 other percomorphs (sea bass, seabream, wrasse, turbot, seahorse and others) | Percomorphs | some | ref | – |
| Cod (*Gadus morhua*) | Gadiformes | ref | ref | yes |
| Northern pike (*Esox lucius*) | Esociformes | – | ref | yes |
| Soldierfish (*Myripristis murdjan*) | Holocentriformes | – | ref | – |
| Atlantic herring (*Clupea harengus*) | Clupeiformes | – | ref | yes |
| Zebrafish (*Danio rerio*) | Cypriniformes | ref (calls) | ref (calls) | yes (calls) |
| Mexican cavefish (*Astyanax mexicanus*) | Characiformes | – | ref | yes |
| Catfish, electric eel, denticle herring | Otomorpha | – | ref | – |
| Elephantfish (*Paramormyrops kingsleyae*) | Bonytongues | – | ref | yes |
| European eel (*Anguilla anguilla*) | Elopomorpha | – | – | **yes (added)** |

### Fish before the teleost genome duplication

| Species | TF | C | OF |
|---|---|---|---|
| Spotted gar (*Lepisosteus oculatus*) | – | ref | yes |
| Bowfin (*Amia calva*) | – | – | yes, if a RefSeq annotation exists |
| Reedfish (*Erpetoichthys calabaricus*) | – | ref | yes |

### Other vertebrates

| Species | TF | C | OF |
|---|---|---|---|
| Human | ref (calls) | ref (calls) | yes (calls) |
| Mouse | ref (calls) | ref (calls) | yes (calls) |
| Chicken | ref | ref | yes |
| Western clawed frog (*Xenopus tropicalis*) | ref | ref | yes |
| Coelacanth | ref | ref | yes |
| Elephant shark (*Callorhinchus milii*) | – | ref | yes |
| Sea lamprey | ref | ref | **no (dropped)** |
| Other mammals, birds and reptiles | many | 10 | – |

### Outgroups

| Species | TF | C | OF |
|---|---|---|---|
| Amphioxus (*Branchiostoma floridae*) | – | – | yes |
| *Ciona intestinalis* | ref | ref | **no (dropped)** |
| Fly, nematode, yeast | ref | ref | – |

## OrthoFinder set: 35 species

12 project species (10 cichlids, betta, stickleback) + 6 ovalentarians + 8 other teleosts (fugu, cod, pike, herring, zebrafish, cavefish, elephantfish, eel) + gar and reedfish + human, mouse, chicken, frog, coelacanth, elephant shark + amphioxus. Plus bowfin and *Amphilophus* if their annotations check out (37).

Species tree for `-s` (choices marked in `orthofinder_species.md`; the eel's position relative to the bonytongues is also debated, so test both). **The file `config/species_tree_orthofinder.nwk` is the one the pipeline uses**: it is checked against the catalog's OrthoFinder species (same set, rooted, fully resolved) and a test keeps it identical to the tree below. To change the tree, change both together.

```
(Bflo,(Cmil,((Lcha,(Xtro,(Ggal,(Mmus,Hsap)))),(Ecal,(Locu,(Aang,(Pkin,((Char,(Drer,Amex)),(Eluc,(Gmor,((Gacu,Trub),(Bspl,(((Esur,(Acen,(Pmar,(Onil,(Nbri,(Pmic,(Sdia,(Hbur,(Mzeb,Acal))))))))),Aoce),(Olat,(Nfur,(Fhet,(Pret,Xmac)))))))))))))))));
```

(Bowfin, if added, joins gar: `(Locu,Amca)`.)

## Changes from the drafts

1. Betta and stickleback promoted to project species (calls in all three methods).
2. European eel added to OrthoFinder: with the bonytongues, it brackets the teleost root, so duplications can be placed at the teleost genome duplication itself.
3. Sea lamprey and *Ciona* dropped from OrthoFinder: they don't serve either aim, and both complicate families. They remain in TreeFam and Compara, whose trees are fixed.
4. The OrthoFinder species tree became configuration (`config/species_tree_orthofinder.nwk`), validated against the species list (v0.8.0).

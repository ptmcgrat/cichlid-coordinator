# OrthoFinder species set (draft)

*Central Project Coordinator knowledge file. Draft, October 15, 2026. Companion to `species_reference.md`.*

## Principles

- **37 core species,** following the three tiers in `species_reference.md`: project cichlids, the ovalentarian relatives needed to define cichlid-specific duplications, and reference species spanning fish and vertebrate evolution.
- **Annotation quality first.** OrthoFinder builds families from our proteomes alone, so a missing gene looks like a loss and a fragmented or duplicated gene model looks like a duplication. Prefer chromosome-level assemblies with **NCBI RefSeq** annotation: the same annotation pipeline as our cichlids (comparable completeness), and our Stages 1–2 fetch them and pick one protein per gene unchanged.
- **No lineage-specific genome duplications** (salmon, carp, sturgeons): they make orthogroups messy without answering our questions.
- **No distant invertebrates or fungi** (fly, worm, yeast). They would merge ancient families and raise the cost, for no gain at the levels we study. Amphioxus and a tunicate are enough to root the vertebrates.
- **Run in alignment mode** (`-M msa`, MAFFT + IQ-TREE or FastTree), with **our species tree** supplied (`-s`), and once with OrthoFinder's own inferred tree for comparison.

Accessions are not listed here: look them up with the command at the end, which also confirms that each RefSeq annotation exists.

## Core species (37)

### Tier 1: project cichlids (10)

| Prefix | Species | Source |
|---|---|---|
| Esur | *Etroplus suratensis* | Authors' annotation (`source=local`) |
| Acen | *Archocentrus centrarchus* | RefSeq (GCF_007364275.1). Possible retained-haplotype duplicates: check before relying on its species-specific counts |
| Pmar | *Pelmatolapia mariae* | RefSeq |
| Onil | *Oreochromis niloticus* | RefSeq |
| Nbri | *Neolamprologus brichardi* | RefSeq |
| Pmic | *Perissodus microlepis* | RefSeq |
| Sdia | *Simochromis diagramma* | RefSeq (scaffold-level, older) |
| Hbur | *Haplochromis burtoni* | RefSeq (scaffold-level, older; fragmentary models) |
| Mzeb | *Maylandia zebra* | RefSeq |
| Acal | *Astatotilapia calliptera* | RefSeq |

### Tier 2: other ovalentarians (6)

These define what is cichlid-specific: a duplication counts as cichlid-specific only if none of them share it.

| Prefix | Species | Group | Why |
|---|---|---|---|
| Olat | *Oryzias latipes* (medaka) | Beloniformes | Model fish; Atherinomorpha |
| Nfur | *Nothobranchius furzeri* (turquoise killifish) | Cyprinodontiformes | Model for aging |
| Fhet | *Fundulus heteroclitus* (mummichog) | Cyprinodontiformes | Model for physiology |
| Pret | *Poecilia reticulata* (guppy) | Cyprinodontiformes | Model for sexual selection and color |
| Xmac | *Xiphophorus maculatus* (platyfish) | Cyprinodontiformes | Model for melanoma, sex determination |
| Aoce | *Amphiprion ocellaris* (clownfish) | Pomacentridae | Damselfishes: a second ovalentarian family, outside Atherinomorpha |

### Tier 3: reference species (21)

| Prefix | Species | Group | Why |
|---|---|---|---|
| Gacu | *Gasterosteus aculeatus* (stickleback) | Percomorphs: Eupercaria | Model for repeated adaptation |
| Trub | *Takifugu rubripes* (fugu) | Percomorphs: Eupercaria | Compact, well-studied genome |
| Bspl | *Betta splendens* (betta) | Percomorphs: Anabantaria | Model for aggression; a third percomorph series |
| Gmor | *Gadus morhua* (cod) | Gadiformes | Teleost outside percomorphs |
| Eluc | *Esox lucius* (pike) | Esociformes | Salmonid lineage without the salmonid duplication |
| Char | *Clupea harengus* (Atlantic herring) | Clupeiformes | Otomorpha, outside the carps and catfishes |
| Drer | *Danio rerio* (zebrafish) | Cypriniformes | Main fish model; ortholog calls |
| Amex | *Astyanax mexicanus* (cavefish) | Characiformes | Model for eye loss; second otophysan |
| Pkin | *Paramormyrops kingsleyae* (elephantfish) | Osteoglossomorpha | Earliest-branching teleosts in our set |
| Locu | *Lepisosteus oculatus* (spotted gar) | Holostei | Pre-teleost-duplication state |
| Amca | *Amia calva* (bowfin) | Holostei | Second holostean; RefSeq annotation: to check |
| Ecal | *Erpetoichthys calabaricus* (reedfish) | Polypteriformes | Earliest-branching ray-finned fish |
| Lcha | *Latimeria chalumnae* (coelacanth) | Sarcopterygii | Lobe-finned fish |
| Xtro | *Xenopus tropicalis* (western clawed frog) | Amphibians | Diploid frog (not *X. laevis*, which is tetraploid) |
| Ggal | *Gallus gallus* (chicken) | Birds | |
| Mmus | *Mus musculus* (mouse) | Mammals | Ortholog calls |
| Hsap | *Homo sapiens* (human) | Mammals | Ortholog calls |
| Cmil | *Callorhinchus milii* (elephant shark) | Cartilaginous fishes | Jawed-vertebrate outgroup to bony fish |
| Pmrn | *Petromyzon marinus* (sea lamprey) | Cyclostomes | Jawless vertebrate; 1R/2R |
| Bflo | *Branchiostoma floridae* (amphioxus) | Cephalochordates | Root: invertebrate chordate before 1R/2R |
| Cint | *Ciona intestinalis* (tunicate) | Tunicates | Closest invertebrate relatives of vertebrates |

Prefix note: *Amia calva* and *Petromyzon marinus* would naturally take `Acal` and `Pmar`, which belong to *A. calliptera* and *P. mariae*; hence `Amca` and `Pmrn`.

## Candidates to add later (with `--assign`)

The core run must use alignment mode for these to be added afterwards.

| Species | Why | Caveat |
|---|---|---|
| *Amphilophus citrinellus* | Second Neotropical cichlid: tests the *A. centrarchus* duplications | Ensembl only (2016 assembly); check annotation quality |
| *Pundamilia nyererei* | Lake Victoria haplochromines | Fragmented 2014 assembly |
| *Pholidichthys leucotaenia* | Sister family of cichlids: the most important missing species | Unannotated (needs an NCBI annotation request, or our own gene models) |
| Ptychochrominae, e.g. *Paratilapia polleni* | Second-earliest cichlid branch | Low-coverage draft only |
| *Anguilla anguilla* (eel) | Elopomorpha: resolves the base of teleosts | Optional |
| *Acanthochromis polyacanthus* | Second damselfish | Optional |

## Species tree for `-s`

OrthoFinder needs a rooted, fully resolved tree. Three choices are made here where relationships are still debated; they are marked so they can be changed:

1. **Percomorph series:** (Eupercaria, (Anabantaria, Ovalentaria)).
2. **Within Ovalentaria:** ((cichlids, damselfishes), Atherinomorpha).
3. **Base of teleosts:** bonytongues first (eels not included).

```
(Bflo,(Cint,(Pmrn,(Cmil,((Lcha,(Xtro,(Ggal,(Mmus,Hsap)))),(Ecal,((Locu,Amca),(Pkin,((Char,(Drer,Amex)),(Eluc,(Gmor,((Gacu,Trub),(Bspl,(((Esur,(Acen,(Pmar,(Onil,(Nbri,(Pmic,(Sdia,(Hbur,(Mzeb,Acal))))))))),Aoce),(Olat,(Nfur,(Fhet,(Pret,Xmac))))))))))))))))));
```

Running OrthoFinder a second time with its own species tree (STAG and STRIDE) shows whether these choices affect the duplication counts.

## Looking up accessions

On `mzebra`, with NCBI's `datasets` tool (already used by Stage 1), for each species:

```bash
for sp in "Oryzias latipes" "Nothobranchius furzeri" "Fundulus heteroclitus" "Poecilia reticulata" \
          "Xiphophorus maculatus" "Amphiprion ocellaris" "Gasterosteus aculeatus" "Takifugu rubripes" \
          "Betta splendens" "Gadus morhua" "Esox lucius" "Clupea harengus" "Danio rerio" "Astyanax mexicanus" \
          "Paramormyrops kingsleyae" "Lepisosteus oculatus" "Amia calva" "Erpetoichthys calabaricus" \
          "Latimeria chalumnae" "Xenopus tropicalis" "Gallus gallus" "Mus musculus" "Homo sapiens" \
          "Callorhinchus milii" "Petromyzon marinus" "Branchiostoma floridae" "Ciona intestinalis"; do
  datasets summary genome taxon "$sp" --reference --assembly-source RefSeq --as-json-lines 2>/dev/null | \
    python3 -c "import json,sys
for l in sys.stdin:
    r=json.loads(l); a=r.get('assembly_info',{}); n=r.get('annotation_info',{})
    print(f\"{r['organism']['organism_name']:<28} {r['accession']:<18} {a.get('assembly_level',''):<12} \"
          f\"{a.get('assembly_name',''):<22} annotated: {n.get('release_date','none')}\")" || echo "$sp: none"
done
```

Each line gives the RefSeq accession, assembly level, assembly name and annotation date. A species with no line has no RefSeq reference annotation, and needs another source (Ensembl, or the authors') or should drop to the candidates list.

## Run outline

1. Add the 27 new species to `species.tsv` (`source=ncbi`), run Stages 1–2: one protein per gene, as for the cichlids.
2. Put the 37 proteomes in one folder (in the reorganized layout: `Inputs/Proteins/proteomes/`).
3. Run, for example: `orthofinder -f <proteomes> -M msa -T iqtree -s species_tree.nwk -t 32 -a 8 -o Results/OrthoFinder2026`.
4. Convert its gene trees and hierarchical orthogroups into our family format, so Stages 6–8 make calls exactly as for TreeFam and Compara; and keep its own duplication mapping as an independent duplication analysis.

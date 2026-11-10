# Curated functional and independent-UPR readouts

The primary measurements are fixed mean-expression readouts for ER protein folding, ER-associated degradation, ER targeting and integrated stress response signaling. Independent UPR definitions use the curated Reactome parent and a GO parent term. Expression of these genes does not measure biochemical pathway activity. The same memberships separately annotate perturbation-target functions in `target_functional_annotations.tsv`; these descriptive target strata are small or absent and are not alternative readout definitions.

## Definitions and provenance

| Readout | Authoritative definition | Canonical / measured genes |
|---|---|---|
| ER folding | [GO:0034975](https://amigo.geneontology.org/amigo/term/GO:0034975) | 13 / 9 |
| ERAD | [GO:0036503](https://amigo.geneontology.org/amigo/term/GO:0036503) | 96 / 67 |
| ER targeting | [GO:0045047](https://amigo.geneontology.org/amigo/term/GO:0045047) | 38 / 34 |
| ISR | [GO:0140467](https://amigo.geneontology.org/amigo/term/GO:0140467) | 37 / 20 |
| Reactome UPR | [R-HSA-381119.4](https://reactome.org/content/detail/R-HSA-381119) | 92 / 69 |
| GO ER-UPR | [GO:0030968](https://amigo.geneontology.org/amigo/term/GO:0030968) | 64 / 41 |

ER targeting specifically covers signal-directed co/post-translational protein targeting/import into the ER; general ER-to-Golgi vesicular trafficking is outside its scope. ISR annotations are propagated from its kinase-branch descendants; the annotated gene set is not a manually selected list of kinase or ATF4 targets.

`authoritative_raw/` preserves the complete API responses and official HGNC symbol snapshot. QuickGO used ontology release 2026-09-27 and annotation snapshot 2026-07-28 10:08. Human UniProtKB protein gene products, including electronic evidence, were propagated through `is_a` and `part_of`; `NOT`, RNAcentral and ComplexPortal annotations do not define individual protein-gene members. Reactome database98 participants were expanded to UniProt protein reference entities. The HGNC complete-set snapshot was last modified on2026-11-05 and downloaded on2026-11-06; its HTTP provenance and SHA256 are retained. HGNC UniProt IDs, approved symbols and unique previous/alias symbols resolve nomenclature; unresolved or ambiguous records remain visible in the mapping tables. No canonical gene maps to multiple columns in the fixed assay universe.

`locked_memberships.tsv` fixes all sixteen definitions: four full functional modules, their Hallmark-disjoint and between-functional-module-unique versions, and full/disjoint versions of both UPR parents. Disjoint definitions remove the complete original113-member Hallmark snapshot before intersecting measured genes. Canonical mapping represents94 of those genes in the assay universe: the original90 plus the older assay symbols H2AFX, IARS, TARS and TTC37. These four do not belong to the new curated real readouts, but are also excluded from random-module pools. Module-unique versions remove all members shared with the other three full functional modules. No genes are added to replace excluded or unmeasured members. `hierarchy_and_overlap.tsv` reports gene-set overlap and set containment, rather than asserting that every set is a distinct biological pathway; `curated_direct_hierarchy.tsv` preserves direct GO relations and Reactome subpathways. The original measured Hallmark90 is retained unchanged as a reference.

Memberships were locked at 2026-11-06 05:37:18 UTC, before their response scoring, in `membership_lock_manifest.tsv`. The source-only random-module lock is `scoring_lock_manifest.tsv`. Its eight families each have100 size-matched modules using the same9,775 source non-targeting cells and fixed5×5 expression/detection bins. All30,200 member draws matched exact bins. The tested readout and full Hallmark overlap are excluded from each sampling pool; all seeds and draws are retained. No validation responses enter matching.

## Measurements and outputs

All four studies retain the existing6,780-gene CP10k denominator. The metadata-only four-study intersection has6,783 genes; the extra three do not change primary normalization. The fixed panels are438 K562/RPE1 genes,30 Adamson genes and19/19/17 Frangieh genes in Control/IFNγ/co-culture. K562 extraction also retains the source union needed for the older external panels, giving467 source genes.

Primary real readouts recompute every log term using the library excluding the targeted gene, when it belongs to the denominator, and remove that gene from readout numerators. Target-specific controls undergo the identical operation before global SD and control means are computed. Targets measurable only outside the denominator are not subtracted from that denominator. A separate fixed-library sensitivity uses exact own-member numerator removal for both real and random readouts. These random modules calibrate that fixed-library sensitivity; they do not calibrate the exact-denominator primary analysis.

K562/RPE1 retain matched-GEM controls and target-cell-weighted contrasts with global target-specific control SDs. Adamson retains the original global63(mod) anchor; barcode capture lanes are technical groups, not biological replicates. Frangieh uses pooled non-targeting cells within each condition. No gene or module is removed for response magnitude. Empty post-exclusion readouts or nonpositive control SD are unestimable; fewer than three remaining members is flagged as a few-transcript proxy.

`scored/{study}_gene_effects.tsv.gz` contains all real readouts in both normalization modes and all800 null readouts in fixed-library mode. Companion files contain real-readout GEM contrasts, member-level expression contrasts, target-specific control moments, expression/detection covariates and measurement accounting. `{study}_cell_cache.rds` retains real per-cell scores, raw real-readout counts, target counts, control scores and null sufficient sums; it avoids an800-null-by-all-cell matrix. `f_positive` uses groups with at least five target cells; the primary mean contrast retains all observed groups.

## Reproduction

Use the preserved R environment records and packages `data.table`, `Matrix`, `hdf5r` and `digest`; authoritative downloads use Python's standard library and the Reactome REST utility. `prepare_metadata.R` reads feature annotations only. `fetch_authoritative.py`, `freeze_memberships.py` and `prepare_nulls.R` reconstruct definitions and source-only controls. `extract_counts_scores.R PROJECT STUDY` streams each original matrix once; `summarise_scores.R PROJECT STUDY` derives moments from the compact cache. Studies are `K562`, `RPE1`, `Adamson` and `Frangieh`; `PROJECT` is the parent research-project directory. Source locations are explicit in `prepare_metadata.R` and the companion reconstruction inventory. Configure paths on a new machine and preserve historical locks separately from a new reproduction run.

`verify_Hallmark_reference.R` requires both normalization modes to reproduce the existing438-gene K562/RPE1 reference vectors and exact common libraries/target support before new readouts are released. `independent_numerical_QA.R` checks direct raw-count cell scores and exported GEM reaggregation. Execution logs and session information are retained. Inference belongs to the separately frozen `../specificity/` protocol; source-only alternate ranking choices are frozen in `../ranking/alternative_definitions/` before validation-table release.

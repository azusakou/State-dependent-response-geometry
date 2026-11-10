# Hallmark identifier-coverage sensitivity

This technical sensitivity adds four unambiguous HGNC identifier matches to the unchanged primary90-gene Hallmark readout:

| Assay symbol | Approved HGNC symbol |
|---|---|
| H2AFX | H2AX |
| IARS | IARS1 |
| TARS | TARS1 |
| TTC37 | SKIC3 |

These four were already determined by the frozen full113-member Hallmark snapshot and HGNC mapping before functional-readout outcomes were examined. Together with the original90, they give the94 measured canonical identifiers in `fixed_canonical_members.tsv`. Here, *canonical* refers to gene nomenclature, not canonical ER-UPR activity or mechanism.

`technical_addendum.txt` was locked at2026-11-06 06:02:16 UTC before94-response computation. Its manifest SHA256 is `7e78841ce2b5859b299343c63d2ce0b05fdad7613960b43cb85417fe0c4fd4e3`. `execution_lock_manifest.tsv` additionally fixes the code, narrow-count inputs and descriptive statistics specification before scoring. Earlier primary definitions, source locks, all80 functional/parent tests and all ranking lists remain unchanged.

None of the four genes was present in the real-readout raw-count caches. Partial control-target caches were insufficient to recover treated scores. `extract_four_alias_counts.R` therefore reads only these four count features in one pass per source and aligns them to the existing cached cells. Existing6,780-gene library totals and assigned-target counts are reused without modification. Source and RPE1 dense storage requires reading the relevant compressed row selections; Adamson CSR storage is streamed once, and Frangieh CSC storage permits direct gene-column reads. Where prior control-target counts overlap, their maximum absolute difference is zero. Input feature rows, file sizes, retained cells and elapsed times are recorded in `narrow_extraction_accounting.tsv`.

`score_canonical94.R` uses the existing control anchors, fixed gene cohorts and both normalization modes. `exact_Lminus_target` recomputes all log terms using the library excluding the target if it belongs to the6,780-gene denominator; `fixed_L6780` retains that denominator. Both remove the own targeted member from the numerator and apply identical operations to target-specific controls before calculating SD and matched-control means. The original90 readout is independently reconstructed by this same adapter before94 outputs are released.

`reference90_reconstruction_QA.tsv` verifies every retained cell and all original gene effects in both modes. Maximum gene-effect differences are below (10^{-13}), with identical target support. The primary90 outputs are not overwritten.

For each study, `*_gene_effects.tsv.gz`, `*_GEM_effects.tsv.gz` and `*_control_moments.tsv.gz` contain the94-readout sensitivity, with `axis=Hallmark_UPR94_canonical` and the same mode names and response units as the companion module analysis. `*_cell_scores.rds` contains compact real-cell scores and metadata. Its control rows are unassigned base scores; target-specific control references are supplied in the moment and GEM tables. Original raw90 counts remain in the companion caches; `*_additional_raw_counts.rds` stores only the four added features.

Reproduce from the parent research-project directory with R packages `data.table`, `Matrix`, `hdf5r` and `digest`:

```sh
Rscript paper19opp-revision2/upr_robustness/modules/canonical_Hallmark94/extract_four_alias_counts.R PROJECT
Rscript paper19opp-revision2/upr_robustness/modules/canonical_Hallmark94/score_canonical94.R PROJECT
```

Configure the documented raw-input paths for a relocated reproduction while preserving the original historical manifests. This is a post-primary technical robustness analysis with a definition locked before its own outcomes. It introduces no94-based ranking, new formal test family or new biological hypothesis. Only previously fixed cohorts and candidate lists are projected onto these responses. Existing90-member random modules are not a size-matched94-gene specificity benchmark.

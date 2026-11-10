# Revision2 analysis and manuscript

This snapshot supplies the manuscript, final figures,188 indexed Supplementary Data tables, compact frozen inputs, analysis code, statistical results and execution records. It preserves the original scientific choices. It excludes source correspondence and manuscript-editing histories.

## Read the manuscript

The current sources are `manuscript_revised_red.tex` and `supplementary_revised_red.tex`; their PDF review copies have17/57 pages. All required final figures and CAS template dependencies are supplied. Run from this directory:

```sh
tectonic -X compile manuscript_revised_red.tex
tectonic -X compile supplementary_revised_red.tex
```

References are embedded in each source. Existing red marks accumulate all changes relative to the final first revision. The curated-readout results and original main figure are now NoteS13/FigureS15; section3.5/Figure4 presents the independent HepG2 evaluation. `README_compile.txt` gives the conventional LaTeX alternative.

## Analysis map

| Directory | Analysis |
|---|---|
| `reconstruction/` | Common RNA measurements, study-held-out prediction, independent candidate/axis evaluation, source-frozen RPE1 tests and lineage validation. Begin with its README. |
| `upr_robustness/` | Frozen quintiles, ranking metrics/baselines, functional modules, alternative definitions, matched random modules and global multiplicity checks. |
| `independent_validation/` | Unused HepG2 screen: frozen readout, exact-composition scoring, literal lists, all baselines/nulls, uncertainty and support flags. Begin with its README. |
| `final_adjustment/` | Subsequent metadata-only Jurkat/X-Atlas coverage search. No new RNA outcomes were evaluated. |
| `supporting_analyses/` | Supporting sampling simulation and intervention-choice checks. |
| `supplementary_data/` | Final indexed numerical and figure-source tables, including every negative and unavailable result. |

RPE1 shares the Replogle study with the discovery experiment. HepG2 is from a separate study with overlapping investigators. Previously examined module/ranking comparisons are robustness analyses; they are not untouched independent validation. The primary HepG2 ordering endpoint uses417 genes; complete literal Top10/20 primary utility remains unavailable. The>=10-cell results stay secondary. These distinctions are not changed by publication of the code.

## Inputs, dependencies and execution

Run the historical reconstruction scripts from the repository root or with the explicit project/output arguments specified in their component READMEs. Run the independent HepG2 commands from this directory. The `paper19opp-revision1/GitHub` symlink points to the root legacy snapshot, whose scientific inputs/code match the source archive. On platforms that do not materialize Git symlinks, create an equivalent directory alias before running those historical paths.

The main R analysis uses `data.table`, `Matrix`, `hdf5r` and `digest`; particular analyses additionally use packages listed in their component READMEs and recorded sessionInfo files. Python source-download and metadata scripts use the standard library. No large model training is introduced in the revision2 UPR analysis.

Large original matrices and29 source-derived inputs/downloads over10MiB are excluded from this source repository; their byte sizes and SHA256 identities are in `LARGE_INPUTS.tsv`. The original reconstruction scripts and source registries are preserved, including their historical local raw-data paths. Supply those public data at corresponding paths or adapt input-path arguments without changing cohorts, gene sets, rankings, metrics or thresholds. The complete research archive accompanying the manuscript retains large processed intermediates. Downloads and recomputation can be substantial; a complete offline rerun is not claimed from this compact snapshot alone.

For independent HepG2 reproduction:

```sh
# Start here: cd paper19opp-revision2
python3 independent_validation/download_parallel.py
python3 independent_validation/run_when_ready.py
Rscript independent_validation/numerical_QA.R
Rscript independent_validation/create_figure.R
Rscript independent_validation/export_and_supplement.R
```

The first command downloads the checksum-pinned public matrix. `run_when_ready.py` performs extraction, frozen metrics and shared-GEM sensitivity; it requires the original matrix and frozen source inputs. The supplied compact gene effects and ranking inputs already support `analyze_holdout.R` for downstream metric reproduction. To preserve historical locks/results, reproduce in a copy. See `independent_validation/README.txt` for support and inference details.

Run `python3 tools/verify_revision2.py` from the repository root to check the compact snapshot. `GITHUB_SNAPSHOT.sha256` hashes every current revision file other than itself; existing source locks preserve their original role and are not rewritten to conceal omissions. Dataset licences, accession provenance and investigator attribution remain as recorded in source metadata and the root README.

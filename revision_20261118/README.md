# Revision analyses (2026-11-18)

Added after the four-team review. All scripts read-only; outputs listed.

| Script | Output | What it adds |
|---|---|---|
| `v11_replication_check.py` | `v11_suSie_check.tsv` | GTEx v11 fine-mapping (SuSiE) release-level check at 17q21.32: liver SKAP1 retains two credible sets containing the v8 lead (48,127,316) and the risk lead rs1860862 (48,127,474, PIP 0.015); muscle RP5-890E16.5 retains one set; SKAP1 has none in v11 pancreas, CBX1 none in v11 muscle |
| `dryrun_894_b2.py` | (module builder) | 100 expression-decile-matched random modules of 894 genes (seed 20261118), same procedure as the 862-gene matched null |
| `w3_incremental` re-run (894) | `incremental_result_894.tsv`, `incremental_null_894.tsv` | Bulk outcome layer with the corrected-annotation programme: dAUC +0.0001, p=0.14 |
| `w3_singlecell` re-run (894) | `singlecell_IR_test_894.tsv`, `singlecell_null_matched_894.tsv` | Single-cell layer: malignant-minus-normal -0.053, upper-tail p=1.00, lower-tail p<0.01 |
| `a4` re-run (894) | `correlation_894.tsv`, `null_normal_894.tsv` | Cross-cancer proliferation correlations -0.032 (tumour) / +0.472 (normal); adipocyte module +0.13 / +0.18 |

Manuscript: the unsupported "352 testable pairs (liver 0.922, muscle 0.874)" sentence was
replaced by the v11 fine-mapping statement above; no backing file for the 352 count existed.

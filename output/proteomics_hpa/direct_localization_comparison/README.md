# Direct comparison of HPA localization proportions: MSI2 324 versus 328

**No broad localization category differs significantly between 324 and 328 hit sets.** Across all three matched pairs, no nominal p-value is below 0.05 and no BH FDR is below 0.10. The lowest FDR is **0.4981**, for Vesicular/endolysosomal membership in the full-cluster pair.

## The mitochondrial example

For the original EGFP-filtered cluster sets, the exact table is:

| Original hit set | Mitochondrial annotation | No mitochondrial annotation |
|---|---:|---:|
| 324 | 36 | 201 |
| 328 | 37 | 203 |

The proportions are **36/237 = 15.19%** and **37/240 = 15.42%**, respectively. Two-sided Fisher p = **1.000**, BH FDR = **1.000**, conditional odds ratio (324 relative to 328) = **0.983**, exact 95% CI **0.578–1.670**. The difference is **-0.23 percentage points**. These data provide no evidence of differing mitochondrial-annotation proportions between these two hit sets; they do not establish equivalence.

## Inputs and method

This analysis reuses the exact original Reactome input sets and the fixed eight broad HPA categories from the [broad exploratory ORA](../broad_exploratory_ora/README.md). No hit selection, HPA mapping or confidence rule is changed. Multi-label membership is preserved, with one gene counted at most once within each category. Primary confidence means Enhanced, Supported or Approved.

For each matched pair and each category, compare category membership versus nonmembership between 324 and 328. Use **two-sided Fisher exact tests**, because the question is whether proportions differ in either direction. Odds ratio >1 means higher odds of category membership among 324 hits; <1 means higher odds among 328 hits. The confidence intervals are unadjusted exact 95% conditional odds-ratio intervals.

BH correction is across the **eight broad categories separately within each matched pair**, giving three families and 24 tests. A supplementary `BH_FDR_all_24` column pools all 24 tests. This is one between-set test per category, not separate one-sided 324 and 328 ORAs. The measured-proteome background defines the eligible original sets but does not enter these direct 2x2 tables; the rows contain the hit sets being compared.

### Shared-gene handling

Full-cluster and EGFP-filtered sets are disjoint and use their original counts directly. The original isoform-specific sets share **EIF4A2**. Treating this same gene as an independent observation in both rows is inappropriate for the usual Fisher sampling model.

Original-set membership and proportions are preserved in the `full_*` columns (172 and 300 genes). For that pair's statistical test, shared membership is held fixed and inference uses only the nonshared genes, **171 versus 299**. This is a comparison of the isoform-exclusive portions of the existing sets, not an independent-sample test of the two overlapping full sets. The `test_*` columns record the actual tables, and the audit identifies the shared gene and every original membership. No original gene-set file is altered. Plots show original-set descriptive proportions, while the specific-pair statistical results use the documented nonshared comparison.

## Results

| Matched pair | Original 324 / 328 sizes | Tested 324 / 328 sizes | Lowest nominal p | Lowest BH FDR |
|---|---:|---:|---:|---:|
| Full clusters | 421 / 609 | 421 / 609 | 0.0623 | 0.4981 |
| EGFP-filtered clusters | 237 / 240 | 237 / 240 | 0.1810 | 0.8444 |
| Isoform-specific subsets | 172 / 300 | 171 / 299 | 0.2352 | 0.8357 |

The smallest p-value comes from Vesicular/endolysosomal membership in the full clusters: **30/421 (7.13%)** for 324 versus **65/609 (10.67%)** for 328, odds ratio **0.642**, 95% CI **0.394–1.027**, nominal p **0.0623**, FDR **0.4981**. This is not significant.

### Nucleolus and Mitochondrial comparisons

| Pair | Category | Original 324 overlap | Original 328 overlap | Test odds ratio | Nominal p | BH FDR |
|---|---|---:|---:|---:|---:|---:|
| Full clusters | Nucleolus | 34/421 | 49/609 | 1.004 | 1.0000 | 1.0000 |
| Full clusters | Mitochondrial | 59/421 | 74/609 | 1.178 | 0.3959 | 0.7902 |
| EGFP-filtered | Nucleolus | 18/237 | 16/240 | 1.150 | 0.7249 | 0.8444 |
| EGFP-filtered | Mitochondrial | 36/237 | 37/240 | 0.983 | 1.0000 | 1.0000 |
| Isoform-specific | Nucleolus | 13/172 | 26/300 | 0.864 | 0.7312 | 0.8357 |
| Isoform-specific | Mitochondrial | 19/172 | 40/300 | 0.810 | 0.5632 | 0.8357 |

For the specific pair, displayed overlaps use the original descriptive denominators, while odds ratios/p-values use the nonshared test denominators described above. Full exact test counts and confidence intervals are provided in the tables.

These direct comparisons do not support a broad localization difference between the isoform-associated hit sets. They also show why nominal mitochondrial ORA enrichment relative to the measured proteome need not imply different mitochondrial proportions between 324 and 328: the EGFP-filtered sets have nearly identical proportions. This is not evidence of identical localization or a formal equivalence test, and does not resolve finer functional pathway differences.

## Deliverables and validation

- `direct_comparison_all.tsv`: all 24 comparisons, original and test 2x2 counts, proportions, percentage-point differences, odds ratios, 95% confidence intervals, p-values, both FDR scopes, and contributing genes.
- `direct_cluster.tsv`, `direct_hits.tsv`, `direct_specific.tsv`: individual pair tables with per-pair BH correction.
- `nucleolus_mitochondria_comparisons.tsv`: the two biological hypotheses.
- `original_membership_and_shared_gene_audit.tsv` and `contributing_genes_long.tsv`: original memberships, shared-gene treatment and category contributors.
- `analysis_summary.tsv`, `input_checksums.tsv`, `sessionInfo.txt` and `validation_checks.tsv`: summary and reproducibility records. All checks passed, including an independent two-sided conditional-hypergeometric reconstruction of Fisher p-values and verification of unchanged inputs.
- `figures/`: three paired PDF/PNG plots of original-set proportions. `significant_FDR005.tsv` is header-only because no comparison passes that threshold.

Reproduce from the repository root:

```sh
Rscript analysis/proteomics_hpa/run_direct_localization_comparison.R
```

HPA nonmembership means no retained annotation to the category, not proof of biological absence. Annotations are gene-level priors from other cell contexts. These exploratory gene-membership tests do not model sample-level proteomics uncertainty or dependence among interacting proteins. All earlier mapping, historical-input reconstruction and exploratory-selection limitations remain applicable.

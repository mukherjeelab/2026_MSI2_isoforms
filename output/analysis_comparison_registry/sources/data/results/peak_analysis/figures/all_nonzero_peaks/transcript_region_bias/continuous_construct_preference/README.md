# Continuous construct-preference analysis

This directory contains the primary continuous analysis of whether relative
324WT-versus-328WT CLIP enrichment differs by transcript-region annotation.
It is separate from the older threshold-based delta-group analysis in the
parent directory.

The `quantile_normalized_log2cpm/` subdirectory contains a separate
normalization sensitivity analysis. It converts raw counts to full-library
CPM, applies `log2(CPM + 0.5)`, quantile-normalizes the eight sample columns,
and then recalculates construct enrichment and delta-delta. The primary results
described below are unchanged.

## Peak universe and score

The analysis includes all 16,046 merged peaks with a raw count greater than
zero in at least one of the eight libraries. This filter is independent of
delta-delta sign, preferred construct, and transcript-region annotation.
No additional abundance or direction-dependent filter is applied.

For each peak:

`delta_delta = 324WT replicate-mean log2(CLIP/SMI) - 328WT replicate-mean log2(CLIP/SMI)`

Each construct-specific enrichment is the existing CPM-normalized ratio
`log2((mean CLIP CPM + 0.5)/(mean SMI CPM + 0.5))`. Positive delta-delta
favors 324WT and negative delta-delta favors 328WT.

## Statistical tests

The overall five-region Kruskal-Wallis test has H = 132.15, df = 4, p = 1.348e-27, and epsilon-squared = 0.007989.
Pairwise Wilcoxon rank-sum tests compare locations and pairwise two-sample
Kolmogorov-Smirnov tests compare complete distributions. Each family of ten
pairwise p-values is adjusted separately using Benjamini-Hochberg.

For 3' UTR versus intron, the median difference (3' UTR minus intron) is -0.7322, Cliff's delta is -0.1592 (small), Wilcoxon BH q = 3.868e-15, and KS BH q = 6.254e-34.

Small p-values should be interpreted together with median differences and
rank-based effect sizes because the peak counts are large. These are
peak-level exploratory tests; peaks from the same gene or transcript may not
be statistically independent.

## Files

- `continuous_delta_delta_by_region.{pdf,png}`: violin/boxplot distributions
  for all five transcript-region categories, with n and median in labels.
- `construct_enrichment_scatter_by_region.{pdf,png}`: density-binned 324WT
  versus 328WT enrichment facets with common axes and the y = x diagonal.
- `delta_delta_3UTR_vs_intron_density.{pdf,png}`: focused density comparison.
- `delta_delta_3UTR_vs_intron_ecdf.{pdf,png}`: focused ECDF comparison.
- `continuous_delta_delta_peaks.tsv`: all peak-level values used in figures
  and tests.
- `delta_delta_region_summary.tsv`: n, median, quartiles, mean, SD, and sign
  percentages by transcript region.
- `delta_delta_region_overall_tests.tsv`: overall Kruskal-Wallis result and
  epsilon-squared effect size.
- `delta_delta_region_pairwise_tests.tsv`: pairwise Wilcoxon and KS tests, BH
  q-values, median differences, and Cliff's delta.
- `sessionInfo.txt`: R and package versions used for this run.

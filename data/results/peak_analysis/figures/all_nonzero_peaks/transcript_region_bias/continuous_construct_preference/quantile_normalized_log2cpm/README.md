# Quantile-normalized log2-CPM sensitivity analysis

This directory contains a normalization sensitivity analysis of the same
16,046 construct-independent counted peaks used by the primary continuous
construct-preference workflow. It does not replace the original CPM-based
analysis.

## Normalization

1. Subset the raw 43,794-peak matrix to peaks with a count above zero in at
   least one of the eight libraries (16,046 peaks).
2. Divide each raw count by its full mapped-library size and multiply by one
   million to obtain CPM.
3. Transform each sample as `log2(CPM + 0.5)`.
4. Quantile-normalize the eight transformed sample columns together using
   `preprocessCore::normalize.quantiles`.
5. For each construct, subtract the mean normalized SMI abundance from the
   mean normalized CLIP abundance. Delta-delta is 324WT enrichment minus
   328WT enrichment.

This filter is independent of delta-delta sign, preferred construct, and
transcript-region annotation. Quantile normalization maps sample ranks to a
shared reference distribution. Tied zero-count peaks remain tied and receive
their averaged reference rank, so sparse libraries do not have exactly
identical empirical marginals after normalization. The procedure is still
deliberately aggressive and should be interpreted as a sensitivity analysis
rather than automatically as a superior normalization.

## Results

The overall Kruskal-Wallis test has H = 62.881, df = 4, p = 7.189e-13, and epsilon-squared = 0.003671.
For 3' UTR versus intron, the median difference is -0.2228, Cliff's delta is -0.1223 (negligible), Wilcoxon BH q = 3.809e-09, and KS BH q = 6.461e-32.
Compared with the primary CPM-pseudocount delta-delta values, the quantile-normalized values have Pearson r = 0.8924, Spearman rho = 0.8841, and 92.85% sign agreement.

Because sample sizes are large, p-values should be interpreted alongside
median differences and Cliff's delta. These are peak-level exploratory tests
and do not model dependence among peaks from the same gene or transcript.

## Files

- `quantile_normalized_log2_cpm.tsv`: the normalized 16,046 by 8 matrix.
- `normalization_sample_summary.tsv`: sample distributions before and after
  quantile normalization.
- `normalization_comparison_to_primary.tsv`: correlation, sign agreement,
  and overall medians versus the primary CPM-pseudocount delta-delta values.
- `continuous_delta_delta_peaks.tsv`: normalized construct enrichments,
  delta-delta values, and transcript annotations for every eligible peak.
- `continuous_delta_delta_by_region.{pdf,png}`: five-region continuous plot.
- `delta_delta_3UTR_vs_intron_density.{pdf,png}` and
  `delta_delta_3UTR_vs_intron_ecdf.{pdf,png}`: focused comparisons.
- `delta_delta_region_summary.tsv`, `delta_delta_region_overall_tests.tsv`,
  and `delta_delta_region_pairwise_tests.tsv`: summaries and exploratory tests.
- `sessionInfo.txt`: R and package versions.

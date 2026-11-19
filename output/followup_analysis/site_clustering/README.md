# Reference MSI2 within-feature site clustering

This is a new distance analysis of the existing reference summit assignments. Introns and 3′UTRs are analyzed separately. No construct-specific site selection is used.

## Primary numerical results

| Feature type | Features with >=2 sites | Summits | Observed median NN (nt) | Null median (95% interval; nt) | Observed/null |
|---|---:|---:|---:|---:|---:|
| intron | 3,799 | 14,598 | 644 | 2,892 (2,795.1–3,000.1) | 0.223 |
| 3UTR | 5,317 | 22,952 | 252 | 241 (236.0–247.0) | 1.046 |

| Feature type | Neighbor threshold | Observed fraction | Null mean (95% interval) | Enrichment | One-sided clustering P |
|---|---:|---:|---:|---:|---:|
| intron | <=50 nt | 2.12% | 2.81% (2.44%–3.17%) | 0.756 | 1.000000 |
| intron | <=100 nt | 6.90% | 5.32% (4.85%–5.79%) | 1.297 | 0.001996 |
| intron | <=250 nt | 28.80% | 11.62% (11.07%–12.23%) | 2.479 | 0.001996 |
| intron | <=500 nt | 45.64% | 19.38% (18.69%–20.07%) | 2.355 | 0.001996 |
| 3UTR | <=50 nt | 5.10% | 14.83% (14.28%–15.39%) | 0.344 | 1.000000 |
| 3UTR | <=100 nt | 12.96% | 26.63% (25.97%–27.25%) | 0.487 | 1.000000 |
| 3UTR | <=250 nt | 49.66% | 51.18% (50.44%–51.86%) | 0.970 | 1.000000 |
| 3UTR | <=500 nt | 76.91% | 72.71% (72.15%–73.21%) | 1.058 | 0.001996 |

At 250 nt, intron enrichment is 2.37x, 2.46x, and 2.51x in the 2-site, 3–4-site, and >=5-site strata, respectively. Intron median NN distance is substantially below its null. The 3′UTR profile differs: it is depleted at 50 and 100 nt, close to the expectation at 250 nt, and modestly enriched at 500 nt. Both feature types are depleted at 50 nt. Thus the intron excess is scale-dependent, not evidence of excess at every short-distance threshold.

## Eligibility and assignment

- Input: `../tables/summit_normalized_positions.tsv.gz`, created by `../scripts/site_spatial.R`. Source feature metadata: `data/results_round2/features.tsv.gz` (paths relative to the project root where applicable).
- Original longest matching feature assignments are retained, including strand-aware transcript-oriented offsets. No peak is reassigned. A reference summit may independently belong to each feature type in different transcript models, as in the original analysis; feature types are never pooled.
- Introns: length strictly >500 nt. 3′UTRs: original minimum length >=100 nt. Retain only exact assigned feature IDs containing >=2 summits after the length filter.
- Distances use transcript coordinates. For spliced 3′UTRs, gaps across genomic introns are excluded from the coordinate distance; intron distances use the assigned intron interval.
- Primary weighting: one nearest-neighbor distance per summit. Features containing more summits contribute more distances; the exact same site counts and weighting apply to every null draw. Site-count strata are 2, 3–4, and >=5.

## Distance definitions

For sorted within-feature positions x1,...,xn, consecutive gaps are x[i+1]-x[i]. The nearest-neighbor distance for an internal summit is the smaller adjacent gap; the two endpoint summits use their sole adjacent gap. With two summits, both nearest-neighbor entries equal the single consecutive gap. For >=3 summits these distributions are generally different, not merely edge-case equivalents. No cross-feature distances are computed. Threshold fractions include equality (distance <= threshold). Equal-position ties, if present, have distance zero and are retained.

## Null and inference

- 500 iterations; R seed 26092201. The script calls set.seed once and processes introns then 3UTRs in a fixed, sorted row order. RNG/software details are saved in session_info.txt.
- Independently randomize each assigned site on the half-integer grid 12.5:(L-12.5), implemented as 12.5 + floor(runif(n)*(floor(L-25)+1)). This is the same end-margin and sampling convention as the current spatial-background adaptation.
- Exact feature identity, length, strand, and site count are preserved. Positions are sampled in each feature’s transcript orientation. Sampling allows replacement/coincident sites; no pairwise exclusion or extra spacing rule is imposed. This preserves the existing adaptation, rather than adding the historical PARPIPE inter-site spacing logic.
- Observed positions are retained unchanged, including 15 intron and 113 3UTR sites outside the null end margins (see margin diagnostics).
- Each iteration recalculates all within-feature nearest-neighbor and consecutive distances, then the all-feature and site-count-stratified statistics.
- For the median-distance statistic, null_median is the median of the 500 iteration-level medians; the 95% null interval is their 2.5th–97.5th percentiles. For threshold fractions, the null expectation is the mean of the 500 fractions; intervals are their 2.5th–97.5th percentiles. These are null intervals, not confidence intervals on observed measurements.
- observed_null_ratio divides observed median distance by the median of null medians, or observed threshold fraction by the mean null fraction. A distance ratio <1 or threshold ratio >1 indicates clustering relative to this null.
- One-sided clustering P: (1 + count(null statistic <= observed median distance))/501 for medians; (1 + count(null fraction >= observed fraction))/501 for thresholds. Ties count as equally extreme. Minimum possible P is 1/501 = 0.001996. Both tail P values and a conservative doubled-smaller-tail P are also supplied; P=1 for the clustering direction can accompany a significant depletion.
- q_BH_clustering adjusts the 40 tests (2 feature types x 4 strata including All x 5 metrics) together, separately for nearest-neighbor and consecutive-gap outcomes. The feature types remain separate outcomes throughout.

## Interpretation limits

The immediate inference is clustering relative to uniformly randomized positions conditioned on exact feature lengths and site counts. It is not evidence of cooperative binding. The null does not condition on the previously observed donor-proximal positional preference, local sequence, accessibility, mappability, or the spacing/resolution imposed during peak calling. Any of these can influence the observed distance distribution. In particular, short-distance depletion is reported as a departure from this null, without assigning a biological mechanism.

## Visualization

- Red: observed nearest-neighbor density. Transparent cyan: all 500 independently calculated null density curves. Black: pointwise median null density. Separate intron and 3UTR panels, sharing x limits but with separate y scales clearly labeled.
- Display-only reflected Gaussian KDE at the nonnegative boundary, fixed bandwidth 25 nt for all observed/null curves and both feature types. All tests and summaries use unsmoothed distances. KDE source grids are saved so the plots can be regenerated without drawing new random positions.
- The broad view uses a shared upper limit covering at least the 99th percentile of both observed and pooled-null distances within each feature type (rounded up to 1,000 nt). It is not the absolute maximum. A separate 0–1,000 nt zoom reveals local structure. Neither view renormalizes after truncation; exact excluded tail fractions and bandwidth are saved in MSI2_density_display_parameters.tsv.
- Threshold plot: observed fraction in red, null mean in dark blue, pointwise 95% null interval shaded blue. Columns show All, 2, 3–4, and >=5-site strata. Numerical thresholds are the primary results.

## Reproduction

From the project root:

```sh
Rscript output/followup_analysis/site_clustering/analyze_site_clustering.R
Rscript output/followup_analysis/site_clustering/plot_site_clustering.R
```

Running only the plotting script reuses saved distance matrices and never regenerates randomizations. Requirements: R, data.table, ggplot2 and scales. The analysis script checks source hashes and writes numerical summaries before plotting is run.

## Saved files

- `MSI2_site_clustering_summary.tsv`: all-feature primary metrics, feature/summit counts, null intervals, ratios and P values.
- `MSI2_site_clustering_by_sitecount.tsv`: primary metrics stratified by site count.
- `MSI2_site_clustering_thresholds.tsv`: all threshold results including strata.
- `MSI2_site_clustering_null_iterations.tsv.gz`: all 500 iteration statistics for both distance definitions and each stratum.
- `MSI2_consecutive_site_distance_summary.tsv`: supplementary consecutive-gap statistics. Fractions here use gaps as their denominator, not summits.
- `MSI2_site_clustering_all_statistics.tsv`: complete combined statistics.
- `MSI2_observed_nearest_neighbor_distances.tsv.gz`, `MSI2_observed_consecutive_distances.tsv.gz`: individual observed distances and feature/peak identifiers.
- `MSI2_site_clustering_included_summits.tsv.gz`, `MSI2_site_clustering_features.tsv.gz`: exact included assignments, feature metadata and transcript segments.
- `MSI2_site_clustering_eligibility.tsv`, `MSI2_site_clustering_margin_diagnostics.tsv`: eligibility counts and end-margin diagnostics.
- `cache/intron_site_distances.rds`, `cache/3UTR_site_distances.rds`: complete randomized position, nearest-neighbor, and consecutive-gap matrices; columns are iteration 1–500. Position and NN rows map to site_map; consecutive-gap rows map to gap_feature_id and ascending within-feature gap rank. Randomized peak labels are exchangeable bookkeeping, not biological identities.
- `MSI2_site_distance_density_inputs_broad.tsv.gz`, `MSI2_site_distance_density_inputs_zoom.tsv.gz`, `MSI2_site_clustering_threshold_plot_inputs.tsv`: exact plotting inputs.
- Requested PNG/PDF plots, source scripts, input hashes, final RNG state and session information are saved alongside these outputs.

## Figure caption

Reference MSI2 nearest-neighbor summit distances within introns >500 nt and spliced 3′UTRs >=100 nt containing at least two assigned summits. Red curves show observed distances; cyan curves show all 500 feature-matched randomizations, preserving exact feature lengths, strands and site counts with the existing 12.5 nt end-margin convention; black curves show the pointwise median null density. Densities use a reflected 25 nt Gaussian bandwidth for display only. Numerical inference uses raw summit distances. Intronic sites show an excess of neighbors at 100–500 nt relative to this null, with the strongest reported enrichments at 250 and 500 nt; the 3′UTR pattern differs. These comparisons assess local clustering and do not establish cooperative binding.

## Independent validation

`validate_site_clustering.R` independently reconstructs nearest-neighbor distances using full pairwise distance matrices within every feature, for observed positions and null iterations 1, 250 and 500. It also reconstructs consecutive gaps and checks all 500 iterations for valid positions, the half-integer grid and exact site counts. All checks passed; see `validation.tsv`. Original source MD5 hashes remain unchanged. `output_manifest.tsv` records SHA-256 hashes of the generated files.

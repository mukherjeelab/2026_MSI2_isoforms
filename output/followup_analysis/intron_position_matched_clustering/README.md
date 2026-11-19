# Intron clustering after matching positional preference

Two new reference-MSI2 intron tests extend the existing nearest-neighbor analysis. All introns remain strictly >500 nt, original feature assignments and strand-aware coordinates are retained, and no 3′UTR analysis is changed.

## Main comparison

| Null model / cohort | Introns | Summits | Observed median NN (nt) | Null median (95% interval; nt) | Observed <=250 nt | Null mean <=250 nt (95% interval) | 250 nt enrichment | Clustering P at 250 nt |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Original uniform, full introns | 3,799 | 14,598 | 644 | 2,892 (2,795.1–3,000.1) | 28.80% | 11.62% (11.07%–12.23%) | 2.479x | 0.001996 |
| Position-matched, full introns | 3,799 | 14,598 | 644 | 2,829 (2,728.2–2,929.3) | 28.80% | 12.09% (11.50%–12.64%) | 2.382x | 0.001996 |
| Uniform within first 20% | 1,089 | 3,475 | 410 | 935 (877.5–1,002.5) | 34.82% | 22.75% (21.21%–24.32%) | 1.530x | 0.001996 |

The 250 nt excess persists after matching the genome-wide relative-position preference: enrichment changes from 2.48x to 2.38x. It also persists when restricting both observed and randomized sites to the first 20% of an intron (1.53x). Observed nearest-neighbor medians remain shorter than expected in both new tests. These results support clustering beyond the pooled, 1%-binned positional preference represented by the primary null. They do not establish cooperative binding or rule out other spatial covariates.

## Site-count-stratified 250 nt enrichment

| New null | Site-count stratum | Introns | Summits | Enrichment | Clustering P |
|---|---|---:|---:|---:|---:|
| Position-matched, full introns | 2 sites | 1,586 | 3,172 | 2.286x | 0.001996 |
| Position-matched, full introns | 3-4 sites | 1,285 | 4,346 | 2.378x | 0.001996 |
| Position-matched, full introns | >=5 sites | 928 | 7,080 | 2.401x | 0.001996 |
| Uniform within first 20% | 2 sites | 556 | 1,112 | 1.458x | 0.001996 |
| Uniform within first 20% | 3-4 sites | 359 | 1,185 | 1.424x | 0.001996 |
| Uniform within first 20% | >=5 sites | 174 | 1,178 | 1.643x | 0.001996 |

At 100 nt, the position-matched full-intron cohort has a 1.22x excess, whereas the first-20% cohort has a 0.81x ratio (depletion rather than enrichment). The conclusion is distance-scale dependent; it should not be described as excess at every short-distance threshold. Full 50/100/250/500 nt results are provided in the summary and stratum tables.

## Primary position-matched null

1. Fit a site-weighted empirical histogram of relative intron positions in the same 100 bins used in the existing positional analysis. Its source is every reference summit assigned to an eligible intron >500 nt: 20,364 summits in 9,565 introns, including single-site introns. This uses the genome-wide eligible positional preference, not a separate distribution estimated for each intron.
2. Preserve the original clustering cohort: 3,799 exact introns and 14,598 summits in introns with >=2 assigned sites. Preserve feature ID, length, strand and site count. Positions are already in transcript orientation, with 0% at the donor-side intron beginning.
3. For each site, sample a relative-position bin with probability equal to its empirical site fraction, then sample uniformly within that bin. Convert u to the recipient intron’s base-centered coordinate floor(u*L)+0.5.
4. Reject and redraw positions outside [12.5, L-12.5]. This conditions the empirical position density on the same half-integer grid and end-margin support as the existing PARPIPE-style adaptation. It does not clip positions onto boundaries. Empirical preferences are retained in expectation, subject to this per-intron conditioning, not as exact fixed bin counts in every realization.
5. Independently generate all sites with replacement, as in the prior null, allowing coincident positions. Do not add minimum pairwise spacing. Recalculate raw within-intron nearest-neighbor distances for each of 500 iterations. Fixed R seed: 26092301.

## Secondary first-20% null

- Select original assigned summits whose center-based relative coordinate is strictly <0.20. Retain only introns with >=2 summits remaining within the zone: 1,089 introns / 3,475 summits.
- Preserve each exact intron, length, strand and observed number of sites in its first-20% zone. Analyze nearest neighbors only among the retained zone sites. Sites outside the zone are not candidate neighbors.
- Draw independently and uniformly on half-integer coordinates beginning at 12.5 nt and strictly below 0.2*L. The formula is 12.5 + floor(runif(n)*ceiling(0.2*L-12.5)). There is no extra 12.5 nt exclusion at the artificial 20% boundary. The true intron end margin is automatically satisfied.
- Run 500 iterations, fixed R seed 26092302. Distances remain in original nt coordinates, not rescaled zone units. Strata use the number of sites inside the zone, not total sites across the full intron.
- This is a different cohort from the full-intron comparisons; its enrichment ratio should not be used to estimate the percentage of the original signal explained by localization. Uniform sampling inside the zone also does not account for a finer positional gradient within that zone.

## Statistics and interpretation

- Primary metric: one nearest-neighbor distance per summit within the same feature. Internal summits use the smaller adjacent gap; endpoints use their single adjacent gap. Site-weighting is identical between observed and null for each cohort.
- Medians: compare observed median to the median of 500 iteration medians. Thresholds: compare the observed fraction with neighbor distance <=50, <=100, <=250, or <=500 nt to the mean fraction across iterations. The 95% null intervals are the 2.5th/97.5th percentiles of the iteration-level statistic.
- Empirical one-sided clustering P=(1+number of null statistics at least as clustering-extreme as observed)/501. Smaller medians and larger threshold fractions are the clustering directions. Ties are included. The smallest resolvable P is 1/501=0.001996; this is not a claim that the true tail probability equals exactly 0.002.
- Both tail P values and a conservative doubled-smaller-tail P are saved. BH adjustment includes all 40 new tests (2 models x 4 strata including All x 5 metrics). Original uniform results retain their original adjustment family and are included as a historical comparison.
- The position preference is fitted to the observed data and treated as fixed. The P values are conditional Monte Carlo results under this plug-in null; uncertainty in estimating the positional density is not propagated. The null reproduces the pooled 1%-binned shape, not site-count-specific, length-specific or intron-specific positional densities.
- Sequence motifs, mappability, local accessibility, peak-calling resolution and finer-scale intron architecture remain unmodeled. Persistence supports additional clustering beyond the modeled global donor-proximal preference, not a mechanism or cooperative binding.
- Original observed sites were not altered to enforce the null margins: 15 full-intron sites and 4 first-20% sites lie outside the randomization support; see cohort_diagnostics.tsv.

## Validation

`validate_position_matched.R` independently uses pairwise distance matrices for every feature to reconstruct observed nearest-neighbor distances and null iterations 1, 250 and 500. It checks all 500 randomizations for half-integer coordinates, allowed support and exact per-feature counts. Full-intron IDs, peak IDs, observed positions and observed distances match the previous uniform analysis exactly.

The validator also integrates the fitted histogram analytically over each recipient intron’s allowed grid to calculate the expected 1%-bin probabilities after margin conditioning. Mean simulated probabilities agree with those expectations within Monte Carlo error (maximum absolute standardized deviations 2.30 for the position-matched null and 1.94 for the first-20% null; validation threshold 6). Input file hashes are unchanged.

## Figures

- `intron_position_matched_thresholds.png/pdf`: original uniform, position-matched and first-20% threshold curves, each stratified by site count. Shading is the pointwise 95% null interval.
- `intron_position_matched_density_zoom.png/pdf`: nearest-neighbor density over 0–1,000 nt, all 500 cached null curves in cyan, null median in black, observed in red. Shared axes; the third panel uses the zone-restricted cohort. Display-only reflected Gaussian KDE with fixed 25 nt bandwidth; no numerical test uses smoothed distances. Densities are not renormalized to the displayed range.
- `intron_position_preference_check.png/pdf`: global empirical positional distribution, simulated mean and analytic expected profile; this checks whether the primary null actually follows the fitted preference.

## Reproduction

Run from the project root:

```sh
Rscript output/followup_analysis/intron_position_matched_clustering/analyze_position_matched.R
Rscript output/followup_analysis/intron_position_matched_clustering/validate_position_matched.R
Rscript output/followup_analysis/intron_position_matched_clustering/plot_position_matched.R
```

The plotting script reads existing caches and never regenerates nulls. R dependencies: data.table, ggplot2 and scales. Session information and input hashes are saved.

## Tables and caches

- `intron_position_matched_summary.tsv`: all-feature results for both new models.
- `intron_position_matched_by_sitecount.tsv`: site-count-stratified results.
- `intron_position_matched_all_statistics.tsv`: all new tests, both tails, ratios and adjusted P values.
- `uniform_vs_position_matched_vs_first20.tsv`: new statistics alongside the original intron-uniform results.
- `intron_position_matched_null_iterations.tsv.gz`: every iteration-level statistic for each model and stratum.
- `empirical_position_source_summits.tsv.gz`, `empirical_intron_position_distribution.tsv`: exact empirical source and bin probabilities.
- `full_intron_analysis_summits.tsv.gz`, `first20_analysis_summits.tsv.gz`: exact analyzed site assignments and metadata.
- `first20_zone_feature_counts_before_minimum.tsv`: zone counts before requiring >=2 zone sites.
- `cache/position_matched_distances.rds`, `cache/first20_uniform_distances.rds`: all randomized positions, all nearest-neighbor distances, all positional profiles, observed positions/distances, row-to-feature/peak maps, model probability vectors and final RNG states. Matrix columns correspond to iterations 1–500.
- `positional_null_diagnostics.tsv`, `cohort_diagnostics.tsv`, `position_matching_validation.tsv`, `validation.tsv`: cohort and positional checks.
- `density_zoom_plot_inputs.tsv.gz`, `threshold_plot_inputs.tsv`: exact plotted values.
- `input_provenance.tsv`, `session_info.txt`, `output_manifest.tsv`: provenance and file hashes.

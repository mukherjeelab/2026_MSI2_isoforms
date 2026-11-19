# Intron summit combined null-cloud and percentile restyle

## Figure caption

Reference MSI2 summit positions within introns longer than 500 nt, shown in transcript-oriented normalized 0-100% intron coordinates. Red: observed MSI2 profile. Transparent blue curves: all 500 cached randomized profiles from the existing PARPIPE-style background procedure. Light blue band: 10th-90th percentiles; darker blue: 25th-75th percentiles. Black line: null median. No gray outer envelope is drawn. Left, each represented intron contributes total weight one; right, each summit contributes equally. Profiles use unsmoothed 1% bins. The observed excess toward the intron beginning is consistent with donor-proximal enrichment across multiple relative-position bins; these bins do not specify a fixed nucleotide distance from the splice site. This is a restyle of the existing analysis, with no new randomizations.

## Reproduction and provenance

Run from the project root: `Rscript output/followup_analysis/figure_restyles/plot_intron_null_cloud_percentiles.R`.

- Observed profiles: `output/followup_analysis/intron_gt500/intron_gt500_profiles.tsv`.
- Cached null array: `output/followup_analysis/intron_gt500/PARpipe_intron_gt500_nulls.rds`.
- Original generating code: `output/followup_analysis/scripts/intron_gt500.R`.
- All 500 draws are displayed in each panel, with no subset selection or smoothing.
- Existing tables, caches, analytical figures, assignments, >500 nt filter, weighting, and numerical results are unchanged.
- Blue overlap shows the cached null cloud; it is not a p-value color scale.
- The null remains the prior independent half-integer positional-grid adaptation with 12.5 nt end margins, not a new implementation of the full historical PARPIPE script.
- See plot_provenance_cloud_percentiles.tsv for source MD5 hashes, array-slice mapping, seed, software versions, and validation.

## Checks

Verified 100 bins and 500 draws per weighting scheme; cached medians and 10/25/75/90 percentiles match the saved profile table; every observed/null profile sums to one; source MD5 hashes are identical before and after plotting.

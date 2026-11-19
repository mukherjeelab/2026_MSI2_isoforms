# CLIP-only 324 minus 328 delta by transcript region

Four descriptive violin/boxplots compare two fixed peak sets with and without quantile normalization. Library-depth CPM scaling is retained in both, as requested. SMI values are not used for normalization or delta calculation. Existing source files and earlier analyses are unchanged.

## Peak sets

- **16,046 peaks:** original construct-independent universe with at least one raw count across the eight original libraries. Region membership matches the supplied delta-delta figure exactly.
- **4,385 peaks:** existing QC set in `data/results/peak_analysis/figures/QC/VST_counts_for_QC.tsv`; verified to match total raw count >=10 across the eight original libraries exactly. Only its peak IDs are reused, not its VST values. This set contains 3,913 3′UTR, 221 CDS, 197 Other, 42 Intron and 12 5′UTR peaks.
- Thus SMI contributes to the historical membership definitions, but is excluded from the new abundance calculation and quantile normalization. No new CLIP-positive filter or delta cutoff was applied.
- The nearby 4,315-row omniCLIP score-validation file is not a 4k unique-peak set: it contains 3,220 unique peak IDs, repeated across construct comparisons.

## Calculation

For each of the four CLIP libraries:

`CPM = raw peak count / full-library primary mapped alignments * 1,000,000`

`x = log2(CPM + 0.5)`

Without quantile normalization:

`delta = mean(x for 324 CLIP replicates) - mean(x for 328 CLIP replicates)`

With quantile normalization:

1. Apply `preprocessCore::normalize.quantiles` to the peak-by-four-CLIP-library matrix of x values, with its standard handling of ties.
2. Compute the same difference of replicate means on that normalized matrix.

There are two biological replicates per construct. Positive delta indicates higher 324 CLIP abundance; negative indicates higher 328 abundance. This is a log-scale abundance contrast, not a difference of linear CPM values and not CLIP/SMI enrichment. Replicates are averaged **after** log transformation in both branches; the calculation is not log2 of the arithmetic mean CPM. This consistent order isolates the effect of adding quantile normalization across the two branches.

Quantile normalization is fitted separately within each peak set, across its four CLIP columns. Consequently, changing from 16k to 4k changes both membership and the quantile-normalization target; the 4k normalized values are not simply a subset of the 16k normalized values. No eight-library normalization matrix or SMI-derived factors are reused. Full mapped-library denominators are not recomputed from peak subsets.

## Interpretation and sparse-count diagnostics

These are descriptive sensitivity plots, with no new significance tests or binding calls. Quantile normalization imposes a shared rank-based abundance scale; it can remove global abundance shifts and is not proof that the resulting direction is biologically more accurate.

The 16,046-peak set includes 6,212 peaks with zero counts in all four CLIP libraries. Their CPM-only deltas are zero, explaining large zero-centered densities and several zero regional medians. The 4,385-peak set includes 174 such peaks. Standard tied-rank quantile normalization in that smaller set gives these 174 rows nonzero deltas despite their absence of CLIP counts; these values should not be interpreted as measured binding differences. All these rows are flagged in the peak-level output. No rows were silently removed or reset to zero after normalization.

The 4k set has only 42 intron peaks and 12 5′UTR peaks. The corresponding violin shapes should be read with those small counts in mind.

## Figures and files

Each figure is provided as PNG and PDF, with the same x-axis range, region colors and order as each other. Labels show peak number and median; white boxes show median/IQR, black dots mark medians and the red reference line is zero. Violin widths are scaled equally by region and do not encode peak count.

- `CLIP_delta_peaks_16046_CPM_only.png/pdf`
- `CLIP_delta_peaks_16046_quantile_normalized.png/pdf`
- `CLIP_delta_peaks_4385_CPM_only.png/pdf`
- `CLIP_delta_peaks_4385_quantile_normalized.png/pdf`
- `*_peak_deltas.tsv.gz`: exact per-peak annotations, construct means, deltas and zero-CLIP flags used for plotting.
- `*_raw_CLIP_counts.tsv.gz`, `*_CLIP_CPM.tsv.gz`, `*_log2_CLIP_CPM.tsv.gz`, `*_QN_log2_CLIP_CPM.tsv.gz`: all calculation matrices, separately for each dataset.
- `CLIP_delta_region_summary.tsv`: region counts, medians, quartiles, means and zero-CLIP counts.
- `CLIP_delta_diagnostics.tsv`, `normalization_sample_summary.tsv`: dataset and library summaries.
- `CLIP_library_metadata.tsv`: exact four libraries, constructs, replicates and CPM denominators.
- `input_provenance.tsv`, `session_info.txt`, `validation.tsv`, `output_manifest.tsv`: source hashes, software versions, independent checks and output hashes.

## Reproduce

From the project root:

```sh
Rscript output/followup_analysis/clip_only_delta/plot_clip_only_delta.R
```

Requirements: R, data.table, ggplot2 and preprocessCore. The script verifies membership, unique annotations, four CLIP-only sample columns, tied-rank behavior, a known quantile-normalization example, and replicate-pair equivalence of delta. Independent Python checks reconstructed CPM and every saved delta from the matrices. Input MD5 hashes were unchanged after execution. All four PNGs were visually inspected.

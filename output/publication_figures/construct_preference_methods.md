# Construct preference by transcript region

Reproduce from the project root with `Rscript output/publication_figures/plot_construct_preference.R`.

The input is the original 9,369 retained peaks in `data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/relaxed_filter_continuous_delta_delta_peaks.tsv`. The original selection rule and per-peak values are unchanged. The plotted value is MSI2-324 replicate-mean log2(CLIP/SMI) minus MSI2-328 replicate-mean log2(CLIP/SMI), with the original 0.5 CPM pseudocount. Peaks depleted in both constructs were excluded by the source analysis.

5' UTR is merged into Other. Counts: CDS 398; 3' UTR 8,271; Intron 290; Other 410. Other comprises 28 5' UTR, 240 Ambiguous UTR/CDS, 138 Other exon/noncoding, and 4 Intergenic peaks. These are strand-aware representative-summit annotations. Ambiguous UTR/CDS means overlap with more than one of 5' UTR, CDS, and 3' UTR annotation classes; Other exon/noncoding means exonic overlap without any of those primary classes. Intergenic means no same-strand transcript overlap after the preceding annotation rules.

Each non-reference category is compared with the full 3' UTR distribution using a two-sided two-sample Kolmogorov-Smirnov test. R's stats::ks.test is run with exact=FALSE, simulate.p.value=TRUE, B=99999, seed 20260928. Conditional Monte Carlo simulation preserves pooled tied values. BH correction covers exactly the three comparisons. Plot labels are BH-adjusted q-values; 3' UTR is labeled Reference. The horizontal zero line is only a visual guide, not the statistical null used in these comparisons.

Results (raw p; BH q; KS D): CDS 0.00999; 0.014985; 0.08189175. Intron 0.00001; 0.00003; 0.16474554. Other 0.15145; 0.15145; 0.05608370. Intron reaches the Monte Carlo minimum p of 1/(99999+1); the tail probability is not resolved more finely. KS D is the maximum absolute difference between empirical cumulative distribution functions.

These are exploratory peak-level comparisons, not biological-replicate tests. Correlation among peaks from the same gene can affect inferential calibration. Failure to reject for Other is not evidence of identical distributions. Results apply to the selected peak set.

Statistics, category counts/medians, Other annotation counts, source data, PNG (600 dpi), and vector PDF are saved beside this file. No sample counts are displayed on the plot.

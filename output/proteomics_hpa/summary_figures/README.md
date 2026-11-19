# HPA localization analysis: summary figures

These figures visualize the completed analyses without rerunning or changing any statistical test. The conclusion is **no statistically significant HPA localization enrichment or between-set difference detected**, not proof that the two isoforms have identical cofactors or equivalent localization profiles.

## Figure guide

1. **`01_analysis_overview`** — The different questions and analysis choices in one overview. Each small dot is the minimum category FDR within one correction family; the large dot is the minimum across families in that row. Grey spans show the range of family minima, not confidence intervals. The shaded region is FDR <0.10; dotted/dashed lines mark 0.05/0.10. No result crosses either threshold. Includes all four ranked GSEA analyses and their four exact-zero-exclusion counterparts, eight original HPA ORAs, eight annotated-universe ORAs, 12 exact-Reactome-input ORAs, six broad-category ORAs and three direct-comparison families. Forty-five families represent related analyses of the same experiment, not independent replications. Minimum FDR is not an effect size or a metric for choosing a preferred analysis.
2. **`02_direct_category_proportions`** — Red/blue connected points show original 324/328 proportions for every broad category and all three matched pairs. Multi-label categories overlap, so percentages need not sum to 100%. These descriptive proportions retain all original memberships, including shared EIF4A2 in the specific sets.
3. **`03_direct_odds_ratios`** — The direct-comparison effect estimates and exact unadjusted 95% confidence intervals, on a log odds-ratio axis. Every interval includes 1. The specific-pair panel correctly uses 171/299 nonshared genes for inference; other panels use the original disjoint sets. Intervals convey the differences still compatible with the data. They are not simultaneous/multiplicity-adjusted intervals and do not demonstrate equivalence.
4. **`04_ranked_sensitivity_comparison`** — Collapsed GSEA NES values for all genes versus exact-zero exclusion at both HPA confidence levels. Negative NES points toward 328 and positive toward 324; no result is significant. The high-confidence nonzero Sperm structures estimate is absent because it falls below five genes.
5. **`05_nucleolus_mitochondria_focus`** — Descriptive proportions for the two independently motivated hypotheses, with numerator/denominator labels. Bars are fractions of genes, not means of biological replicates. The direct confidence intervals are shown in Figure 3. The nearly identical EGFP-filtered mitochondrial proportions (15.2% versus 15.4%) give p=1 and FDR=1.

Each figure is available as a high-resolution PNG and vector PDF. **`HPA_localization_summary_all_figures.pdf`** contains all five figures as a single five-page document. For presentation, Figures 1, 2 and 3 provide the overall conclusion, observed proportions and uncertainty; Figure 5 focuses on the specific biological hypotheses.

## Plot data and reproducibility

`all_analysis_families.tsv` and `method_summary.tsv` contain the FDR summaries. `direct_comparison_plot_data.tsv`, `ranked_sensitivity_plot_data.tsv` and `hypothesis_plot_data.tsv` supply the plotted values. `figure_data_checks.tsv` records checks against the saved results. `sessionInfo.txt` records R/package versions.

Run from the repository root:

```sh
Rscript analysis/proteomics_hpa/plot_analysis_summary.R
```

The original analysis folders retain complete statistics, mappings and methods. HPA annotation is prior gene-level localization information, not a localization measurement in this experiment. An absent annotation does not prove biological absence. These figures summarize localization analyses only; they do not negate the separate Reactome functional-pathway findings.

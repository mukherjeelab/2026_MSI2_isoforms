# Additional direct CLIP comparisons

Both versions retain the original denominator (all eligible aligned reads), replicate averaging, colors and y-axis range. No SMI is used.

The intron/3UTR figure is a display subset: no renormalization and no change to the original six-category BH-adjusted q values.

The five-category figure sums intergenic and other fractions within each library before calculating means, effects and exact library-level permutation tests. All five categories are then BH-adjusted together. Other category fractions remain unchanged.

The exact two-sided 2-versus-2 library tests have minimum P = 1/3. The q values are exploratory; a lack of significance is not evidence of equivalence.

Reproduce from the project root: `Rscript output/followup_analysis/read_annotation_comparisons/additional_figures/plot_additional.R`.

Per-library values, means and statistics are saved beside each PNG/PDF. Original figures and source tables are unchanged.

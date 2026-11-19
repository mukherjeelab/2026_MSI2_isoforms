# Unstacked read-annotation comparisons

The two figures reuse the exact fractional aligned-read category table behind the original stacked plot: `output/followup_analysis/tables/read_annotation_summary.tsv`. The annotation method gives each eligible aligned read total weight one, split equally among its distinct overlapping categories. Each sample’s six category fractions sum to one. No reads or peak annotations are recalculated.

## Version 1: CLIP versus pooled SMI

For each category, calculate the unweighted arithmetic mean of the four SMI **library fractions** (two 324 SMI and two 328 SMI libraries). Do not pool their raw reads or weight libraries by sequencing depth. For each isoform, average its two CLIP library fractions and divide by the shared SMI mean. Bars show that fold ratio; white points show each CLIP replicate divided by the same reference. The dashed line at one represents SMI, and each panel prints the reference fraction. There are no SMI bars.

Category statistics compare the two original CLIP library fractions with all four original SMI library fractions. The four SMI values remain separate statistical observations; their average is not treated as four identical replicates or an error-free observation. There are 12 tests (two isoforms x six categories).

White points only show CLIP replicate dispersion around a common denominator. They are not confidence intervals and do not depict uncertainty in the SMI mean. SMI standard deviations and ranges are saved in `pooled_SMI_reference.tsv`.

## Version 2: direct 324 versus 328

Bars show the unweighted mean of the two CLIP library fractions for each isoform within each category; white points show both replicates. No SMI normalization is used. Panel labels report the mean 324-minus-328 difference in percentage points and adjusted P value. Six within-category tests compare the two 324 libraries against the two 328 libraries. No category is tested against a different category.

## Statistical method and limitations

Use the difference in mean library fractions as the test statistic. Enumerate every allocation of whole libraries to groups while preserving group sizes, and calculate the exact two-sided P as the fraction of allocations with an absolute mean difference at least as large as observed, including ties. The observed allocation is included; no Monte Carlo correction or seed is required. Reads are never treated as independent biological replicates and fractional counts are not rounded into binomial read counts.

For 2 CLIP versus 4 SMI libraries there are 15 possible allocations, so the minimum possible P for this absolute-difference test is 1/15 = 0.0667. For 2 versus 2 CLIP libraries there are six allocations and complementary assignments have identical absolute differences, giving minimum two-sided P = 2/6 = 0.3333. Thus these designs cannot attain nominal P <0.05 with this test, even before adjustment. This is a resolution limit, not evidence of equivalence.

Apply Benjamini-Hochberg adjustment separately across the 12 CLIP-versus-SMI tests and six direct-isoform tests. Plots show adjusted values (`q_BH`); raw P values, effect sizes, group means, sample numbers, and every enumerated null statistic are saved in the tables. No tests have q <0.05.

The label-permutation analysis assumes exchangeability of library fractions under the null. Pooling all four SMI libraries, as requested, mixes construct backgrounds and does not preserve the original matched CLIP–SMI pairing; the tests are exploratory rather than a paired analysis. The direct isoform test is also unpaired; shared A/B labels alone were not taken as proof of cross-construct pairing. Category fractions are compositional, so an increase in one category necessarily accompanies decreases elsewhere; these plots do not establish absolute binding changes.

## Main descriptive effects

- Intron: 324 CLIP 42.36%, 328 CLIP 20.80%; difference +21.56 percentage points. Fold over pooled SMI: 1.47x and 0.72x.
- 3′UTR: 324 CLIP 18.93%, 328 CLIP 39.03%; difference -20.11 percentage points. Fold over pooled SMI: 1.42x and 2.92x.

## Files

- `read_categories_CLIP_vs_pooled_SMI.png/pdf`: reference-normalized, nonstacked replicate means.
- `read_categories_324_vs_328.png/pdf`: direct fractional-composition comparison.
- `source_library_fractions.tsv`: exact source values with group labels.
- `pooled_SMI_reference.tsv`: equally weighted SMI mean, SD, range and n for each category.
- `replicate_mean_fractions.tsv`: CLIP means, SD, ranges, reference means, folds and percentage-point effects.
- `CLIP_plot_inputs.tsv`: all individual CLIP points and shared references.
- `CLIP_vs_pooled_SMI_statistics.tsv`, `324_vs_328_statistics.tsv`: all inferential and effect-size results.
- `exact_permutation_distributions.tsv`: every label-assignment statistic.
- `validation.tsv`: independent Python reconstruction of the permutation tests; means, fraction sums and BH adjustment were also independently checked.
- `input_provenance.tsv`, `session_info.txt`, `output_manifest.tsv`: source and output hashes and software versions.

## Reproduce

From the project root:

```sh
Rscript output/followup_analysis/read_annotation_comparisons/analyze_read_annotation.R
```

Requires R, data.table, ggplot2 and scales. The source file hash is checked before and after execution. Both PNG figures were visually inspected. Existing figures are preserved.

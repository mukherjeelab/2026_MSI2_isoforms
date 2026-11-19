# MSI2 isoform-associated proteins: HPA localization results

**No compartment passes BH FDR < 0.05.** The ranked analyses and exploratory specific-set analyses do not establish a shift between nuclear and cytoplasmic/local membrane-associated environments for MSI2 324 versus 328.

The analysis uses HPA 25.1 and all resolvable genes from the supplied quantitative comparison, without filtering on differential p-values. Of 1,677 tested protein rows, 1,654 (98.63%) map to HPA; collapsing 13 duplicate-gene pairs by median log2 fold change gives **1,641 ranked genes**. All 23 unresolved rows remain documented. Primary-confidence localization is available for 1,439 unique ranked genes.

## Ranked enrichment

| Analysis | Categories tested | Significant at FDR <0.05 | Smallest FDR |
|---|---:|---:|---:|
| Fine HPA, Enhanced/Supported/Approved | 41 | 0 | 0.857 |
| Collapsed HPA, Enhanced/Supported/Approved | 15 | 0 | 0.827 |
| Fine HPA, Enhanced/Supported | 29 | 0 | 0.781 |
| Collapsed HPA, Enhanced/Supported | 15 | 0 | 0.741 |

Primary collapsed results for the central biological questions:

| HPA compartment | Measured genes | NES | FDR |
|---|---:|---:|---:|
| Nucleoplasm, including nuclear speckles/bodies | 626 | -0.968 | 0.951 |
| Nucleoli | 122 | -1.079 | 0.827 |
| Nuclear membrane | 25 | +0.763 | 0.968 |
| Cytosol | 599 | -0.953 | 0.951 |
| Endoplasmic reticulum | 105 | -1.127 | 0.827 |
| Golgi apparatus | 68 | -0.888 | 0.951 |
| Vesicular system | 168 | -1.167 | 0.827 |
| Plasma membrane/cell junctions | 160 | -0.869 | 0.968 |
| Actin-filament associated | 35 | -0.681 | 0.968 |
| Microtubule associated | 54 | -0.952 | 0.951 |
| Centrosome | 51 | -0.912 | 0.951 |
| Intermediate filaments | 21 | +0.716 | 0.968 |
| Mitochondria | 186 | -1.137 | 0.827 |

Positive NES points toward 324-associated proteins and negative NES toward 328-associated proteins, but **none of these directions constitutes significant enrichment**. Nucleoplasm and Cytosol both have weak negative scores; they do not show an opposing nuclear/cytoplasmic pattern. In the fine analysis, Nuclear speckles also has a weak negative score (NES -0.904; FDR 0.857), while Nuclear bodies is near the positive null expectation (NES +0.999; FDR 0.857). The broad Nucleoplasm group is retained despite its size.

The Enhanced/Supported sensitivity analysis agrees on the absence of significant enrichment. Four collapsed categories change NES sign: actin-filament associated, Golgi, intermediate filaments and plasma membrane/cell junctions. This reinforces that their primary-analysis direction should not be treated as a robust biological preference. The sensitivity comparison is saved category by category.

## Exploratory isoform-specific ORA

The historical notebook rules reproduce **171 324-specific genes** and **298 328-specific genes** in the mapped experimental universe. These are exploratory classifications: the rule uses k-means clusters, a twofold cutoff and nominal p <0.5. One gene with conflicting duplicate-row classifications is withheld from both hit sets and retained in the 1,641-gene background.

No localization is significant in any ORA family. The most notable nonsignificant observation among the requested compartments is ER annotation in the 328-specific set: **28/298 hits**, compared with **105/1,641 background genes**, odds ratio **1.70**, nominal p **0.0166**, collapsed FDR **0.249** (fine-grained FDR **0.680**). This does not support a corrected-significant ER association.

There are no significant GSEA leading edges to interpret. The significant-leading-edge table is therefore empty by design; leading edges for every tested category are still supplied for transparent inspection. Running plots show selected nonsignificant patterns from each direction and explicitly label them as such.

## Additional sensitivity checks

Both requested follow-up checks retain the no-significant-enrichment conclusion. Removing all **579 exact-zero gene scores** leaves **1,062 genes** for the four GSEA reruns. No category reaches FDR <0.05; minimum FDRs across the four collections range from **0.701 to 0.865**. Categories falling below five retained genes are recorded as untested.

Restricting ORA to the **1,439 genes with at least one primary-confidence HPA localization** leaves **155 324-specific** and **272 328-specific** genes. None of the eight ORA families yields a significant category. The 328-specific ER observation weakens from odds ratio 1.70 and collapsed FDR 0.249 to odds ratio **1.62** and FDR **0.406**. The all-measured/mapped universe remains primary.

See [the sensitivity report](sensitivity/README.md) for methods, category-level comparisons, restricted hit sets and six additional PDF/PNG figure pairs. These checks show robustness of the FDR conclusion to the two specified choices, not absence of all possible annotation or processing biases.

## Interpretation limits

This analysis evaluates whether the cofactor ranking is enriched for proteins **annotated by HPA to localize to** each compartment. HPA is gene-level prior knowledge from selected human cell lines, not a direct localization measurement in this experiment and generally not isoform-specific. The lack of significant enrichment does not prove absence of a localization difference.

The workbook's original QC, model and imputation details are not available. All supplied tested rows were retained through the identifier audit; 579 ranked genes have exactly zero log2 fold change. The upstream reason for these tied values is unknown. Methods, source hashes, complete mapping, gene sets, results and validation are documented in [README.md](README.md).

Main figures: [collapsed NES](figures/NES_collapsed_primary.pdf), [fine NES](figures/NES_fine_primary.pdf), and [cross-analysis comparison](figures/analysis_tile_collapsed.pdf).

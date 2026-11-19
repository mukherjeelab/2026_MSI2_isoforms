# Exploratory ORA of eight broad HPA localization categories

**No category passes BH FDR <0.05 or <0.10.** The smallest FDR is **0.18952**, for Mitochondrial annotations in the original EGFP-filtered 328 hit set: **37/240 hits** versus **186/1,638 background genes**, odds ratio **1.5274**, nominal p **0.02369**.

## Fixed inputs and exploratory scope

This follow-up uses all three paired Reactome inputs from the preceding matched-input analysis, without reselecting or reclassifying any gene. The Entrez-mapped hit sizes remain 421/609 for the full 324/328 clusters, 237/240 for the EGFP-filtered clusters, and 172/300 for the isoform-specific subsets. Every test uses the same **1,638 measured Entrez genes**. EIF4A2 remains in both original isoform-specific sets. Unknown HPA annotations and the ambiguous MATR3 mapping are handled exactly as before; original hit and background memberships are retained.

Inputs were copied unchanged to `unchanged_hit_memberships.tsv` and `unchanged_measured_background.tsv`; input checksums are recorded. These inputs were reconstructed from the current Reactome notebook code, as documented in the [matched-input report](../reactome_matched_ora/README.md), because historical saved Reactome input objects were unavailable.

This is an **additional exploratory analysis** motivated by the user's independent Reactome hypotheses. The eight-category mapping was fixed and saved before calculating enrichment. It is a user-requested aggregation of HPA terms, not an official new HPA ontology and not a replacement for the earlier analyses. Correcting within eight categories does not adjust for the full collection of exploratory analyses performed across this project.

## Broad category definitions

Only HPA 25.1 **Enhanced, Supported and Approved** annotations are used, as in the previous matched-input ORA. Uncertain annotations are excluded.

| Broad category | Included HPA terms or previously documented HPA parent groups |
|---|---|
| Nuclear/nucleoplasmic | Nucleoplasm, nuclear bodies, nuclear speckles, kinetochore, mitotic chromosome, nuclear membrane |
| Nucleolus | Nucleoli, nucleoli fibrillar center, nucleoli rim |
| Cytosolic | Cytosol, cytoplasmic bodies, aggresome, Rods & Rings |
| Mitochondrial | Mitochondria only |
| ER/Golgi | Endoplasmic reticulum and Golgi apparatus |
| Vesicular/endolysosomal | HPA vesicular group: vesicles, endosomes, lysosomes, lipid droplets and peroxisomes |
| Plasma membrane/junctional | Plasma membrane and Cell Junctions |
| Cytoskeletal | Actin-filament associated, microtubule associated, intermediate filaments, centrosome/centriolar satellite, and primary-cilium groups |

The Cytoskeletal category includes cleavage furrow, focal adhesion sites, cytokinetic bridge, microtubule ends, midbody, midbody ring, mitotic spindle, basal body, primary cilium, its tip and transition zone, in addition to the named structural groups. Sperm-specific HPA terms are retained in the mapping audit but excluded from this eight-category scheme. A gene with only excluded annotations remains in the measured background and any original hit set.

The vesicular label is intentionally broad: lipid droplets and peroxisomes follow the previously used HPA vesicular grouping, and should not be interpreted as exclusively endolysosomal structures. Nucleolus is not absorbed into Nuclear/nucleoplasmic. Mitochondria is not combined with any other category. Genes annotated to multiple fine terms within one broad category count once; genes annotated to different broad categories remain in each applicable category. A gene can therefore occur in both Nuclear/nucleoplasmic and Nucleolus if HPA independently annotates both.

`fine_to_broad_mapping.tsv` records all 49 terms, their assigned broad category or explicit exclusion, the prior HPA parent and source URL. `hpa_broad_gene_sets.tsv` records the full gene sets; `measured_gene_broad_membership.tsv` records measured memberships with the contributing fine terms.

## Statistical analysis

For each of the six original hit sets, conduct eight one-sided Fisher exact enrichment tests (`alternative='greater'`) against the fixed background. The 2x2 cells are overlap, hit-set size minus overlap, background category count minus overlap, and the remaining background. Odds ratios are R's conditional maximum-likelihood estimates. Nominal p-values were independently checked against the hypergeometric upper tail.

**BH correction is applied across exactly eight categories separately for each original hit set**: each pair's 324 and 328 analyses have separate correction families. This yields six families and 48 tests in total. All eight categories contain at least five measured genes and less than 95% of the background, so no category is removed. No additional annotation-coverage or fold-change filter is applied to hit or background genes.

## Results by input pair

| Original pair | 324 minimum FDR | 328 minimum FDR | Categories below FDR 0.10 |
|---|---:|---:|---:|
| Full clusters | 0.240 | 0.663 | 0 |
| EGFP-filtered clusters | 0.254 | 0.190 | 0 |
| Isoform-specific subsets | 0.792 | 0.513 | 0 |

### Nucleolus and Mitochondrial hypotheses

The background counts are **122/1,638 for Nucleolus** and **186/1,638 for Mitochondrial**. These category memberships, overlaps and nominal p-values are identical to the preceding matched-input analysis; their BH FDRs change because the correction family is now the eight broad categories.

| Pair | Isoform | Category | Hit overlap | Odds ratio | Nominal p | BH FDR |
|---|---|---|---:|---:|---:|---:|
| Full clusters | 324 | Nucleolus | 34/421 | 1.127 | 0.3181 | 0.6510 |
| Full clusters | 328 | Nucleolus | 49/609 | 1.146 | 0.2690 | 0.6804 |
| EGFP-filtered | 324 | Nucleolus | 18/237 | 1.025 | 0.5054 | 0.9831 |
| EGFP-filtered | 328 | Nucleolus | 16/240 | 0.871 | 0.7309 | 0.9861 |
| Isoform-specific | 324 | Nucleolus | 13/172 | 1.018 | 0.5242 | 0.7917 |
| Isoform-specific | 328 | Nucleolus | 26/300 | 1.227 | 0.2186 | 0.5830 |
| Full clusters | 324 | Mitochondrial | 59/421 | 1.399 | 0.0300 | 0.2404 |
| Full clusters | 328 | Mitochondrial | 74/609 | 1.132 | 0.2410 | 0.6804 |
| EGFP-filtered | 324 | Mitochondrial | 36/237 | 1.493 | 0.0317 | 0.2538 |
| EGFP-filtered | 328 | Mitochondrial | 37/240 | 1.527 | 0.0237 | 0.1895 |
| Isoform-specific | 324 | Mitochondrial | 19/172 | 0.966 | 0.5938 | 0.7917 |
| Isoform-specific | 328 | Mitochondrial | 40/300 | 1.256 | 0.1375 | 0.5500 |

The mitochondrial nominal enrichments in both EGFP-filtered sets and the full 324 cluster do not survive correction. The broad HPA analysis therefore provides no FDR-significant support for either localization hypothesis. This does not refute specific Reactome rRNA-processing or mitochondrial-translation pathway results: localization categories test broader, different gene memberships and biological questions. No claim is made here that the independently reported Reactome findings have been revalidated.

## Deliverables and reproduction

- `broad_hpa_ora_all.tsv`: all 48 tests with overlap/background counts, GeneRatio/BgRatio, odds ratios, p-values, BH FDR and contributing source symbols, Entrez IDs and HPA IDs. Six `ora_genes_*.tsv` files provide the individual input-set results.
- `contributing_genes_long.tsv`: one contributor per category/set with the underlying fine HPA annotations.
- `nucleolus_mitochondria_hypotheses.tsv`: dedicated table for the two independently motivated hypotheses.
- `analysis_summary.tsv`, `fine_to_broad_mapping.tsv`, gene-set/membership files, unchanged input copies, checksums and software versions support the audit.
- `figures/`: three paired PNG/PDF dot plots, one per original input pair.
- `validation_checks.tsv`: exact input preservation, 8-test BH families, Fisher/hypergeometric agreement, preserved Nucleolus/Mitochondrial counts and p-values, and multi-label membership checks. All passed. Significant-result files are header-only because no result passes the indicated threshold.

Run from the repository root after the matched-input analysis:

```sh
Rscript analysis/proteomics_hpa/run_broad_exploratory_ora.R
```

HPA remains gene-level prior localization knowledge, not evidence of localization in the experimental cell type or of a particular protein isoform. Missing annotations remain unknown. Earlier fine-grained, collapsed and sensitivity outputs are preserved separately.

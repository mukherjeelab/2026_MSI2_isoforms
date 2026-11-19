# HPA ORA using the exact existing Reactome input sets

All three paired Reactome input definitions in `code/Proetomics.Rmd` were tested against fine-grained and collapsed HPA localization categories. **No category reaches BH FDR <0.05 or <0.10 in any of the 12 analyses.**

The lowest FDR is **0.27799**, for **Endoplasmic reticulum in the original 328-specific set**, collapsed analysis: **28/300 hits**, **105/1,638 background genes**, conditional odds ratio **1.685**, nominal p **0.01853**. This is not significant after correction.

## Input sets: no redefinition

The original notebook's clean-data, heatmap/clustering and hit-selection code was executed verbatim, rather than implementing new selection rules. The executed code fragments and all six original symbol vectors are saved in `../reactome_input_audit/`. No archived historical Reactome result object or saved historical gene vectors were found, so these are exact reconstructions from the current notebook code and installed identifier database, not a claim to recover an unknown earlier interactive R session.

| Notebook pair | Original 324 symbols | Original 328 symbols | Mapped 324 Entrez genes | Mapped 328 Entrez genes |
|---|---:|---:|---:|---:|
| `genes_324_cluster` / `genes_328_cluster` | 421 | 611 | 421 | 609 |
| `genes_324_hits` / `genes_328_hits` | 237 | 241 | 237 | 240 |
| `genes_324_specific` / `genes_328_specific` | 172 | 301 | 172 | 300 |

The full-cluster pair comes from the notebook's three-condition, row-scaled k-means (k=3, seed 123), taking cluster 2 for 324 and cluster 3 for 328. The EGFP-filtered pair uses those respective clusters plus the notebook's matching isoform/EGFP fold change >=2 and nominal p <0.5. The isoform-specific pair uses the union of clusters 2 and 3, the original 324/328 fold-change direction with >=twofold effect and 0 < nominal p <0.5. These thresholds describe the original inputs; they were not tuned for HPA.

**EIF4A2 is in both original isoform-specific sets and remains in both analyses.** The previous HPA workflow excluded a gene with conflicting protein-row directions; that rule is not used here. The other two pairs have no intersection between their mapped 324 and 328 hit sets. See `between_set_overlap_preserved.tsv`.

## Exact measured universe and identifier mapping

The notebook defines its background as `unique(prot_df$gene)` and maps symbols with `clusterProfiler::bitr(..., fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)`. That function is reused directly. It yields **1,643 original symbols and 1,638 unique Entrez IDs**; all original symbols and mapping failures are recorded. The installed `org.Hs.eg.db` version is **3.20.0**.

**Every Fisher test uses the same 1,638-gene background.** This is the experimentally measured, mapped universe supplied to the Reactome function. It is not restricted to the subset having Reactome pathway annotations or HPA localization annotations. Enrichment software can internally intersect supplied universes with its own annotation coverage; such a Reactome-specific subset is not substituted for the experiment's supplied universe here.

The HPA 25.1 reference and annotation files are the same frozen downloads used in the earlier analysis. HPA genes are connected to the original Entrez IDs through `org.Hs.eg.db` Ensembl mappings, with concordant exact symbols to resolve multiple candidates and exact-symbol fallback where no Ensembl match exists. The complete bridge, candidate IDs, mapping method and coverage flags are saved.

- **1,637/1,638** measured Entrez genes resolve uniquely to HPA genes.
- **1,437** have at least one primary-confidence HPA localization.
- **201** have no usable primary localization, including one unresolved HPA mapping. They remain in the background and any original hit sets.
- **MATR3 (Entrez 9782)** maps to two candidate HPA Ensembl genes and cannot be uniquely resolved by symbol; its annotation is left unknown. Its original set memberships are retained. No annotations are propagated arbitrarily to a selected gene.

This Entrez-based input reconstruction differs from the earlier workbook-Ensembl/HPA-first workflow's 1,641-gene universe. It also preserves original set overlaps and does not add newly rescued symbol-free workbook rows. These changes are required to match the existing Reactome inputs rather than redefine them.

## HPA categories and tests

The analysis uses **Enhanced + Supported + Approved** localization assignments; Uncertain is excluded. Multi-localization membership is preserved: one Entrez gene can count in several categories, but only once within a given category. Fine categories and the 15 collapsed groups use the previously frozen HPA hierarchy mapping.

As in the preceding HPA analysis, categories must contain **at least five measured background genes** and **less than 95% of the measured universe**. This retains 41 fine terms and all 15 collapsed compartments. All 49 fine terms remain in the output, with ineligible categories explicitly labeled and statistics set to NA. No arbitrary 500-gene upper limit is imposed on large compartments such as Nucleoplasm. No filters are applied to the hit sets themselves beyond the original notebook's Entrez mapping.

For each pair, each isoform and each resolution, a one-sided Fisher exact test (`alternative='greater'`) tests enrichment using:

| | In category | Outside category |
|---|---:|---:|
| In original hit set | a | n-a |
| Other measured genes | K-a | N-n-K+a |

Here N=1,638, n is the unchanged mapped hit-set size, K is the background category size, and a is the overlap. The nominal p-value is independently checked against the hypergeometric upper tail. Odds ratios are Fisher's conditional maximum-likelihood estimates.

The main `BH_FDR` is corrected **separately for each original input set and resolution**, yielding 12 families: three pairs x two isoforms x two resolutions. No 324 and 328 tests are pooled. `BH_FDR_both_resolutions` additionally supplies correction across the fine and collapsed tests together within each original set (six families), for readers who prefer that scope. Primary interpretation uses the explicitly labeled per-resolution correction.

## Results

Minimum BH FDR in each analysis:

| Original pair | 324 fine | 324 collapsed | 328 fine | 328 collapsed |
|---|---:|---:|---:|---:|
| Full clusters | 0.469 | 0.451 | 0.598 | 0.378 |
| EGFP-filtered clusters | 0.617 | 0.476 | 0.971 | 0.355 |
| Isoform-specific subsets | 0.965 | 0.794 | 0.760 | 0.278 |

The 328-specific ER observation is the smallest corrected p-value overall. For EGFP-filtered 328 hits, Mitochondria has 37/240 hits versus 186/1,638 background genes, odds ratio 1.527, p=0.02369 and collapsed FDR=0.355. For the full 328 cluster, ER has 49/609 hits versus 105/1,638 background genes, odds ratio 1.520, p=0.02517 and collapsed FDR=0.378. These are nonsignificant descriptive observations, not established localization preferences.

The matched-input comparison therefore does not uncover significant HPA localization enrichment for any of the existing Reactome gene sets. Reactome pathways and HPA localizations address different annotation questions; this result does not establish or refute the original Reactome pathway findings.

## Deliverables

- **`hpa_ora_all_exact_reactome_sets.tsv`**: complete table plus 12 individual `ora_*.tsv` tables. Columns include overlap count, hit-set/background sizes, GeneRatio and BgRatio as counts and decimals, odds ratio, nominal p, BH FDR, status, contributing original symbols, Entrez IDs and HPA genes.
- `analysis_summary.tsv`: per-family counts, minimum p/FDR and number below 0.05/0.10.
- `original_symbol_hit_memberships.tsv`, `exact_mapped_hit_memberships.tsv`, `original_measured_universe.tsv`: exact auditable inputs and mapping failures.
- `measured_universe_hpa_mapping.tsv`, `orgdb_entrez_ensembl_bridge.tsv`, `measured_gene_localization_membership.tsv`: mapping and multi-label category memberships.
- `hit_mapping_QC.tsv`, `universe_QC.tsv`, `unannotated_or_unmapped_universe_genes.tsv`, `between_set_overlap_preserved.tsv`: coverage and overlap audit.
- `figures/`: six paired PDF/PNG dot plots, one fine and one collapsed comparison for each original pair. X is GeneRatio, size is overlap, color is FDR. All eligible categories are shown.
- `validation_checks.tsv`: exact membership, universe inclusion, numerical Fisher/hypergeometric and BH checks. `input_checksums.tsv` and `sessionInfo.txt` preserve provenance. `significant_localizations.tsv` contains a header and zero rows because no category is significant.

Run from the repository root:

```sh
Rscript analysis/proteomics_hpa/audit_reactome_inputs.R
Rscript analysis/proteomics_hpa/run_reactome_matched_ora.R
```

The audit loads the original notebook's packages and reconstructs inputs; the HPA script uses frozen local annotation files and does not rerun Reactome enrichment or download new data. Earlier HPA results remain separate and unchanged.

HPA annotations are gene-level prior localization knowledge, not a measurement of location in this experiment and generally not protein-isoform-specific. Unknown annotation remains unknown. The exploratory notebook thresholds and all limitations of reconstructing inputs from code rather than archived historical objects apply.

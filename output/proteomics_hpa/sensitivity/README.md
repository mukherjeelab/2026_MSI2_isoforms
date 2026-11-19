# Sensitivity analyses: exact-zero ranks and incomplete localization annotation

Both requested sensitivities preserve the primary conclusion: **no HPA localization category passes BH FDR <0.05**. These checks address two specific analysis choices; they do not establish that incomplete annotation or unknown upstream processing has no other effect.

## 1. Excluding exact-zero gene scores from GSEA

Removed the **579 genes with exactly zero gene-level log2(324/328)** after the original median duplicate-gene collapse. This leaves **1,062 genes**, with their original scores and relative ordering unchanged. No near-zero tolerance, p-value filter, re-estimation of fold changes, or additional HPA-coverage filter was introduced. The excluded genes, retained ranking and RNK file are saved.

All four HPA collections were rerun with the original weighted pre-ranked fgseaMultilevel settings: seed 20260924, gseaParam=1, eps=0, sampleSize=1001, nPermSimple=10000, serial execution. Gene sets were re-intersected with the reduced ranking. Eligibility remains at least five retained genes and less than 95% of the current ranked universe. BH correction is recalculated within each collection across its eligible categories.

| HPA collection | Primary categories | Nonzero categories | Nonzero minimum FDR | Significant |
|---|---:|---:|---:|---:|
| Fine, primary confidence | 41 | 40 | 0.865 | 0 |
| Collapsed, primary confidence | 15 | 15 | 0.701 | 0 |
| Fine, Enhanced/Supported | 29 | 24 | 0.799 | 0 |
| Collapsed, Enhanced/Supported | 15 | 14 | 0.720 | 0 |

Seven collection/category combinations fall below five measured genes: primary fine Peroxisomes; high-confidence fine Cytokinetic bridge, Endosomes, Focal adhesion sites, Peroxisomes and Primary cilium tip; high-confidence collapsed Sperm structures. They are explicitly marked untested with NA statistics, not treated as negative results. They remain in the complete comparison table.

Among categories tested in both runs, only the primary collapsed Actin-filament associated group changes NES sign (-0.681 to +0.663; sensitivity FDR 0.956). This weak nonsignificant direction change does not support an actin-related preference. No nuclear/cytoplasmic shift becomes significant after removing zeros. The significant-leading-edge table is empty by design, and complete leading edges remain available.

## 2. ORA restricted to genes with primary-confidence HPA annotation

The sensitivity universe contains the **1,439 measured/mapped genes with at least one Enhanced, Supported or Approved localization**. The original 1,641-gene measured/mapped universe remains primary. The restricted universe is fixed for all four localization collections, including Enhanced/Supported sets, so this analysis isolates the requested background restriction rather than adding a different high-confidence universe rule.

Hit sets are intersected with that same 1,439-gene universe. No gene is reclassified, and exact-zero exclusion is **not** applied to ORA. The one gene with conflicting duplicate classifications stays in the background but not either hit set.

| Classification | Original | Annotated universe | Removed |
|---|---:|---:|---:|
| 324-specific | 171 | 155 | 16 |
| 328-specific | 298 | 272 | 26 |
| Shared/non-specific | 1,171 | 1,011 | 160 |
| Conflicting duplicate classifications | 1 | 1 | 0 |
| Total background | 1,641 | 1,439 | 202 |

The original one-sided Fisher enrichment test, category eligibility and per-collection/per-isoform BH correction are retained. All eight ORA families were rerun. No category is significant; minimum FDRs range from 0.406 to 1.000. Eligible category counts remain unchanged because every member of any tested HPA confidence set is in the primary-annotated universe.

The previously noted 328-set ER result becomes weaker:

| ER metric, primary-confidence collapsed ORA | All mapped | Annotated only |
|---|---:|---:|
| Hit overlap | 28/298 | 28/272 |
| Background overlap | 105/1,641 | 105/1,439 |
| Conditional odds ratio | 1.704 | 1.624 |
| Nominal p | 0.0166 | 0.0271 |
| BH FDR | 0.249 | 0.406 |

The annotated-only background does not uncover a significant association that was hidden by inclusion of unknown-localization genes. A protein annotated to at least one compartment may still have unrecorded additional localizations, so this restriction does not eliminate all annotation incompleteness.

## Saved outputs and reproduction

- `gsea_nonzero_all.tsv` and four individual collection tables contain complete results, including untested terms, HPA/retained sizes, ES/NES, p/FDR and leading-edge identifiers/proteins.
- `gsea_primary_vs_nonzero.tsv` compares category size, eligibility, ES, NES, p, FDR and direction to the primary runs.
- `ora_annotated_all.tsv` and eight individual tables contain the restricted-universe ORA results and contributing genes/proteins.
- `ora_primary_vs_annotated.tsv` compares counts, odds ratios and p/FDR against the primary ORA.
- Ranked genes, excluded exact-zero genes, annotated background, restricted/excluded classifications, hit-set counts, expanded leading edges and analysis summary tables are included.
- Six paired PDF/PNG figures show four GSEA comparisons and two primary-confidence ORA comparisons. High-confidence ORA results are fully tabulated. Asterisks mark FDR <0.05; none occur.
- `validation_checks.tsv` records rank/universe checks, independent running-sum ES calculations, Fisher/hypergeometric agreement and checksums confirming existing primary outputs were not modified during the run. `sessionInfo.txt` records software versions.

Run from the repository root after the primary analysis:

```sh
Rscript analysis/proteomics_hpa/run_sensitivities.R
```

The script reads the frozen primary gene sets, scores and classifications and writes only into this sensitivity directory. It uses no new downloads. The original HPA 25.1 provenance, mapping procedures, exploratory specificity threshold (nominal p <0.5) and biological caveats remain applicable. HPA annotations are gene-level prior localization knowledge, not localization measurements in this experiment.

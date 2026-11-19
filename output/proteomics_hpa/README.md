# MSI2 324 versus 328: HPA localization enrichment

## Result

The [direct 324-versus-328 proportion comparison](direct_localization_comparison/README.md) tests all eight broad categories within each original Reactome pair. No comparison is significant; the lowest FDR is 0.498. The EGFP-filtered mitochondrial example gives 36/237 versus 37/240, two-sided Fisher p = 1.00.

The additional [eight-category exploratory ORA](broad_exploratory_ora/README.md) preserves all three original Reactome input pairs and their background. No result passes FDR 0.10; the minimum is 0.190 for Mitochondrial annotations in the EGFP-filtered 328 set.

For the follow-up using all three **exact original Reactome input pairs**, see [reactome_matched_ora/README.md](reactome_matched_ora/README.md). That separate analysis preserves the notebook's original Entrez-mapped universe and set overlaps; it does not redefine the inputs using the HPA-first classifications below.

**No HPA localization was significantly enriched at BH FDR < 0.05** in any of the four ranked GSEA analyses or eight exploratory ORA analyses. These data do not establish a nuclear-versus-cytoplasmic or membrane-associated cofactor shift between MSI2 324 and 328. This is a failure to detect enrichment, not evidence that the isoforms have identical interaction environments.

See [RESULTS.md](RESULTS.md) for the concise biological summary, [gsea_all.tsv](gsea_all.tsv) and [ora_all.tsv](ora_all.tsv) for every category, including categories not tested. Figures are in [figures](figures), with 14 paired PDF/PNG outputs.

## Sources and scope

- Proteomics: `data/FredHutchProteomics.xlsx`, Sheet1, 1,677 rows and 32 columns. Every row has a finite positive 328/324 fold change and a valid comparison p-value. The workbook provides gene symbols, Ensembl IDs and labels that often contain UniProt accessions.
- HPA **25.1**, using Ensembl **109**, as stated on the downloaded [HPA subcellular data page](https://www.proteinatlas.org/humanproteome/subcellular/data). The 49-term lookup and gene annotations were downloaded on 2026-09-24. The downloaded archives and provenance HTML are retained unchanged in `data/hpa_25.1/`. `source_manifest.tsv` records URLs, SHA256 hashes, file sizes and source file modification times. The complete HPA reference `proteinatlas.tsv.zip` was also downloaded to distinguish absent localization annotation from absent gene mapping and to rescue exact UniProt/symbol matches.
- Exploratory classifications: the existing `code/Proetomics.Rmd` notebook. The original notebook and source workbook were not modified.
- **Upstream quality filtering cannot be independently reconstructed from this repository.** No QC pass/fail field, peptide thresholds, identification FDR, differential-model specification, or missing-value/imputation provenance is supplied with this table. We therefore treat the complete supplied comparison table as the available tested universe, and introduce no new significance, abundance, contaminant or EGFP-enrichment filter for GSEA. All 1,677 rows remain in the protein-level audit.
- Other comparison columns containing -1/1 are categorical flags, not usable continuous differential test statistics. The quantitative fold-change column is used directly; the model is not refitted.

## Identifier mapping

`protein_identifier_mapping.tsv` preserves each workbook row number, protein label, source identifiers, candidate matches, mapping method, rank inclusion and localization coverage. `unmapped_or_excluded_proteins.tsv` explicitly lists the 23 excluded rows.

1. Strip Ensembl version suffixes and tokenize semicolon/comma-delimited IDs. Only human `ENSG` IDs are considered for human Ensembl matching.
2. Accept a unique HPA Ensembl match if no exact-symbol match conflicts with it. Multiple HPA Ensembl candidates can be resolved only when an exact source symbol selects one of them.
3. If no HPA Ensembl candidate exists, use a unique exact HPA symbol; if neither exists, try the label's exact UniProt accession against the complete HPA reference.
4. Conflicting identifiers or unresolved multiple candidates are excluded from gene-level tests and retained in the audit. No fuzzy symbol matching or automatic alias replacement is used. A unique source human Ensembl ID absent from HPA would remain in GSEA as an unannotated gene; none required that fallback here.
5. Collapse multiple mapped rows for the same gene by the **median log2(324/328)**. This avoids selecting an extreme effect or the smallest p-value. Save every contributing protein label and source row; min/max row effects are included in the gene ranking. Sort equal scores by Ensembl ID without jitter.

Mapping QC:

| Quantity | Count |
|---|---:|
| Tested protein rows | 1,677 |
| Rows mapped to an HPA gene | 1,654 (98.63%) |
| Rows with any HPA localization, including Uncertain | 1,473 (87.84%) |
| Rows without any HPA localization, including unresolved rows | 204 (12.16%) |
| Rows with primary-confidence localization | 1,452 (86.58%) |
| Rows with Enhanced/Supported localization | 1,026 (61.18%) |
| Input rows with multiple IDs | 85 |
| Input rows with multiple human IDs | 82 |
| Unresolved ambiguous/conflicting rows | 3 |
| Other unresolved rows | 20 |
| Genes represented by two mapped protein rows | 13 |
| Unique genes in GSEA and ORA background | 1,641 |
| Unique genes with any localization | 1,460 |
| Unique genes with primary-confidence localization | 1,439 |

An absent HPA localization is **unknown**, not evidence that a protein cannot occupy that compartment. The 202 genes without primary localization remain in the ranked universe and measured HPA-mapped ORA background. This choice is documented because annotation coverage can influence enrichment tests.

## HPA gene sets

`hpa_gene_location_long.tsv` preserves the multi-label gene/location/confidence structure, including Uncertain annotations for audit. Primary sets use the union of **Enhanced, Supported and Approved**. Sensitivity sets use **Enhanced and Supported** only. Neither uses Uncertain annotations.

The collapsed mapping is parsed from the archived [HPA subcellular overview](https://www.proteinatlas.org/humanproteome/subcellular), not inferred from GO parent terms. `collapsed_location_mapping.tsv` records all 49 fine terms, 15 collapsed compartments, original HPA parent, GO term and source URL. Nucleoplasm includes nuclear bodies, speckles, kinetochore and mitotic chromosome; Cytosol includes aggresome, cytoplasmic bodies and Rods & Rings; other groups follow the published hierarchy. The vesicular system includes vesicles, endosomes, lysosomes, lipid droplets and peroxisomes. Sperm structures are retained for completeness, but are not emphasized in the interpretation or selected running plots.

Each of the four collections is saved as both a long TSV and an Ensembl-ID GMT (`gene_sets_*`). Multiple fine terms within a collapsed group count a gene only once. Genes can belong to multiple different groups.

## Ranked GSEA

- Ranking: `-log2(f_MSI2_328_WT /MSI2_324_WT)`, so positive favors 324 and negative favors 328. No test-statistic column is available. No p-value cutoff is applied.
- All 1,641 resolved unique genes are ranked. 595 source rows have exactly unit fold change; after mapping/collapsing, **579 genes have zero score** (35.28%). They remain in the ranking. The supplied table does not explain whether these values reflect upstream imputation; we do not infer that they are biological equality measurements.
- Weighted, two-sided pre-ranked GSEA using `fgsea::fgseaMultilevel`, `gseaParam=1`, `eps=0`, `sampleSize=1001`, `nPermSimple=10000`, seed 20260924 and a serial BiocParallel backend. Software versions are saved in `sessionInfo.txt`. See the [fgsea documentation](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html).
- Test a set if it contains **at least 5 measured genes** and **less than 95% of the ranked universe**. No arbitrary fixed upper size limit. Large groups such as Nucleoplasm are retained. In this dataset no category reaches the 95% limit.
- BH correction is applied within each confidence-by-resolution collection, across eligible localization categories. Untested categories have NA statistical results and explicit status. The correction is not pooled across the four complementary analyses.
- Positive NES indicates the positive/324 end of the rank list; negative NES indicates the negative/328 end. Direction alone is not evidence for enrichment. NES near or below 1 in magnitude and high FDR should not be described as an established compartment preference.
- Complete results contain HPA size, measured size, ES, NES, nominal p-value, BH FDR, Monte Carlo `log2err`, direction and leading-edge IDs, symbols and protein labels. `leading_edges_all.tsv` expands every leading edge with rank and original protein labels. **`leading_edges_significant.tsv` has a header and zero rows because no category is significant.**

## Exploratory ORA

The notebook's original classification is reconstructed before HPA mapping: average the supplied normalized values by original gene symbol for EGFP, WT324 and WT328; scale each gene across these three condition means; run R k-means with k=3, seed 123 and default `nstart=1`; use clusters 2 and 3; require at least a twofold 324/328 effect and **0 < nominal p < 0.5**. Positive effects are 324-specific and negative effects 328-specific. This permissive historical threshold is reproduced, not presented as a validated differential-interaction call. No additional EGFP fold-change filter is added to these specific sets; the notebook contains a different EGFP-filtered candidate analysis that is not the source of these classifications.

The reproduced cluster centers confirm cluster 2 is highest in 324 and cluster 3 highest in 328. Centers, original-symbol cluster assignments and per-protein classifications are saved. At the mapped-gene level, any qualifying row assigns a gene to a specific set unless rows support both directions. One such conflicting gene is excluded from both hit sets but remains in the background.

- 324-specific: **171 genes**.
- 328-specific: **298 genes**.
- Shared/non-specific: **1,171 genes**. This label means not selected by the historical rule, not evidence of equal binding.
- Conflicting duplicate classifications: **1 gene**.
- Background: **1,641 unique, actually tested, successfully HPA-mapped genes**, including those with no localization annotation. No unmeasured HPA genes enter the background.
- One-sided Fisher exact test (`alternative='greater'`) for enrichment, with the hit set nested within the background. Counts form the standard 2x2 table: hit/in-location, hit/out-of-location, nonhit/in-location, nonhit/out-of-location. Odds ratios are R's conditional maximum-likelihood estimates.
- Same minimum size 5 and maximum <95% coverage rule; separate BH correction within each confidence-by-resolution-by-isoform family (eight analyses). Results report all counts, odds ratio, p-value, FDR and contributing IDs/symbols/proteins.

## Figures and interpretation

- Four NES dot plots: fine/collapsed crossed with primary/high confidence. Point size is measured gene count, color is BH FDR, and an asterisk-shaped point marks FDR <0.05 when present. No significant points occur here.
- Eight running enrichment plots: for each primary resolution, the two smallest nominal p-values per direction among the prespecified biologically relevant compartment terms listed in the analysis script. The selection table is saved. All are explicitly labeled nonsignificant; these plots illustrate rank patterns rather than significant findings.
- Two localization-by-analysis tile plots include primary GSEA, high-confidence GSEA and each primary ORA. Color is signed -log10(FDR), red for 324 and blue for 328. Asterisks denote FDR <0.05; grey means not tested. ORA signs indicate the tested hit set, not a signed depletion test. Values on this plot are evidence summaries, not directly comparable effect sizes across GSEA and ORA.

## Reproduction and validation

The requested exact-zero GSEA and annotated-background ORA follow-ups are saved separately in [sensitivity/README.md](sensitivity/README.md). Both preserve the no-significant-enrichment conclusion; the results above remain the primary analyses. Run `Rscript analysis/proteomics_hpa/run_sensitivities.R` after the primary workflow to reproduce the additional four GSEA and eight ORA runs, comparison tables and six figure pairs.

Run from the repository root with R and the packages in `sessionInfo.txt`:

```sh
Rscript analysis/proteomics_hpa/run_analysis.R
Rscript analysis/proteomics_hpa/validate_results.R
python3 analysis/proteomics_hpa/record_provenance.py
```

The scripts use frozen local inputs and require no network for reruns. Downloads can be reobtained from `source_manifest.tsv`, but current URLs may change release; compare checksums and version metadata before replacing frozen files.

`validation_checks.tsv` records independent reconstruction of median rank scores and weighted running-sum ES, a hypergeometric cross-check of Fisher p-values, BH-family checks, universe checks and figure pairing. The PNG plots and representative rendered PDF pages were visually inspected. No original data or pre-existing analysis outputs were overwritten.

## Biological limits

HPA supplies gene-level prior localization knowledge, mainly from immunofluorescence in selected human cell lines. It generally is not isoform-specific and does not demonstrate localization in this experiment's cell type or condition. Any future significant result should be phrased as enrichment for proteins **annotated by HPA to localize to** a compartment. Multi-label overlap, incomplete annotation, correlated proteins and the proteomics experiment's detectability/missingness affect interpretation. Gene-set randomization does not model sample-level biological replication. The exploratory specific sets use a permissive historical classifier. The present analysis supports no FDR-significant localization shift.

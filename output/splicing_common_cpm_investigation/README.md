# MSI2 landmark enrichment with common-CPM normalization

This analysis places CLIP and SMI coverage on a common primary-mapped-read CPM scale. Results from the original normalization are retained in `output/splicing_first_pass/` for methodological comparison. The workflow uses exported count tables and does not require access to FASTQ or BAM files.

## Common normalization

Every sample uses `raw_count / primary_mapped_reads * 1e6`. CLIP/SMI enrichment is calculated within the true construct/replicate pair with a primary pseudocount of 0.1 CPM and a sensitivity pseudocount of 0.5 CPM.

## Implicit zeros

The sparse bin table stores positive rows only. Missing sample/site/bin combinations are exact zeros. Every site-weighted denominator includes the full condition-independent retained-site universe.

## Replicate regions

Regional enrichment is calculated independently for the four matched CLIP/SMI pairs. Region means average the specified bins within each site and then weight retained sites equally. Bins are not treated as biological replicates.

## Side of splice

Donor and acceptor site tables contain paired exonic and intronic regional delta-delta values. Side difference is exonic minus intronic delta-delta, so negative values indicate a stronger 328-over-324 effect on the exonic side. Gene values are medians across sites assigned to the primary gene ID. Paired Wilcoxon tests, paired standardized mean differences, rank-biserial effects, and gene-cluster bootstrap intervals are reported.

## Coverage dependence

Coverage is condition-independent: whole-window overlap counts are converted to CPM per sample and summed across all eight samples. Balanced quartiles are defined within landmark type without using construct direction or CLIP enrichment; site ID breaks exact coverage ties deterministically. Coverage-stratified metaprofiles and primary regional summaries use site-weighted delta-delta. A pooled-read comparator within each quartile diagnoses the aggregate-weighting reversal.

## Limitations

Only two biological replicates exist per construct. Gene-cluster bootstrap intervals quantify variation across genes and do not create additional experimental replication. This analysis measures CLIP enrichment near landmarks, not PSI or transcript-level differential splicing. Annotated transcript 3' ends are not validated cleavage or PAS coordinates.

## Reproduction

```bash
Rscript analysis/investigate_landmark_normalization_and_coverage.R --input data/results --old-output output/splicing_first_pass --output output/splicing_common_cpm_investigation --bootstrap 1000 --pseudocount 0.1 --sensitivity-pseudocount 0.5 --seed 328324
```

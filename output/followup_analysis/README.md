# MSI2 follow-up analysis

The [HTML report](report.html) brings together coverage, normalization, spatial organization, motif, and external-RBP comparisons. Existing analyses remain in their original locations and are linked from the report. `analysis_comparison.tsv` summarizes the scientific comparisons; `new_analysis_configurations.tsv` records each calculation configuration. Audit tables distinguish the input state from the completed outputs.

The endogenous origin of the supplied MSI2 reference BED has not been confirmed by an experiment or accession. Analyses using that reference should be interpreted with this provenance limitation.

## Rebuild from the repository root

Prerequisites: R with data.table, ggplot2, readr, dplyr, tidyr, tibble, preprocessCore, ragg, GenomicRanges, GenomicAlignments, Matrix, jsonlite, rtracklayer, Biostrings, Rcpp; a C++ compiler; Python3 standard library. Existing input exports are required. Raw BAMs are not required for these follow-ups because every eligible aligned record's position/CIGAR and annotation membership is retained in the round-two exports.

1. `python3 output/followup_analysis/scripts/audit.py` initializes the pre-analysis audit.
2. `Rscript output/followup_analysis/scripts/raw_count_quantile.R --output output/followup_analysis/raw_count_quantile`
3. `Rscript output/followup_analysis/scripts/region_intron.R`
4. `Rscript output/followup_analysis/scripts/read_sensitivities.R`
5. `Rscript output/followup_analysis/scripts/plot_reads.R`
6. `python3 output/followup_analysis/scripts/prepare_regions.py`
7. `Rscript output/followup_analysis/scripts/site_spatial.R`
8. `Rscript output/followup_analysis/scripts/motif.R`
9. `Rscript output/followup_analysis/scripts/supplemental.R` followed by `Rscript output/followup_analysis/scripts/final_statistics.R` and `Rscript output/followup_analysis/scripts/UTR_motif_background.R`
10. `Rscript output/external_rbp_splice_comparison/analyze_external.R`
11. `Rscript output/external_rbp_splice_comparison/finish_external.R`
12. `Rscript output/followup_analysis/scripts/quality_checks.R` and `Rscript output/followup_analysis/scripts/final_presentation.R`
13. `python3 output/followup_analysis/scripts/verify_outputs.py`
14. `python3 output/followup_analysis/scripts/build_report.py`
15. `python3 output/followup_analysis/scripts/update_registry.py`

The scripts write only their designated follow-up/extension directories, except the final registry step, which appends entries and inserts an addendum. Original historical analyses are not replaced. Spatial sparse-matrix caches can be reused only with identical inputs; the input hash manifest documents this run. To analyze changed inputs, use a new output directory or invalidate those caches explicitly. When rebuilding, preserve historical outputs or write to a separate output directory.

## Inputs and provenance

- Original featureCounts matrix: `data/results/peak_analysis/peak_count_matrix.tsv`. Fixed nonzero-count universe: 16,046 peak IDs.
- Read-level source: data/results_round2/*/read_annotations.tsv.gz; eight genome libraries, 2,180,719 eligible aligned records. Each read's original record weight is 1 after upstream collapsing. Same-strand M/= /X block overlap, MAPQ>=10, primary unique alignments; N and D are not coverage.
- Genomic annotation: exact GENCODEv26 primary assembly release used by the project. Newly retrieved from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.primary_assembly.annotation.gtf.gz
- Sequence: UCSC hg38.2bit from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
- Chain: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
- PARpipe scripts: immutable commit 79b1aabfc985bdf750fa304bd1d217442da940b7, files scripts/Spatial.pl, scripts/Spatial.R, scripts/annotate.pl. Archived under data/reference_followup/. The repository tree records the commit. The half-integer 25-nt-margin spatial randomization is reproduced, not the obsolete pipeline.
- External RBP sources and conversion QC: output/external_rbp_splice_comparison/; untouched originals and metadata in data/external_rbp_controls/original/.

Read-overlap reconstruction deliberately uses every peak listed for each aligned record. Historical featureCounts did not enable -O, so ambiguous multi-peak reads can differ; see reconstructed_count_QC.tsv. Equal-depth sensitivity compares to its own reconstructed full-depth baseline, not an assertion of exact equality to the historical counts.

## Seeds and statistical units

| Procedure | Seed | Randomizations or sampling depth |
|---|---:|---|
| Exact read subsampling | 260921 | 50 independent draws per library; CLIP depth 93,659 and SMI depth 257,457 |
| Spatial site nulls | 260925 | 500 draws |
| External splice nulls | 260923 | 500 draws |
| Dinucleotide motif null | 260924 | 200 randomizations per 51-nt sequence |

Resampling uncertainty is technical, rather than biological. Feature- and peak-level tests can involve correlated genes or transcripts; their statistical unit does not represent independent biological replication.


## Outputs

`tables/` retains data for every new figure, including replicate and iteration data. `cache/` retains read-category memberships, per-read spatial matrices, resampling arrays, site windows and null distributions. `raw_count_quantile/` holds the exact new normalization variant. `figures/` contains matching PNG and PDF files. Original full-depth/scaled analyses and their established alternative normalizations remain linked from the report.

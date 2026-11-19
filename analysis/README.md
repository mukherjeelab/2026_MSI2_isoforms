# Downstream analysis scripts

Run scripts from the repository root. Detailed definitions, input requirements, and saved results are linked below.

| Scripts | Purpose and documentation |
|---|---|
| `visualize_peak_deseq.R`, `visualize_all_nonzero_peaks.R` | [Fixed-peak figures and quality control](../data/results/peak_analysis/figures/README.md) |
| `analyze_continuous_construct_preference.R` | [Continuous transcript-region comparison across 16,046 supported peaks](../data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/continuous_construct_preference/README.md) |
| `analyze_quantile_normalized_construct_preference.R` | [Quantile-normalized log2-CPM sensitivity analysis](../data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/continuous_construct_preference/quantile_normalized_log2cpm/README.md) |
| `analyze_construct_biased_peak_regions.R` | [Threshold-based and selected-peak transcript-region analyses](../data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/README.md) |
| `analyze_landmark_splicing_first_pass.R` | [Site-weighted landmark analysis using the original normalization](../output/splicing_first_pass/README.md) |
| `investigate_landmark_normalization_and_coverage.R` | [Common-CPM landmark profiles, regional contrasts, and coverage sensitivity](../output/splicing_common_cpm_investigation/README.md) |
| `diagnose_clip_smi_components.R` | [CLIP and SMI component decomposition](../output/splicing_clip_smi_diagnostic/README.md) |
| `compare_round2.R` | [Read-selection, peak-preference, and standardized feature-profile comparisons](../output/round2_comparisons/README.md) |
| `analyze_site_enrichment_distributions.R` | Site-level enrichment distributions; results in `output/splicing_site_enrichment_distributions/` |
| `proteomics_hpa/` | [Protein identifier mapping, localization enrichment, and sensitivity analyses](../output/proteomics_hpa/README.md) |
| `count_fastq_target.py`, `align_smjydd_*.py`, `summarize_smjydd_min500.py`, `check_smjydd_candidate_junctions.py` | [Construct-sequence checks and read classification](SMJYDD_ANALYSIS_NOTES.md) |
| `prepare_smjydd_igv.py`, `convert_smjydd_igv_bam.py`, `realign_l2_igv.py` | [Construct-reference alignment views and IGV resources](../outputs/SMJYDD_IGV/README.md) |

Additional scripts are stored with their result bundles under `output/followup_analysis/`, `output/external_rbp_splice_comparison/`, and `output/publication_figures/`. These locations are retained to preserve existing report links and provenance. See the [reproducibility guide](../docs/REPRODUCIBILITY.md) before rerunning analyses.

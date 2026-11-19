# Reproducibility guide

This guide describes the available computational workflows and the inputs needed to rerun them. Analysis-specific READMEs provide the detailed normalization, filtering, statistical, and interpretation requirements.

## Working environment

Run commands from the repository root, or open `2026_MSI2_isoforms.Rproj` in RStudio. The project contains R, Python, and Snakemake workflows; it does not provide a single locked environment for every analysis. Saved `sessionInfo.txt` files record R and package versions for individual runs. The CLIP workflow includes Conda environment specifications under [`code/clip_peak_analysis/envs/`](../code/clip_peak_analysis/envs/).

Install the dependencies of the workflow being reproduced. For example, the read/feature comparison uses `data.table`, `ggplot2`, and `ragg`; the continuous construct-preference analysis uses `dplyr`, `ggplot2`, `readr`, `ragg`, `tibble`, and `tidyr`. Broader follow-up analyses require additional Bioconductor packages and a C++ compiler, as listed in their [methods](../output/followup_analysis/README.md). Proteomics localization dependencies and input provenance are documented in the [HPA analysis](../output/proteomics_hpa/README.md).

## Workflow entry points

These are independent entry points, not a single ordered pipeline. Scripts can overwrite their saved outputs; preserve a checkout of the version being reproduced and inspect each script's input and output arguments first.

| Workflow | Command from the repository root | Details |
|---|---|---|
| Audit saved fixed-peak counts | `Rscript code/clip_peak_analysis/scripts/audit_results.R` | [Checks and inputs](../code/clip_peak_analysis/README.md) |
| Continuous preference, all 16,046 supported peaks | `Rscript analysis/analyze_continuous_construct_preference.R` | [Methods](../data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/continuous_construct_preference/README.md) |
| Quantile-normalized log2-CPM sensitivity | `Rscript analysis/analyze_quantile_normalized_construct_preference.R` | [Methods](../data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/continuous_construct_preference/quantile_normalized_log2cpm/README.md) |
| Selected-peak publication figure | `Rscript output/publication_figures/plot_construct_preference.R` | [Selection and tests](../output/publication_figures/construct_preference_methods.md) |
| Read-selection and feature-profile comparison | `Rscript analysis/compare_round2.R` | [Methods](../output/round2_comparisons/README.md) |
| Common-CPM landmark profiles | `Rscript analysis/investigate_landmark_normalization_and_coverage.R` | [Full command and parameters](../output/splicing_common_cpm_investigation/README.md) |
| CLIP/SMI component decomposition | `Rscript analysis/diagnose_clip_smi_components.R` | [Full command and prerequisite outputs](../output/splicing_clip_smi_diagnostic/README.md) |
| Proteomics localization | `Rscript analysis/proteomics_hpa/run_analysis.R` | [Workflow and follow-up commands](../output/proteomics_hpa/README.md) |

For the ordered coverage, spatial, motif, and external-RBP workflows, follow the [follow-up reproduction instructions](../output/followup_analysis/README.md). For upstream fixed-peak counting and paired DESeq2, follow the [Snakemake instructions](../code/clip_peak_analysis/README.md), including configuration of the external BAM location and SLURM profile.

## Analysis definitions

Several saved analyses address related questions with different inputs. Their peak universes and normalization choices must be kept distinct:

- **43,794 merged fixed regions:** the complete featureCounts matrix after source-peak filtering and same-strand merging.
- **16,046 supported peaks:** the original fixed-peak universe with a nonzero raw count in at least one of eight libraries; used for the continuous transcript-region analysis.
- **16,711 supported peaks:** the separately extracted round-two full-genome universe, with a different annotation procedure.
- **9,369 selected peaks:** the source set for the figure in `output/publication_figures/`, after the source analysis excludes peaks depleted in both constructs.

Construct preference is generally reported as MSI2-324 enrichment minus MSI2-328 enrichment, but enrichment calculations differ among workflows. Consult the corresponding methods for replicate averaging, CPM denominators, pseudocounts, and transcript-region assignment. Within-construct significance alone is not a test of a between-construct difference. Peak- and gene-level analyses retain their documented statistical units and limitations.

## Inputs stored outside Git

| Input or artifact | Availability |
|---|---|
| Raw CLIP/input BAMs | External to this repository; configure the accessible BAM location before running featureCounts. |
| `data/reference_followup/hg38.2bit` | Download from [UCSC hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit) to this path for sequence-dependent follow-ups. |
| `output/analysis_comparison_registry.zip` and `output/analysis_comparison_registry_with_followups.zip` | Generated distribution bundles; the unpacked registry is versioned in `output/analysis_comparison_registry/`. |
| `tmp/` and local build environments | Not versioned. Install required software separately; construct-sequence/IGV scripts may reference project-local installations. |

Source URLs, input hashes, and reference versions are recorded in the relevant analysis manifests and methods. Before reusing a spatial cache, verify that the inputs match its manifest; changed inputs require cache invalidation or a separate output directory.

## Historical records

The portable [analysis comparison registry](../output/analysis_comparison_registry/README.md) is a saved evidence snapshot. Its copied sources, manifests, and build metadata describe the state at the time it was assembled. Later documentation edits do not update that historical snapshot. Earlier notebooks, normalization variants, and preliminary summaries remain available for provenance and methodological comparison.

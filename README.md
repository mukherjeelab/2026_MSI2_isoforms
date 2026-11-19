# Isoform-Specific Localization Diversifies Human MSI2 Function

Data, analysis code, and supporting results for the manuscript by Kathryn Walters, Nadine Koertel, Amber Baldwin, Marcin Sajek, and Neelanjan Mukherjee.

**Manuscript status:** submitted for peer review.

**Preprint:** [bioRxiv, version 2](https://www.biorxiv.org/content/10.64898/2026.05.16.725683v2) · **DOI:** [10.64898/2026.05.16.725683](https://doi.org/10.64898/2026.05.16.725683)

This repository brings together reporter-assay analyses, proteomics, and RNA-binding analyses of MSI2 isoforms, with particular emphasis on MSI2-324 and MSI2-328. Saved data tables, figure source data, workflow parameters, and supporting methodological comparisons accompany the analysis code.

## Contents

| Directory | Contents |
|---|---|
| [`code/`](code/README.md) | Reporter-assay and proteomics notebooks; fixed-peak CLIP counting and DESeq2 workflow |
| [`analysis/`](analysis/README.md) | Scripts for CLIP enrichment, transcript-region comparisons, landmark profiles, proteomics localization, and construct-sequence checks |
| [`data/`](data/) | Experimental input tables, reference annotations, and processed sequencing results |
| [`figures/`](figures/) | Figure files associated with the analysis notebooks |
| [`output/`](output/) | Analysis reports, statistical tables, figure source data, and methodological comparisons |
| [`outputs/`](outputs/) | Construct-sequence validation tables, alignments, and IGV resources |

## Analysis guide

| Analysis | Code and methods | Results |
|---|---|---|
| Dual-luciferase reporter assays | [`DualLuciferase.Rmd`](code/DualLuciferase.Rmd) | [`figures/`](figures/) |
| Protein interactions and pathway enrichment | [`Proetomics.Rmd`](code/Proetomics.Rmd) | [`figures/`](figures/) |
| Fixed-peak CLIP enrichment | [Workflow and experimental design](code/clip_peak_analysis/README.md) | [Counts, DESeq2 results, and QC](data/results/peak_analysis/README.md) |
| Construct preference by transcript region | [Continuous analysis](data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/continuous_construct_preference/README.md) | [Figure index and source tables](data/results/peak_analysis/figures/README.md) |
| Selected-peak construct-preference figure | [Selection, annotation, and statistical methods](output/publication_figures/construct_preference_methods.md) | [`publication_figures/`](output/publication_figures/) |
| Landmark enrichment and normalization | [Common-CPM methods](output/splicing_common_cpm_investigation/README.md) | [CLIP and SMI component analysis](output/splicing_clip_smi_diagnostic/README.md) |
| Read selection and standardized coverage profiles | [Methods and findings](output/round2_comparisons/README.md) | [HTML report](output/round2_comparisons/report.html) |
| Coverage, spatial organization, and motif controls | [Methods and reproduction](output/followup_analysis/README.md) | [HTML report](output/followup_analysis/report.html) |
| Proteomics subcellular-localization analysis | [HPA mapping and enrichment methods](output/proteomics_hpa/README.md) | [Results summary](output/proteomics_hpa/RESULTS.md) |
| Construct-sequence validation | [Analysis notes](analysis/SMJYDD_ANALYSIS_NOTES.md) | [Sample-specific alignments](outputs/SMJYDD_sample_specific_alignment/README.md), [IGV guide](outputs/SMJYDD_IGV/README.md) |

The guide identifies analysis resources by topic. Supporting sensitivity analyses and historical comparisons are retained alongside the manuscript-related work; directory names such as `first_pass`, `round2`, and `followup_analysis` preserve existing paths and provenance. The [analysis registry](output/analysis_comparison_registry/README.md) documents relationships among saved configurations. A directory name or inclusion in this repository does not establish a correspondence to a particular submitted figure panel.

## Reproducibility and data availability

See the [reproducibility guide](docs/REPRODUCIBILITY.md) for prerequisites, workflow entry points, reference downloads, and the distinction between inspecting saved results and rerunning upstream processing. Run commands from the repository root. Download the repository and open HTML reports locally to view their linked figures.

Processed CLIP counts, sample annotations, result tables, and many reference inputs are versioned here. Raw CLIP BAMs are stored outside the repository. The downloaded hg38 genome sequence, generated ZIP bundles, local software environments, and temporary build files are excluded from Git. Consequently, a checkout supports inspection of the saved analyses but requires additional inputs and software for some reruns.

## Citation

Walters K, Koertel N, Baldwin A, Sajek M, Mukherjee N. *Isoform-Specific Localization Diversifies Human MSI2 Function*. bioRxiv. 2026. [doi:10.64898/2026.05.16.725683](https://doi.org/10.64898/2026.05.16.725683). See the [submitted preprint version](https://www.biorxiv.org/content/10.64898/2026.05.16.725683v2).

When referring to a repository analysis, also record the Git commit and the relevant script or result path so that the exact computational version can be identified.

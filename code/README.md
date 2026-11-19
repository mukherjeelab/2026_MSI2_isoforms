# Analysis notebooks and CLIP workflow

| Resource | Purpose |
|---|---|
| [`DualLuciferase.Rmd`](DualLuciferase.Rmd) | Reporter-assay processing, statistical comparisons, and figures |
| [`Proetomics.Rmd`](Proetomics.Rmd) | Protein-interaction comparisons, clustering, and pathway enrichment; original filename retained |
| [`SpliceAnalysis.Rmd`](SpliceAnalysis.Rmd) | Splice-junction and isoform-related analyses |
| [`clip_peak_analysis/`](clip_peak_analysis/README.md) | Fixed-peak preparation, featureCounts, paired DESeq2, sample metadata, environments, and result audit |
| [`IsoformCLIP.Rmd`](IsoformCLIP.Rmd) | Historical unfinished notebook, retained for provenance; use the scripted CLIP workflows for reproducible downstream analyses |

Open notebooks from the repository root and review their input paths and package requirements before knitting. Some paths refer to the original computational environment. Rendered notebook files and associated figure directories are retained where available.

The [`analysis/` guide](../analysis/README.md) lists standalone downstream scripts. See the [reproducibility guide](../docs/REPRODUCIBILITY.md) for input availability and environment requirements.

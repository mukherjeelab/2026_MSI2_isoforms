# MSI2 analysis comparison registry

Open `report.html` in a browser. Everything it links is inside this directory; no network or repository checkout is needed. The package contains 104 completed parameter-specific analysis variants across 9 families, one proposed interaction analysis, and 101 figure records (PNG/PDF pairs counted once).

Repository commit: `ace0f5236b33493d7ea61a664cafd93e69e79449`. Analysis files include uncommitted/untracked work; the copied source files and `source_manifest.tsv` hashes define the actual evidence snapshot. Dates are output mtimes, not guaranteed execution dates.

## Contents

- `analysis_registry.tsv`: one row per saved methodological configuration, plus the proposed interaction model. Both original repository paths and portable copied paths are provided.
- `figure_registry.tsv`: all relevant existing PNG/PDF pairs, analysis links, caveats and gallery captions.
- `analysis_relationships.tsv`: parent, exact changes and held-constant choices.
- `quantitative_sensitivity.tsv`: numbers extracted from saved tables, with source paths.
- `grant_figure_candidates.tsv`: five presentation candidates with sensitivity/qualifying companions.
- `sources/`: copied scripts, results, figures, READMEs, configs and presentation sources. Historical analyses were not modified.
- `source_manifest.tsv`: SHA-256 checksums for copied files.
- `repository_inventory.tsv`: whole-tree discovery inventory and scope decisions.
- `presentation_text_inventory.tsv`: extracted slide/notes text; slide-only numbers are not automatically treated as verified.
- `validation.json`: row-count, provenance-link, image-format and portability checks.
- `build_metadata.json`: commit, working-tree status and count definitions.

## Counting convention

Landmarks, replicates and plot styles are grouped within their common calculation framework. Changing a saved selection, pseudocount, normalization, weighting or statistical unit produces a separate row. Round two contributes 60 completed configurations: 30 peak-preference and 30 depth/contrast configurations. Its fixed-nt and scaled-coordinate panels are grouped, with both explicitly recorded. Table-only sensitivities are counted; hypothetical runs are not. This is an explicit audit convention, not a claim that all rows are independent experiments.

## Rebuild

From the original repository root:

`python3 output/analysis_comparison_registry/build_registry.py`

From this portable directory, rebuild into another directory using only bundled evidence:

`python3 build_registry.py --source-root sources --output /tmp/msi2_registry_rebuild`

The portable rebuild cannot retrieve the original git history; retain this build's metadata. It requires only Python's standard library. It reads and summarizes existing outputs and copies evidence; it does not execute analysis scripts, normalize counts, or fit tests.

## Important provenance gaps

The current script no longer writes the historical jointly-positive outputs, although their saved values remain available. The historical transcript_region_enrichment.tsv also differs from current documentation: seven categories and counted-background comparisons, retained separately in the registry. Original/round-two extraction implementations and BAMs are absent. IsoformCLIP.Rmd is an unfinished notebook, not an additional completed analysis. Gene-biotype claims in the old deck lack a separately saved annotation table. Root Rplots.pdf is unlabelled. The construct-by-assay interaction remains proposed, not run. Site-distribution and gene-stratified proximity work is completed, with sparse-stratum estimation limitations.

Raw input BAMs and the largest sparse extraction tables are not included; all figures and result tables needed to review this report are included. This is an audit package, not a complete upstream reprocessing environment.

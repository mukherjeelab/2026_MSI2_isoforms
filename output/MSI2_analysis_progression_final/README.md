# MSI2 analysis progression

The presentation contains 35 main slides and 24 appendix slides. It follows the investigation from discrete peak classes through continuous preference, normalization, downsampling, spatial diagnostics, motif validation and external RBP controls. The appendix groups 123 completed configurations and one proposed model, and indexes all 129 registered figures.

## Files

- `MSI2_analysis_progression_final.pptx`: editable text and native tables; scientific plots embedded as images.
- `MSI2_analysis_progression_final.pdf`: matching 59-page slide PDF, generated from the final PowerPoint renderings. PDF pages are rasterized; use the PPTX for editing and searchable slide text.
- `slide_source_registry.tsv`: slide-level analysis IDs, figures, result tables, scripts and principal claims.
- `slide_content.json`: editable presentation content and full speaker-note provenance.
- `quantitative_evidence.json`: snapshots of directly cited small source tables.
- `source_sha256.tsv`: hashes of cited source files.
- `prepare_content.py`, `build_presentation.mjs`, `render_saved_plots.R`, `export_pdf.py`: reproducible authoring scripts.

## Revision

Run from the repository root. `prepare_content.py` reads the existing analysis outputs and refreshes content and provenance. `render_saved_plots.R` makes two presentation-only renderings from saved TSV values: full-genome expanded donor/acceptor profiles and a landscape external-RBP comparison. It does not rerun statistical tests or randomizations. All other embedded plots are original saved figures.

The PowerPoint uses the bundled Codex `@oai/artifact-tool` JavaScript runtime. `build_presentation.mjs` includes the runtime and presentation-skill paths used for this export; adjust those paths on a different computer. Run with bundled Node, with `RUNTIME_NODE_MODULES` set to its node_modules directory. It writes private validation receipts and slide renders, then updates the delivery PPTX. Run `export_pdf.py` with Python containing reportlab and pypdf after rendering completes. The generated content JSON can also be edited directly before rebuilding.

## Evidence boundaries

The supplied MSI2 omniCLIP reference is provisionally endogenous based on the user's recollection; no accession or publication was confirmed. This uncertainty is retained in the deck. The formal construct-by-assay interaction remains proposed. Read coverage, equal-site mean log ratios and called-peak localization are different measurements. Technical subsampling and permutation envelopes do not replace biological replication.

All 59 slides were rendered and visually reviewed. Crowded appendix tables were resized before final delivery. Numerical summaries were checked against saved result tables, including the normalization comparison, read fractions, component decomposition, randomization results and conversion QC. No new biological analysis was performed for this presentation.

# CLIP and SMI component diagnostic

This additive analysis leaves the primary common-CPM CLIP/SMI analysis unchanged.

CLIP-only panels show observed CLIP CPM per retained landmark. SMI-only panels show the matched RNA/background CPM distribution. Neither is called enrichment alone.

The exact component decomposition uses site-level log terms with the same 0.1 CPM pseudocount as the primary analysis. The CLIP contribution plus the signed SMI contribution equals the unchanged delta-delta profile to numerical precision.

Every sample/site/bin absent from the sparse count table is treated as zero, and every mean uses the full retained-site denominator. The retained filter is one or more raw whole-window overlaps across all eight samples combined.

## Reproduction

```bash
Rscript analysis/diagnose_clip_smi_components.R --input data/results --primary-output output/splicing_common_cpm_investigation --output output/splicing_clip_smi_diagnostic --pseudocount 0.1
```

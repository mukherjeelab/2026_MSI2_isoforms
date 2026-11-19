# MSI2 landmark enrichment: first-pass findings

## Primary site-weighted result

Across the condition-independent retained-site universe, MSI2-328 has higher CLIP-versus-SMI enrichment than MSI2-324 at every supplied 10-nt bin for all three landmark types. Delta-delta is defined as 324 enrichment minus 328 enrichment, so every primary estimate is negative.

- 5' splice donors: the mean delta-delta across the -250 to +250 nt profile is -0.100. The strongest difference is -0.263 at the bin centered at -5 nt, immediately on the exonic side of the donor boundary. The difference becomes substantially smaller on the intronic side but remains negative.
- 3' splice acceptors: the mean delta-delta across the -250 to +250 nt profile is -0.108. The strongest difference is -0.272 at the bin centered at +15 nt, immediately on the exonic side of the acceptor boundary. The difference decreases with distance into the exon but remains negative.
- Annotated transcript 3' ends: the mean delta-delta across the -50 to +50 nt profile is -0.357. The profile ranges from -0.417 at -45 nt to -0.288 at +45 nt. This is a broad difference across the annotated transcript-end window, not evidence about validated cleavage or polyadenylation sites.

The gene-bootstrap 95% interval excludes zero in the negative direction at all 50 donor bins, all 50 acceptor bins, and all 10 annotated transcript-end bins. These intervals quantify variation across annotated gene clusters. They do not replace experimental-replicate inference.

## Replicate consistency

The replicate-level site-weighted profile shapes are highly consistent. Pearson correlations between replicate A and replicate B range from 0.994 to 1.000 across construct and landmark combinations. Spearman correlations range from 0.952 to 1.000. The transcript-end correlations use only 10 bins and should be interpreted accordingly.

## Pseudocount sensitivity

Repeating the site-weighted analysis with a pseudocount of 0.5 instead of 0.1 preserves the delta-delta sign at every bin. Profile correlations between the two analyses are 0.997 to 1.000. Effect magnitudes shrink: the mean delta-delta changes from -0.100 to -0.036 for donors, -0.108 to -0.040 for acceptors, and -0.357 to -0.149 for annotated transcript 3' ends. The direction and profile shape are robust, but the magnitude is sensitive to the pseudocount because the site-level data are sparse.

## Site-weighted and read-weighted results differ

The secondary read-weighted analysis does not reproduce the primary result uniformly.

- Donor aggregate delta-delta is positive across the full window and averages +0.210.
- Acceptor aggregate delta-delta is mixed across position and averages +0.111.
- Annotated transcript-end aggregate delta-delta changes from positive upstream to negative downstream and averages +0.141 across the full window.

This contrast means the most heavily covered sites and the typical retained site give different summaries. The site-weighted analysis remains the primary biological result because each landmark contributes equally. The aggregate analysis should be interpreted only as the distribution of reads across landmarks.

## Analysis universe

- 101,162 donor sites across 13,128 primary gene IDs
- 100,801 acceptor sites across 12,992 primary gene IDs
- 35,769 annotated transcript 3' ends across 10,172 primary gene IDs

All absent sparse sample/site/bin records were treated as exact zeros. No CLIP-specific or direction-dependent site filter was applied.

## Interpretation boundary

This first pass supports a reproducible difference in CLIP-versus-SMI landmark profiles between MSI2-324 and MSI2-328. It does not estimate PSI, demonstrate differential splicing, or identify validated cleavage or PAS positions. Formal regional and side-of-boundary tests remain intentionally outside this first-pass deliverable.

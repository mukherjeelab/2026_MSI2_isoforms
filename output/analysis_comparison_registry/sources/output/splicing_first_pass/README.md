# MSI2 landmark enrichment first-pass analysis

This directory contains the requested local first-pass analysis of the completed BAM-derived landmark bundle. No FASTQ or BAM file was accessed.

## Inputs

- `data/results/sample_metadata.tsv`
- `data/results/landmark_sites.tsv.gz`
- `data/results/landmark_bin_counts.tsv.gz`
- `data/results/landmark_metaprofile_raw.tsv.gz` (reconciliation only)
- `data/results/landmark_extraction_QC.tsv`

## Implicit zeros

`landmark_bin_counts.tsv.gz` contains nonzero rows only. Every absent sample/site/bin combination was treated as an exact zero. Site-weighted means use the full retained-site denominator, including sites with zero coverage in a given sample/bin.

## Normalization

Counts were divided by the existing project-specific full-library factor copied into `sample_metadata.tsv`. Factors were originally computed within each four-library construct comparison as primary mapped reads divided by the comparison-specific geometric mean (324: 230995.596949; 328: 214594.397648). Unique-mapping percentages were not used for normalization.

## Matched enrichment and delta-delta

CLIP and SMI were matched within construct and biological replicate using `pair_id`. Replicate-level enrichment is `log2((normalized CLIP + 0.1) / (normalized SMI + 0.1))`. The pseudocount is applied after normalization. Construct profiles are the mean of the two replicate-level profiles. Delta-delta is `enrichment_324 - enrichment_328`; positive values favor 324 and negative values favor 328.

## Site weighting and uncertainty

The primary profiles first calculate enrichment for each site and bin, then average over all retained sites. The 95% intervals use 1,000 cluster-bootstrap iterations over primary gene IDs. All sites are assigned to the first listed gene ID for bootstrap clustering; `site_gene_counts.tsv` reports multi-gene annotation counts.

## Read-weighted aggregate analysis

Secondary aggregate profiles sum normalized counts within each sample/bin, add the documented pseudocount to those summed values, then calculate matched CLIP/SMI enrichment. These profiles answer where reads are concentrated and are not the primary site-weighted result.

## Pseudocount sensitivity

The primary pseudocount is 0.1 normalized counts. The site-weighted analysis was repeated with 0.5; results are in `results/qc/pseudocount_sensitivity.tsv` and its summary.

## Orientation and terminology

All positions are transcript-oriented. Negative values are upstream and positive values are downstream. Donor plots label EXON | INTRON, acceptor plots label INTRON | EXON, and transcript-end plots label TRANSCRIPT | DOWNSTREAM. Transcript ends are described only as annotated transcript 3' ends, not cleavage or polyadenylation sites.

## Site filtering

The primary universe is the condition-independent `retained_for_bin_counts == yes` set from extraction. No direction-dependent or CLIP-specific filter was added.

## Limitations

Only two biological replicates are available per construct/assay. Gene-bootstrap intervals describe uncertainty across annotated genes, not experimental-replicate uncertainty. Mean site-level log ratios are strongly shrunk toward zero by site/bin combinations with zero in both matched libraries, as intended by the full retained-site universe. The analysis addresses CLIP enrichment around landmarks; it does not estimate PSI or establish differential splicing. Annotated transcript 3' ends are not validated cleavage or PAS coordinates.

## Reproduction

Run:

```bash
Rscript analysis/analyze_landmark_splicing_first_pass.R --input data/results --output output/splicing_first_pass --bootstrap 1000 --pseudocount 0.1 --sensitivity-pseudocount 0.5 --seed 328324
```

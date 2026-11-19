# MSI2 CLIP RNA-processing landmark extraction

This workflow performs **data extraction and compilation only**. It does not
normalize landmark counts, pool replicates, calculate 324-versus-328
delta-delta statistics, or generate final figures.

It is intentionally separate from the existing peak-DESeq workflow. The bundle
uses absolute Bodhi paths for the existing BAM, STAR-log, and peak-normalization
outputs, so `landmark-extraction/` can be placed anywhere on Bodhi.

## Required configuration

Edit `config/config.yaml` and replace:

```yaml
gtf: /absolute/path/to/annotation.gtf
```

The GTF must use the same genome build and chromosome names as the STAR BAMs.
The eight BAM paths and complete paired sample design are in
`config/samples.tsv`. The workflow validates that every `pair_id` has exactly
one CLIP and one SMI library with matching construct and replicate labels.

The workflow does not infer validated cleavage or PAS sites from transcript
ends. Optional validated annotations can be supplied as one-nucleotide BED6+
files:

```text
chrom  start0  end0  site_name  score  strand  gene_id  gene_name  transcript_ids
```

Columns 7-9 are optional, and `transcript_ids` is semicolon separated. Every
interval must have `end0-start0 = 1`.

- `cleavage_sites_bed`: validated cleavage/polyadenylation coordinates.
- `pas_sites_bed`: validated polyadenylation-signal motif coordinates.

Leave either value blank when the annotation is unavailable. No PAS motif
coordinates are invented. GTF-derived transcript ends are instead labeled
`annotated_transcript_3prime_end`.

## Transfer and run on SLURM

Copy the supplied archive from the local computer:

```bash
scp landmark-extraction.tar.gz \
  <netid>@<bodhi-login-host>:/beevol/home/bartholk/2026_MSI2_CLIP_iso/
```

Then unpack it on Bodhi:

```bash
cd /beevol/home/bartholk/2026_MSI2_CLIP_iso
tar -xzf landmark-extraction.tar.gz
cd landmark-extraction
```

If the preprocessing Snakemake environment already exists, activate it.
Otherwise create the included driver environment once:

```bash
conda env create -f envs/snakemake.yaml
conda activate landmark-snakemake
```

Before submission, confirm one BAM and its index and replace the placeholder
`gtf` in `config/config.yaml`:

```bash
ls /beevol/home/bartholk/2026_MSI2_CLIP_iso/clip-preprocessing/results/alignment/NM2026_0081.mapped.sorted.bam{,.bai}
grep -E 'sjdbGTFfile|genomeFastaFiles' \
  /beevol/home/bartholk/2026_MSI2_CLIP_iso/star_index/genomeParameters.txt
```

The `sjdbGTFfile` line can help identify the GTF used to build the STAR index,
but still verify that the file exists and matches the BAM genome build.

Check the workflow graph without submitting jobs, then run it:

```bash
snakemake --profile profiles/slurm --dry-run
snakemake --profile profiles/slurm
```

Equivalently, the included wrapper accepts normal Snakemake arguments:

```bash
bash run_on_slurm.sh --dry-run
bash run_on_slurm.sh
```

Snakemake creates the rule environment from `envs/extraction.yaml`, containing
Python and Samtools.

## Exact junction support

Junctions are identified directly from `N` operations in each primary BAM
alignment's CIGAR string. They are not estimated from STAR `Log.final.out` and
are not obtained by summing `SJ.out.tab` counts.

Default read eligibility is:

- primary alignment: SAM flags 0x4, 0x100, and 0x800 excluded;
- mapping quality at least 10;
- `NH` absent or `NH <= 1` because the genome STAR workflow requested unique
  mapping;
- read strand matching annotation strand, consistent with the current R1 CLIP
  peak-counting configuration.

`junction_support.tsv.gz` contains every supplied-GTF junction for every sample,
including zero-support rows, plus observed unannotated junctions. An alignment
spanning multiple junctions contributes once to each junction but only once to
the sample QC count `alignments_with_CIGAR_N`.

Two support quantities are reported:

- `unique_primary_alignments_spanning`: eligible primary SAM records.
- `unique_fragments_read_names_spanning`: distinct QNAMEs supporting the
  junction. These are equal for the present single-end BAMs; QNAME counts avoid
  double-counting mates if paired-end BAMs are supplied later.

The local denominator is explicitly defined as eligible primary alignment
reference spans overlapping:

```text
[intron_start - junction_local_flank,
 intron_end + junction_local_flank)
```

The default flank is 250 nt. This is a local genomic-locus denominator, not a
claim that all overlapping reads are assignable to one isoform. Junctions
shared by multiple transcripts retain all transcript IDs.

Junction coordinates `intron_start` and `intron_end` are zero-based, half-open
and exactly match the skipped-reference interval encoded by CIGAR `N`.
One-based inclusive helper columns are also included.

## Landmark definitions and orientation

Landmarks are deduplicated by chromosome, genomic position, strand, and site
type. Associated gene IDs, gene names, transcript IDs, and annotation sources
are retained as semicolon-separated values.

`genomic_position` is the one-based coordinate of the landmark nucleotide:

- `splice_donor_5prime`: first intronic nucleotide in transcript orientation;
- `splice_acceptor_3prime`: last intronic nucleotide in transcript orientation;
- `annotated_transcript_3prime_end`: last annotated transcribed nucleotide;
- external cleavage/PAS annotations: the nucleotide in the one-base BED entry.

Relative coordinates always follow transcript 5-prime to 3-prime orientation.
Minus-strand bins are reversed so negative positions are upstream and positive
positions are downstream.

The half-open primary bins are:

- splice donor and acceptor: `[-250, 250)` in 10-nt bins;
- cleavage, PAS, and annotated transcript 3-prime ends: `[-50, 50)` in 10-nt
  bins.

Therefore each splice site has 50 possible bins and each 3-prime landmark has
10 possible bins.

`raw_overlap_count` is the number of eligible deduplicated primary alignment
records with at least one aligned `M`, `=`, or `X` base overlapping the bin.
Skipped `N` bases do not count as coverage. No normalization is applied.

## Condition-independent coverage filter

Full site-by-bin extraction across every annotation model can be unnecessarily
large. The default filter retains a site when its complete landmark window has
at least one eligible alignment across all eight samples combined:

```text
sum(raw window-overlap counts across every sample) >= 1
```

This filter does not use construct, assay, replicate, fold change, or direction
of enrichment. Change `minimum_total_window_overlap_count` only before running
the extraction, and document the new value.

The workflow never silently discards the unfiltered information:

- `landmark_sites.tsv.gz` includes every valid deduplicated site and records
  whether it was retained.
- `landmark_site_summary.tsv.gz` gives unfiltered per-sample/per-site window
  counts.
- `landmark_extraction_QC.tsv` reports exactly how many sites were annotated,
  on the BAM reference, covered, retained, and removed for each site type.
- `landmark_metaprofile_raw.tsv.gz` gives a raw, unnormalized aggregate for
  every individual sample and site type.

## Sparse bin-count representation

`landmark_bin_counts.tsv.gz` is a compact compressed TSV rather than Parquet,
avoiding a PyArrow dependency on the cluster. With the default
`sparse_nonzero_bins: true`, zero-count sample/site/bin rows are omitted.

These omissions are exact implicit zeros. They can be reconstructed from:

1. sites where `retained_for_bin_counts == yes` in `landmark_sites.tsv.gz`;
2. the site-type bin ranges documented above; and
3. all samples in `sample_metadata.tsv`.

No nonzero site/bin observation is aggregated or removed. Set
`sparse_nonzero_bins: false` to emit the full dense table if storage permits.

For the QC metrics `genes_represented` and
`sites_overlapping_genes_represented_in_BAMs`, a represented gene means that at
least one eligible alignment overlaps the union of its annotated GTF exons in
at least one sample. This is computed independently of landmark-window
coverage and condition labels.

## Outputs to transfer

All transferable products are written under `results/`:

```text
sample_metadata.tsv
landmark_sites.tsv.gz
landmark_site_summary.tsv.gz
landmark_bin_counts.tsv.gz
landmark_metaprofile_raw.tsv.gz
junction_support.tsv.gz
landmark_extraction_QC.tsv
README.md
```

### `sample_metadata.tsv`

Preserves sample, construct, CLIP/SMI assay, biological replicate, matched pair,
comparison, and BAM path. It also reports BAM primary mapped counts, MAPQ-filtered
counts, paired-layout detection, available STAR unique-mapping metrics, and the
existing project's full-library normalization quantities.

The workflow copies existing normalization values but does not replace or apply
them. The current peak analysis uses `full_library` normalization based on
primary mapped alignments and stores exact DESeq2 size factors under
`peak-deseq-analysis/results/peak_analysis/deseq2/*/size_factors.tsv`. If those
results live somewhere else on Bodhi, update the two `existing_*` paths in
`config/config.yaml`; missing optional normalization files are recorded as
`not_available`/`NA`, not fabricated.

### `landmark_sites.tsv.gz`

One row per deduplicated landmark, including the requested fields:

```text
site_id, site_type, chromosome, genomic_position, strand,
gene_id, gene_name, transcript_ids, annotation_source
```

Additional coordinate, window, coverage-filter, and retention fields make the
bin table completely reproducible.

### `landmark_bin_counts.tsv.gz`

One nonzero row per retained sample/site/bin by default:

```text
sample_id, site_id, site_type,
relative_bin_start, relative_bin_end, raw_overlap_count
```

Replicates are never pooled.

### `landmark_metaprofile_raw.tsv.gz`

An unnormalized per-sample/per-site-type aggregate. It includes the retained
site denominator, summed raw overlap count, and mean raw overlap count per
retained site for every relative bin. It is provided for QC/convenience only.

## Interpretation limitation

A junction can support multiple annotated isoforms. A read supporting a shared
junction is evidence for every compatible transcript, not proof of one specific
isoform. Use transcript-unique junctions or a transcript-level model when an
isoform-specific claim is required.

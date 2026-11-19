# Original-library read and feature analysis

This independent workflow starts from the eight original 324/328 CLIP and SMI
**genome-aligned BAMs before any peak restriction**. It prepares per-read
annotations and feature profiles, preserving sample, assay, construct, matched
pair, and biological replicate. It does not realign reads, call new peaks, pool
replicates, run differential tests, or subtract SMI from CLIP.

The source BAMs already reflect repeat subtraction, UMI+insert collapsing, and
unique genome mapping from the existing preprocessing pipeline. “All reads” here
means all eligible alignments in those BAMs, not the original raw FASTQ reads.
Each collapsed BAM record counts once; the original abundance encoded in the
FASTA/read name is not used as a weight.

## Analysis matrix

Every sample gets all ten variants below. Motif filtering occurs **before**
computing the equal-depth targets.

| Mode | All eligible reads | UAG-containing eligible reads | Depth definition |
|---|---|---|---|
| `genome_full` | yes | yes | No subsampling, no peak restriction |
| `genome_equal` | yes | yes | Exact minimum genome count across all eight libraries |
| `peaks_full` | yes | yes | Every read overlapping the supplied peaks |
| `peaks_after_genome_equal` | yes | yes | Peak-overlapping subset of the same `genome_equal` draw |
| `peaks_equal` | yes | yes | Exact minimum peak-overlapping read count across all eight libraries |

Variant names are prefixed `all__` or `uag__`. UAG-only targets are determined
from the UAG-only library sizes, independently of the all-read targets. Those
independent random draws are not necessarily nested; every selected UAG read is,
however, a member of `all__genome_full`.

`peaks_after_genome_equal` does not promise equal peak counts. It exposes the
remaining difference in the fraction of the genome library overlapping peaks.
`peaks_equal` deliberately conditions on that fraction and can remove biological
differences in overall peak enrichment. Treat these as sensitivity analyses;
subsampling does not establish that every source of library bias is removed.

Exact sequential random sampling selects k of n eligible records without
replacement, in constant memory. Seed 1729 is saved in QC; reproducibility
requires unchanged BAM record order and the same Python environment. To assess
sampling variability, rerun with another seed **and a different output_dir**:

```bash
bash run_on_slurm.sh --config seed=2718 output_dir=results_seed2718
```

Do not treat different random draws as biological replicates. The source BAMs
are streamed twice (census, extraction); selected read identities and variant
membership are retained without writing ten extra BAMs per sample.

## Read eligibility and motif

- Single-end primary mapped records, excluding flags 0x4, 0x100, and 0x800.
- MAPQ >=10; NH absent or <=1, consistent with the prior genome mapping.
- Same-strand overlap with annotations and peaks, matching the prior R1 analysis.
- At least one aligned M, =, or X block; N and D do not supply coverage.
- Paired-end records cause an error instead of silently changing the count unit.
- Missing read sequence causes an error because motif membership cannot be known.

UAG is searched as **TAG in the original sequenced read orientation**. Reverse
strand BAM SEQ is reverse-complemented first. The search includes any stored
soft-clipped bases and insertions; it does not require the motif itself to align
inside a feature or peak. Hard-clipped sequence is unavailable. This tests a
motif in the read, not a genomic motif or the peak sequence. These orientation
and CIGAR rules follow the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

## Annotation choices

GTF exon models define introns as consecutive exon gaps. CDS plus annotated stop
codons define the coding category. Terminal exonic sequence before/after the
coding bounds defines 5UTR and 3UTR in transcript orientation. Explicit
five_prime_utr/three_prime_utr records are used when CDS is unavailable; otherwise
noncoding exons are `other`. Internal noncoding intervals of incomplete coding
models are not assumed to be terminal UTRs.

A read is annotated by at least one aligned nucleotide overlapping a feature.
The exclusive summary uses the declared precedence:

**CDS > 5UTR > 3UTR > intron > other.**

This is a bookkeeping convention, not a resolution of isoform identity.
`compatible_classes`, `gene_ids`, and `transcript_ids` preserve all overlapping
annotation memberships. Transcript IDs mean overlapping annotations, not a
whole-read splice compatibility assignment. Exclusive counts sum to the selected
library size; compatible counts may sum to more than that. Other includes
noncoding exons and unassigned/intergenic reads.

Introns with identical coordinates and strand are deduplicated. Identical
spliced UTR segment chains are deduplicated; distinct UTR models remain distinct.
Landmarks are deduplicated by type, position, and strand. Associated genes and
transcripts are retained. A read can contribute to several feature models or
landmarks; feature-bin counts must not be summed to estimate library size.

## Larger windows and scaled profiles

Fixed windows use transcript-oriented half-open coordinates:

- Splice donor and acceptor: **[-1000, +1000) nt**, in 10-nt bins.
- Annotated transcript 3′ end: **[-500, +500) nt**, in 10-nt bins.

Donor anchors are the first intronic nucleotide and acceptor anchors the last
intronic nucleotide in transcript orientation. The end anchor is the final
annotated transcribed nucleotide, not a validated cleavage/PAS annotation.
Minus-strand windows are reversed. Windows at chromosome boundaries are clipped;
actual bin lengths are retained rather than discarding those sites.

Length-scaled intron and 3UTR profiles use **100 bins from 0–100%**, oriented 5′
to 3′. A 3UTR is concatenated across its exons; intervening introns are excluded.
Bin edges are floor(i * feature_length / 100). Features shorter than 100 nt have
some empty bins; those bins have no nucleotide exposure and are excluded from
that bin's feature denominator. Integer edges ensure each nucleotide contributes
to only one bin in a particular feature.

Every feature-bin output contains both the count of reads overlapping that bin
and the number of aligned bases. One read contributes at most once to each
feature/bin count, even if multiple alignment blocks overlap it.

`metaprofiles.tsv` includes:

- Raw sums of read overlaps and aligned bases.
- Overlaps per million **selected** reads in that particular variant.
- Mean depth per feature = mean(aligned bases / actual bin length), including
  zero-coverage features with a nonempty bin.
- That same mean depth per million selected reads.
- Feature and nucleotide denominators for each bin.

The annotation feature universe is fixed across all samples and variants; there
is no separate per-condition coverage filter. Nonzero per-feature rows are
stored sparsely in `feature_bin_counts.tsv.gz`; omitted valid feature/bin rows
are exact zeros. Reconstruct them from `features.tsv.gz` and the bin rules in
this README and saved `analysis_config.yaml`. Zero-width bins are undefined,
not biological zero coverage. No nonzero bins are discarded.

## Supplied peaks

The included original `MSI2_peaks_omniCLIP.bed` is filtered at score >=1 and
merged for true overlaps on the same strand, using the original peak preparation
script. Adjacent peaks remain separate. Original peak membership is exported.
A read is counted once in the peak-overlapping sampling pool even if it overlaps
more than one peak. Per-peak counts count each such overlap, so their sum can
exceed the unique peak-overlapping library size.

`reference_peaks/` additionally contains static annotations and feature-bin
profiles for the merged peak intervals, counting each interval once per
feature/bin. These describe where the supplied peaks lie, using their full
intervals rather than summits. They are not sample-specific new peak calls and
cannot be UAG-filtered from a BED alone. The ten sample-specific analyses count
reads inside those supplied peaks. Thus static peak geometry is kept distinct
from changing read support.

## Configure and run on Bodhi

The sample manifest points to the original eight BAMs on Bodhi. Verify those
paths and the GTF path in `config/config.yaml`; the GTF must match the STAR genome
build and chromosome names. The workflow checks reference dictionaries across
all BAMs and rejects an entirely incompatible GTF. The configured annotation is
GENCODE v26 primary assembly, as in the existing project.

```bash
# After transferring read-feature-analysis.tar.gz to the project on Bodhi:
tar -xzf read-feature-analysis.tar.gz
cd read-feature-analysis
conda env create -f envs/snakemake.yaml
conda activate read-feature-snakemake
bash run_on_slurm.sh --dry-run
bash run_on_slurm.sh
```

The included profile submits SLURM jobs and creates rule environments with
Python and Samtools. Add your cluster account/partition if required. The analysis
per sample starts with a 24 GB memory and 24 h wall-time request; tune these after
a pilot. The sparse accumulator spills to SQLite under that sample's output
directory, limiting counter memory; sufficient disk space is required for the
per-read and feature-bin tables. It is removed after successful extraction.
Annotation models still reside in RAM. Per-sample jobs run independently.

To run just one sample after the shared census and annotation preparation:

```bash
bash run_on_slurm.sh results/NM2026_0081/QC.json
```

All eight census jobs are still required to establish the smallest library. A
zero minimum in any sampling pool stops the workflow for inspection rather than
silently producing empty equal-depth comparisons. Empty peak overlaps after a
genome-wide draw are permitted and normalized values become `NA`.

## Outputs

Shared outputs under `results/`:

- `sample_metadata.tsv`: original labels and sample/BAM provenance.
- `analysis_config.yaml`: settings used for this run.
- `sampling_targets.json`: each sample's full and UAG genome/peak census counts,
  reference dictionaries, and exact equal-depth targets.
- `features.tsv.gz`: coordinate/strand/gene/transcript catalog and bin geometry.
- `reference_peaks/`: filtered/merged BED, source membership, peak annotations,
  sparse peak/feature/bin table, and static peak metaprofiles.
- `README.md`: interpretation and coordinate definitions.

Each sample subdirectory contains:

- `read_annotations.tsv.gz`: every eligible read once, motif status, annotation,
  overlapping peaks, and selection membership in all ten variants.
- `category_counts.tsv`: exclusive and compatible feature-category counts,
  fractions, and counts per million.
- `feature_bin_counts.tsv.gz`: sparse raw feature/bin read and aligned-base counts.
- `metaprofiles.tsv`: aggregate raw and depth-scaled fixed/scaled profiles.
- `peak_counts.tsv.gz`: counts for every merged peak and variant, including zeros.
- `QC.json`: actual selected counts, targets, eligibility settings, and seed.

Copy the entire results directory back for further analysis. Gzip TSV files can
be read directly in R or Python. Join each sample_id to sample_metadata.tsv;
biological replicate labels are never collapsed.

## Validation

```bash
python -m unittest discover -s tests -v
```

Tests cover exact/reproducible sample sizes, the interval index against brute
force, CIGAR-spliced overlap, coding/UTR boundaries including stop codons,
reverse-strand motifs, minus-strand and short-feature scaled bins, clipped
windows, and an end-to-end two-library fixture covering all ten variants. The
unit fixture uses a SAM-text adapter for Samtools. A separate full Snakemake
integration check using real BAM files and the Samtools implementation bundled
with pysam also passed, including a clean second dry-run:

```bash
# Optional local validation environment; not needed for the cluster analysis.
python -m pip install 'snakemake>=9,<10' pysam
python tests/validate_workflow.py /tmp/read-feature-validation
```

A chromosome 21 check using the project's actual GENCODE v26 GTF constructed
16,583 features and confirmed that all profile bins conserve their feature's
nucleotide length. The full biological BAMs have not been processed locally.
A cluster pilot remains necessary for production runtime and full-library
verification.

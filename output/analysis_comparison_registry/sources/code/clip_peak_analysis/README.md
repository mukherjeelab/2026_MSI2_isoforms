# MSI2 fixed-peak CLIP enrichment

This directory contains the reproducible downstream analysis used to count the
eight legacy-processed CLIP/input BAMs in a fixed set of previously called
omniCLIP peaks and run paired DESeq2 comparisons. It does not rerun FASTQ
preprocessing, alignment, peak calling, or T-to-C crosslink analysis.

## Experimental design

The BAMs are counted together with featureCounts. DESeq2 is then run separately
for the two constructs using `~ pair + assay`:

- `MSI2_328WT`: NM2026_0081/0082 (pair A) and NM2026_0083/0084 (pair B)
- `MSI2_324WT`: NM2026_0085/0086 (pair A) and NM2026_0087/0088 (pair B)

Positive `log2FoldChange` means IP enrichment relative to input. Size factors
come from total primary mapped alignments in each BAM (`full_library`
normalization), rather than from the external-peak count matrix.

The peak preparation step filters the 58,555 source peaks at score >= 1, retains
46,108, and merges same-strand overlaps into 43,794 non-overlapping regions.
Original peak provenance and scores remain in `peak_members.tsv`.

## Repository layout

```text
data/MSI2_peaks_omniCLIP.bed               original fixed peak set
data/results/peak_analysis/                 downloaded counting/DESeq2 outputs
code/clip_peak_analysis/config.yaml         analysis parameters and BAM path
code/clip_peak_analysis/deseq_samples.tsv   sample annotations and pairing
code/clip_peak_analysis/Snakefile           reproducible workflow
code/clip_peak_analysis/scripts/            peak preparation and DESeq2 code
```

The BAMs are intentionally not stored in Git. `config.yaml` points to their
current location on Bodhi.

## Run on Bodhi

From the repository root:

```bash
conda env create -f code/clip_peak_analysis/envs/snakemake.yaml
conda activate clip-snakemake

snakemake --snakefile code/clip_peak_analysis/Snakefile \
  --profile code/clip_peak_analysis/profiles/slurm --dry-run

snakemake --snakefile code/clip_peak_analysis/Snakefile \
  --profile code/clip_peak_analysis/profiles/slurm
```

Confirm `differential_analysis.bam_pattern` in `config.yaml` before rerunning.
BAM indexes are not required by this workflow.

## Inspect the downloaded results

Run the dependency-light audit from the repository root:

```bash
Rscript code/clip_peak_analysis/scripts/audit_results.R
```

It validates row/sample identity and count totals, then writes compact summaries
to `data/results/peak_analysis/audit/`.

See `data/results/peak_analysis/README.md` for the current result interpretation
and limitations.

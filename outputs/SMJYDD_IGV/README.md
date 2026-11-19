# View the SMJYDD reads in IGV

All reads in this package are at least 500 nt long. Only the reference(s) assigned to each sample are used. The original data and count tables are unchanged.

## Open a prepared session in IGV Desktop

An IGV application was found at `/Applications/IGV_2.8.13.app`. Automatic launch/inspection was blocked by missing computer-control permissions; sessions have been checked as files but not visually tested in IGV.

1. Open IGV.
2. Choose **File > Open Session**.
3. Select one of the XML files in this directory:

| Session | What it shows |
|---|---|
| `L4_on_324_altUTR.xml` | L4 HEK and H295R against 324, zoomed around the 54-base deletion. Start here. |
| `L2_on_328_ref.xml` | Both L2 samples against 328, around the candidate 63-base deletion. Insertion-containing reads can appear with insertions. |
| `L2_on_346_MSI2.xml` | Both L2 samples against 346. 328-like reads can appear with a large deletion over the inserted region; inspect insertion-containing reads for the three-base junction difference. |
| `L8_on_251_Exon1.xml` | Both L8 samples against 251, initially showing the full reference. |

Keep the XML, FASTA, FASTA index, BAM, BAM indexes, and BED file together. Session paths are absolute paths on this Mac, including explicit BAM index paths, to avoid index resolution errors in the installed IGV version. If moving the folder, regenerate the sessions or load the files manually.

These are construct sequences, not genomic chromosome coordinates. Do not choose hg38 or another human genome: the session loads the supplied custom reference FASTA.

## Navigate to a region

Paste a location into the IGV location/search box:

- L4 deletion: `324_altUTR:790-970`
- L4 upstream difference: `324_altUTR:90-150`
- L2 63-base deletion: `328_ref:650-790`
- L2 insertion and three-base junction difference: `346_MSI2:720-850`
- L8 low-frequency 49-base deletion: `251_Exon1:370-480`

Zoom in to see individual bases. Read rectangles show aligned sequence; black deletion bars/gaps indicate bases absent relative to the displayed reference. Colored base mismatches and insertion markers identify other differences. Click a read for its identifier, CIGAR and tags. Candidate regions are marked in an annotation track; these labels are hypotheses, not validated biological calls.

Right-click an alignment track to expand the reads. IGV can downsample deep tracks, so the number of rows visible is not a read count. Unmapped records are preserved in the BAMs but will not appear at a reference locus. L4 HEK consequently has no mapped rows under these alignment criteria.

## How the display relates to the count tables

Each BAM is an **independent single-reference view**, with one best qualifying alignment per retained read. L2 reads are deliberately aligned once to 328 and separately to 346 so you can inspect both populations on either coordinate system. **Do not add reads across the two views or treat a displayed alignment as a unique construct assignment.** MAPQ here is relative to the single-reference view, not confidence distinguishing 328 from 346.

The previous sample-specific competitive annotation call is preserved as BAM tag `ZC`; the target-window call is tag `ZL`. `RG` identifies the sample. These tags let you cross-reference the saved count tables. Mapping requires at least 80 aligned query bases and 85% identity, matching the previous mapping threshold; reads must independently be >=500 nt. Counts per display are in `view_counts.tsv`. This package does not revise biological conclusions or count assignments.

Reverse-strand sequences and qualities are stored in reference orientation, with flag 16, appropriate soft clipping, and the mapper's CIGAR, NM and MD. Supplementary/secondary alternatives are not included; each reference view keeps one qualifying best alignment per read. This keeps the display readable but is not a comprehensive structural-rearrangement analysis.

## Manual loading / IGV-Web

If a session does not open in the installed IGV version, choose **Genomes > Load Genome from File**, select `references.fasta`, then **File > Load from File** and select the BAM(s) and `regions_of_interest.bed`. Desktop IGV finds neighboring `.bam.bai` indexes automatically.

For IGV-Web at https://igv.org/app, load the reference FASTA and its `.fai` using the Genome local-file option, then load the desired BAM **and its `.bam.bai` together** from local files using Tracks. Use the same coordinates above. Follow the current IGV-Web prompts if menu wording differs.

## Reproduction and validation

1. System Python: `python3 -B analysis/prepare_smjydd_igv.py`.
2. Bundled Python 3.12: run `analysis/convert_smjydd_igv_bam.py`. It uses the project-local pysam installation in `tmp/sequence_browser_py312`.

Sorted BAMs and indexes were checked for record counts, unique read names within each view, minimum read length, CIGAR/query-length consistency, reference identity, and indexed retrieval. See `validation.txt`. The BAMs contain all retained records for that view, including unmapped reads, and do not contain reads shorter than 500 nt. The FASTA contains all four original expected sequences; alignments for each sample are restricted to the appropriate reference view.

Documentation:
- https://igv.org/doc/desktop/UserGuide/sessions/
- https://igv.org/doc/desktop/UserGuide/reference_genome/
- https://igv.org/doc/desktop/UserGuide/tracks/alignments/viewing_alignments_basics/
- https://igv.org/doc/webapp/UserGuide/

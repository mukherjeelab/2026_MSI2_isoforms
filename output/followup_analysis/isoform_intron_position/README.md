# Equal-intron positional profiles: endogenous summits and isoform CLIP reads

Run from the project root:

```sh
Rscript output/followup_analysis/isoform_intron_position/run_analysis.R
```

## Figures

- `endogenous_intron_equal_weight.png/pdf`: the original left panel alone, retaining its exact 20,364 reference summit assignments in 9,565 introns strictly longer than 500 nt and all 500 cached null draws. Equal intron **weight** means each represented intron contributes total weight one; the introns do not have equal physical lengths.
- `isoform_read_positions_all_assigned.png/pdf`: main raw-alignment positional analysis, with separate 324/328 panels. Thin red dashed/dotted lines are replicates A/B; thick red is their arithmetic mean. Each replicate is normalized separately before averaging, so the deeper library does not dominate.
- `isoform_read_positions_reference_introns.png/pdf`: restrict the established read assignments to the exact intron IDs represented in the endogenous >500-nt analysis. Unsupported introns cannot contribute a normalized read distribution and are omitted separately per library.
- `isoform_read_positions_common_four_libraries.png/pdf`: restrict to introns with at least one assigned read center in every CLIP library, controlling the intron population across the isoform/replicate comparison. This subset is defined by support, not by positional enrichment.

All profiles use 100 unsmoothed transcript-oriented bins, with 0% at the donor and 100% at the acceptor. Blue shading represents pointwise null percentiles, not P-value thresholds or biological-replicate confidence intervals.

## Alignment source and eligibility

The BAM paths in `data/results_round2/sample_metadata.tsv` point to `/beevol/home/bartholk/2026_MSI2_CLIP_iso/clip-preprocessing/results/alignment/` on Bodhi. BAM files are not present in this local checkout. This run uses the existing full per-alignment exports `data/results_round2/<sample>/read_annotations.tsv.gz`, containing reference start, strand, and CIGAR for every eligible BAM alignment. It does **not** reopen or rerun BAM filtering. The exported row counts are checked against each original `QC.json` census.

Libraries: 324 A = NM2026_0085; 324 B = NM2026_0087; 328 A = NM2026_0081; 328 B = NM2026_0083. No SMI data, peak selection, UAG selection, equal-depth downsampling, or library pooling is used.

The original BAM workflow had already removed repeat reads, collapsed UMI+insert duplicates, and selected unique mappings. Eligible exported records are single-end primary alignments with MAPQ >=10 and NH absent or <=1. Each collapsed record is counted once; read-name abundance suffixes are not multiplied back into counts. These are processed genome-aligned reads, not unprocessed FASTQ molecules.

## Read position and feature assignment

A read contributes one position: its central **aligned nucleotide**, reconstructed from CIGAR M/= /X blocks. Insertions, deletions, skipped reference intervals, and soft clips do not count toward aligned length. For an even aligned length, select the 5'-central of the two central nucleotides in read/transcript orientation. This prevents the midpoint from falling into a skipped CIGAR N interval. Use the center of that nucleotide in continuous feature coordinates: `(offset0 + 0.5) / length_nt`.

Read centers are strand-matched to the same GENCODE v26 intron definitions as the reference analysis. As with reference summits, choose the longest matching intron per point, break ties by original feature-table order, then retain introns >500 nt. A read is counted once, even when transcript models overlap. A read with no eligible intronic center does not contribute, even if a different part of its alignment touches an intron. Identical-coordinate introns were deduplicated upstream; distinct overlapping models remain available for assignment. Reference-subset filtering is performed **after** assignment and does not reassign reads.

For a library, a center in an intron containing n assigned centers receives weight 1/n. Sum weights in each bin and divide by the number of supported introns. Every library profile sums to one. Arithmetic means of the A/B library profiles also sum to one. This normalization measures positional shape, not total intronic binding strength or expression. It gives one-read introns the same total weight as deeply covered introns, matching the requested equal-intron weighting.

## Randomized background

For every library and cohort, preserve exact assigned intron identities, lengths, strands, and read counts. Independently randomize each point on the prior half-integer grid:

```r
12.5 + floor(runif(n) * (length_nt - 25 + 1))
```

Use 500 iterations and a fixed initial seed of 260922. Iteration-matched A/B null profiles are averaged with equal replicate weights; the plotted blue cloud and percentile bands describe these 500 null **averages**. Individual-library nulls are retained in the RDS. The endogenous figure uses its unchanged original cached randomizations.

This intentionally preserves the existing PARPIPE-style **point-position** null with a 25-nt implicit footprint and 12.5-nt end margins. It is not execution of the complete historical PARPIPE script, a simulation of full CIGAR alignments, or a read-length-aware footprint model. Observed centers near an end are retained just as observed summits were; their counts are reported in `read_assignment_qc.tsv`. The null does not account for mappability, sequence bias, transcript abundance, or correlations between reads. Null intervals therefore describe this positional model, not uncertainty between biological replicates. No isoform-difference P values are claimed from these read-level simulations.

Endogenous summit positions and isoform read centers are different measurements from different experiments. Similar donor-proximal profiles support a descriptive positional similarity, not identical binding specificity, cooperative binding, or a statistically established difference between isoforms.

## Reproducible outputs and checks

- `assigned_read_centers.tsv.gz`: every assigned read, center coordinate, chosen feature, length, relative position, and library identity.
- `candidate_introns.tsv.gz`, `cohort_memberships.tsv.gz`: annotation and exact cohort membership.
- `read_assignment_qc.tsv`: exported alignment totals, assigned-center counts, supported intron counts, and end-margin counts.
- `read_position_profiles.tsv`: observed replicate/mean curves and mean-null percentile plotting inputs.
- `read_position_nulls.rds`: all 500 null curves and observed profiles for every library and cohort, including support counts and seed.
- `replicate_mean_null_curves.tsv.gz`: all plotted null curves in tabular form.
- `positional_zone_summary.tsv`: equal-intron-weighted first 10%, first 20%, and last 10% fractions, null expectations/95% intervals, and observed/null ratios for every replicate and its isoform mean.
- `endogenous_plot_inputs.tsv`: unchanged observed reference curve and cached null summaries.
- `input_manifest.tsv`, `sessionInfo.txt`: source checksums and R package versions.

The script verifies exported alignment census counts and unique read names, contiguous intron coordinates, one valid center per read, one feature per retained read, normalized observed/null profiles, and reference-cache median agreement. Existing analyses and inputs are not changed.

## Results from this run

Both isoforms have donor-proximal enrichment in both replicates. In all assigned introns, the first 10% contains 13.73% of the normalized 324 profile (1.40-fold null) and 15.35% of the 328 profile (1.57-fold null). The endogenous equal-intron summit profile has 16.94% in this zone. These percentages are averages of within-intron positional fractions, not unweighted fractions of all reads.

The reference-intron subset gives 15.69% for 324 (1.58-fold null) and 17.90% for 328 (1.81-fold null). When restricted to the same 1,315 introns supported in all four CLIP libraries, the first-10% fractions are 16.84% for 324 (1.70-fold null) and 16.44% for 328 (1.66-fold null). Thus the shared donor-proximal trend is robust to this feature-population check; the apparent 328-versus-324 ordering in the full cohorts is not robust and should not be presented as an established isoform difference.

An independent Python/pandas check reconstructed all 12 observed replicate curves from saved read assignments and verified their agreement with the R plotting inputs to 1e-12. It also checked all six plotted means against the arithmetic mean of their two replicate curves. All four PNG layouts were visually inspected.

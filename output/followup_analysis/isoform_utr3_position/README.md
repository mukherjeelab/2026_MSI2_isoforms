# Equal-feature 3′UTR positional analysis

Reproduce from the project root:

```sh
Rscript output/followup_analysis/isoform_utr3_position/run_analysis.R
python3 output/followup_analysis/isoform_utr3_position/summarize_validate.py
```

## Figures and counts

Every figure is saved as PNG and PDF.

- `endogenous_utr3_equal_weight`: 26,504 reference summits in 8,869 represented 3′UTR features. Uses the existing summit assignments with equal feature weighting; this is a new equal-feature rendering and 500 new null draws, rather than a reuse of the earlier site-weighted 3′UTR curve.
- `isoform_read_positions_all_assigned`: all eligible assigned read centers, with replicate A/B and their equally weighted average in separate 324/328 panels.
- `isoform_read_positions_reference_utr3`: filter existing read assignments to reference-summit feature IDs, without reassigning reads. A reference feature contributes only where it has read support in that library.
- `isoform_read_positions_common_four_libraries`: the same 1,474 features with read support in all four CLIP libraries.

| Library | All assigned: reads / features | Reference subset: reads / features | Common subset: reads / features |
|---|---:|---:|---:|
| 324 A (0085) | 22,879 / 3,205 | 18,880 / 2,403 | 14,795 / 1,474 |
| 324 B (0087) | 29,589 / 3,679 | 24,366 / 2,713 | 18,343 / 1,474 |
| 328 A (0081) | 77,976 / 7,095 | 69,993 / 5,182 | 44,920 / 1,474 |
| 328 B (0083) | 88,739 / 6,534 | 80,217 / 4,892 | 51,481 / 1,474 |
| Union across libraries | 219,183 / 9,844 | 193,456 / 6,403 | 129,539 / 1,474 |

Feature counts are unique feature IDs, not genes. Features overlap across libraries, so counts are not additive. Distinct transcript feature models may overlap genomically. Read counts are distinct eligible collapsed alignment records, summed across libraries, not pre-collapse molecule abundances. `analysis_counts.tsv` additionally reports each isoform's union across replicates; `endogenous_counts.tsv` reports the reference population.

## Methods

Retain the existing 3′UTR minimum **spliced length >=100 nt**, not the intron-specific >500-nt restriction. Use hg38 GENCODE v26 features from `data/results_round2/features.tsv.gz`. For multi-exon 3′UTRs, concatenate annotated exonic segments and exclude intervening introns; reverse coordinates on the minus strand. Position 0% is the annotated beginning of the 3′UTR and 100% its transcript-oriented end. The script independently reconstructs all original reference feature choices and positions and requires exact agreement within numerical precision.

As in the intron analysis, one point per read is its central aligned nucleotide from CIGAR M/= /X blocks. For even aligned lengths, select the 5′-central nucleotide. Soft clips, insertions, deletions, and skipped reference intervals are excluded from aligned length. Select the longest matching strand-compatible eligible 3′UTR; ties follow original feature-table order. Each read is assigned to at most one 3′UTR here. Intron and UTR analyses are separate, so classifications across those analyses need not be mutually exclusive across overlapping transcript models.

The original BAMs reside on Bodhi; this run uses complete saved per-alignment exports in `data/results_round2/<sample>/read_annotations.tsv.gz`, not reopened BAMs. The exports preserve starts, strand, and CIGAR for all eligible records. Original filtering was single-end primary, MAPQ >=10, NH absent or <=1, after repeat subtraction, UMI+insert collapse, and unique mapping. Counts are checked against each original QC.json census. No SMI data, peak restriction, UAG filtering, CPM scaling, or downsampling is applied here.

Each supported feature contributes total weight one: each assigned center has weight 1/(centers assigned to its feature), then sum per 1% bin and divide by the number of supported features. Zero-support features are omitted because a within-feature positional distribution cannot be estimated. Normalize each replicate separately, then average A/B equally. All 100 bins are retained without smoothing. This measures positional shape, not total binding strength or expression.

For each library/cohort, preserve exact feature identity, spliced length, strand, and number of assigned centers. Draw positions independently on the existing PARPIPE-style half-integer grid `12.5 + floor(runif(n)*(length_nt-25+1))`, giving the same 12.5-nt end margins / 25-nt implicit footprint convention. There are 500 iterations with seed 260922. Replicate-null curves are averaged iteration by iteration; these 500 mean curves are the cyan cloud. The endogenous simulation and the read simulations each initialize the fixed seed independently. Dark/light bands are pointwise 25–75 / 10–90 percentiles; black is the null median. Blue does not encode P values or biological-replicate uncertainty.

This is the previous point-position adaptation, not the full historical PARPIPE algorithm or a simulation of CIGAR alignments with actual read lengths. Observed edge centers remain included; `read_assignment_qc.tsv` reports those within 12.5 nt of an end. The null does not model mappability, sequence bias, transcript abundance, or read correlation. No biological isoform-comparison P values are claimed.

## Findings

The all-assigned 324/328 average profiles put 12.30%/12.49% of their equal-feature-weighted mass in the first 10% of the 3′UTR, versus 9.23%/9.10% null expectation (1.33-fold/1.37-fold). This proximal enrichment is present in both replicates and remains in the reference and common-feature subsets. The final 10% is depleted relative to this null (6.35%/6.93% observed). Thus these read-center profiles do not reproduce the strong distal-end rise seen in the endogenous summit profile. Summits and aligned-read centers from different experiments are different measurements; this descriptive contrast alone does not establish a biological change in binding specificity.

## Saved reproducibility files

`assigned_read_centers.tsv.gz`, `candidate_utr3.tsv.gz`, `feature_segments.tsv.gz`, and `cohort_memberships.tsv.gz` contain assignments, segmented coordinates, and exact cohorts. `read_position_profiles.tsv` and `replicate_mean_null_curves.tsv.gz` contain plot inputs. `read_position_nulls.rds` contains all per-library 500-iteration nulls and support counts. `endogenous_assigned_summits.tsv.gz`, `endogenous_plot_inputs.tsv`, and `endogenous_equal_feature_nulls.rds` preserve the reference analysis. `positional_zone_summary.tsv` contains first-10%, first-20%, and last-10% fractions, null expectations/95% intervals and ratios. Input hashes and software versions are in `input_manifest.tsv` and `sessionInfo.txt`.

The R script checks original alignment totals, unique read names, segmented lengths, reference-assignment agreement, valid centers, and unit-sum profiles. Independent Python validation reconstructs every read's spliced offset, all 12 replicate profiles, six replicate means, and the reference equal-feature curve to tolerance 1e-12. See `validation.txt`. All four PNGs were visually inspected. Earlier analyses remain unchanged.

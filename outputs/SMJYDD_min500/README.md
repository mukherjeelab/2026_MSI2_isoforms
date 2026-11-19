# SMJYDD: reads of at least 500 nt only

All current calculations exclude reads shorter than 500 nt. Exactly 500 nt is retained. Only the assigned sample-specific references are used. There are 10,354 retained reads across the six samples. These results replace the earlier unfiltered summaries.

## Length and annotation mapping

Percentages in this table use **all retained reads in each sample**, including ambiguous and nonqualifying reads. Mapping support means the aligned portion resembles that reference; it does not establish an intact, exact full reference. Recurring variants below are included in these mapping counts.

| Sample | Retained reads | Mean length, nt | Median length, nt | Assigned annotation: mapping support | Ambiguous | No qualifying alignment |
|---|---:|---:|---:|---|---:|---:|
| L2 HEK | 1,431 | 843.02 | 830 | 328: 1,244 (86.93%); 346: 80 (5.59%) | 23 | 84 |
| L2 H295R | 1,038 | 780.54 | 817.5 | 328: 717 (69.08%); 346: 10 (0.96%) | 65 | 246 |
| L4 HEK | 29 | 727.48 | 607 | 324_altUTR: 0 (0%) | 0 | 29 |
| L4 H295R | 382 | 947.98 | 1,013 | 324_altUTR: 282 (73.82%) | 0 | 100 |
| L8 HEK | 4,677 | 811.85 | 745 | 251_Exon1: 4,527 (96.79%) | 0 | 150 |
| L8 H295R | 2,797 | 814.51 | 745 | 251_Exon1: 2,570 (91.88%) | 0 | 227 |

The mapping criteria remain at least 80 aligned query bases and 85% identity. L2 reference preference requires a score advantage of at least 20; shared/near-tied reads remain ambiguous. The minimum **read** length is now 500; this should not be confused with a minimum 500-base **alignment**. Coverage and identity for every retained read are in the per-read file.

## Intact versus deleted target region

The stricter target-window assay requires both flanks and therefore has a different denominator from mapping support.

| Sample | Insertion-containing | Insertion-absent | Other/conflicting | Not callable | Insertion-containing % of resolved |
|---|---:|---:|---:|---:|---:|
| L2 HEK | 68 | 1,154 | 9 | 200 | 5.56% |
| L2 H295R | 4 | 532 | 4 | 498 | 0.75% |
| L4 HEK | 0 | 0 | 0 | 29 | Not estimable |
| L4 H295R | 29 | 234 | 0 | 119 | 11.03% |

For L2 these classes mean 346-like and 328-like at the insertion site, respectively, not exact full-reference matches. For L4, insertion-absent means a deletion within 324, not assignment to 328. L8's assigned reference does not contain the target region, so that assay is not applicable.

In L4 H295R, 29/382 retained reads (7.59%) support retention of the target region and 234/382 (61.26%) support its deletion. Among the 263 resolved reads, those fractions are 11.03% and 88.97%.

## Recurrent differences worth investigating

Coordinates below are 1-based inclusive in the named user-supplied reference. Indels are left-normalized; repeated sequence can allow equivalent shifted coordinates. Counts in this section use a more stringent variant screen: a uniquely preferred allowed reference, at least 90% alignment identity and at least 80% reference span, one best alignment per read. Indels require aligned flanks on both sides. Substitutions require base quality Q20. Terminal ten bases are excluded. These are descriptive sequence candidates, not proof of biological origin, and are not an exhaustive variant discovery analysis.

### 1. A 63-base deletion within 328 in both L2 samples

- Representation: **328_ref bases 689–751 deleted**.
- L2 HEK: **57/1,197 evaluable screened reads (4.76%)**, 35 forward and 22 reverse.
- L2 H295R: **158/701 (22.54%)**, 85 forward and 73 reverse.
- Independent exact 40-base deletion-junction searches found **58 HEK and 168 H295R reads** among all retained reads, corroborating the recurrent junction without relying on the mapper. These searches have different eligibility criteria and should not be added to the mapping counts.

This is an additional recurring shorter 328-like sequence, especially in H295R, and is a strong candidate for follow-up. It was not distinguished from ordinary 328 mapping in the earlier annotation-level summaries.

### 2. The insertion-containing L2 reads usually have a three-base-shorter junction than supplied 346

- Representation: **346_MSI2 bases 750–752 (ATC) deleted**; equivalent placement is possible within the repeated junction.
- HEK: **68/74 evaluable 346-like screened reads (91.89%)**, on both strands.
- H295R: **5/7 (71.43%)**; small denominator.
- Exact deletion-junction support: **60 HEK and five H295R reads**. No retained read has an exact match to the corresponding intact reference window with 20-base flanks.
- Expected local sequence: `...AGGATCATCAGTC...`; common observed sequence: `...AGGATCAGTC...`.

The reads support an insertion-containing 346-like population, but its dominant local sequence differs from the supplied annotation. Check the reference/construct sequence record before treating this as a newly arisen biological event.

### 3. L4 H295R has the previously identified 54-base deletion, plus a recurring upstream difference

- Main deletion: **324_altUTR bases 851–904** in left-normalized representation, equivalent to the previously reported **853–906** deletion. The target-window assay finds **234 deletion-supporting and 29 inclusion-supporting reads** after the length filter. The stricter variant screen identifies the exact normalized 54-base CIGAR event in **228/260 reads (87.69%)**, with 122 forward and 106 reverse reads. Different counts reflect different criteria, not additional reads.
- Upstream difference: **bases 124–125 (GT) deleted**, in **130/256 screened reads (50.78%)**, split 65 forward/65 reverse. This changes local `...CAGGGTTAAAA...` to `...CAGGTAAAA...`. Independent exact local-window matches give 109 deletion-window and 109 reference-window reads.

The upstream two-base difference is recurrent and merits inspection, but its short repetitive context means it should not automatically be called biological. The older one-base difference in the supplied consensus does not capture this read-level mixture cleanly.

A one-C deletion at position 78 also appears in 30 screened reads, but is strongly strand-imbalanced (26 forward/four reverse); it is a lower-priority signal. Variant combinations have not been reconstructed into validated full-length haplotypes; the candidate FASTA contains one change at a time.

### 4. L8 predominantly matches 251, with low-frequency deletion signals

No major recurring alternative dominates either L8 sample in this screen. A **49-base deletion at 251_Exon1 bases 404–452** occurs in 29/4,393 screened HEK reads (0.66%), supported on both strands. An independent exact deletion-junction search finds **54/4,677 HEK reads (1.15%)** and **14/2,797 H295R reads (0.50%)**. This is a low-frequency candidate requiring confirmation, not an established isoform.

Other smaller, low-frequency deletions are listed in `recurring_variants.tsv`, including a four-base deletion at 102–105 in 23 HEK and 20 H295R screened reads. Terminal A differences from the earlier consensus comparison were excluded from this internal-variant screen and are not strong evidence of a biological alternative.

### 5. L4 HEK remains unresolved

Only 29 reads survive the length filter, and none meet the mapping criteria against 324. The present data do not support naming a specific biological alternative for this sample.

## Reproduction and validation

Run from the project root:

```bash
python3 -B analysis/summarize_smjydd_min500.py
python3 -B analysis/check_smjydd_candidate_junctions.py
```

Dependencies and mapping settings are inherited from the preceding analysis (mappy/minimap2 2.31, map-ont; biopython 1.83). Raw FASTQs are unchanged. All 10,354 retained read records are unique within samples and are >=500 nt; counts reconcile, and means were independently recalculated from the per-read table. Candidate junction searches inspect sequence strings in both orientations, independently of mapping.

- `sample_summary.tsv`: retained lengths, annotation mapping and locus counts.
- `read_assignments.tsv`: retained reads only, with reference calls, coverage, locus calls, and screened events.
- `recurring_variants.tsv`: recurrent normalized alignment events supported by at least three screened reads; includes counts, denominators and strand support. This is a discovery table, not a validated call set.
- `exact_junction_validation.tsv`: independent exact candidate-junction counts in retained reads.
- `candidate_alternatives.fasta`: single-deletion candidate reference sequences for follow-up. These are modified reference hypotheses, not assembled or validated full-length biological sequences.

Standing settings are recorded in `analysis/SMJYDD_ANALYSIS_NOTES.md`. Candidate variants may reflect biological variation, PCR products/artifacts, sample mixtures, or sequencing/alignment effects; this analysis does not establish their source.

# SMJYDD comparison against four expected references

Analysis date: 2026-10-16. Platform supplied by user: Plasmidsaurus Standard Premium PCR (Oxford Nanopore).

## Findings

The low exact-target counts are not explained solely by sequencing errors. Most callable reads in L2 and L4 H295R support an allele lacking the target-containing insertion. Many other reads do not provide sufficient alignment/coverage for a target-locus call.

| Sample | All reads | Target-containing, spanning | Target-absent, spanning | Other/conflicting locus | Not callable | Target-containing % among present/absent calls |
|---|---:|---:|---:|---:|---:|---:|
| L2 HEK | 5,332 | 75 | 1,192 | 8 | 4,057 | 5.92% |
| L2 H295R | 5,427 | 4 | 555 | 5 | 4,863 | 0.72% |
| L4 HEK | 4,524 | 0 | 0 | 0 | 4,524 | Not estimable |
| L4 H295R | 4,216 | 30 | 237 | 0 | 3,949 | 11.24% |
| L8 HEK | 5,476 | 0 | 0 | 0 | 5,476 | Not estimable |
| L8 H295R | 5,047 | 0 | 0 | 0 | 5,047 | Not estimable |

“Target-containing” allows substitutions and small indels; it does not require the original 57 bases to be error-free. “Target-absent” requires alignment-supported sequence on both sides of the locus. “Not callable” includes insufficient coverage, failed mapping, and failed flank checks; it does not imply biological absence. L8 reads support a different ending represented by 251_Exon1 and cannot be evaluated with the shared-flank target assay. These are read proportions from PCR products, not estimates of original cell/transcript fractions or unique starting molecules.

### Reference and consensus comparison

- **324_altUTR:** 1,070 nt, contains the target. The insertion between the shared flanks is 54 nt because the target's initial TCA overlaps the left flank.
- **328_ref:** 830 nt, lacks the insertion. Both L2 supplied consensuses match this sequence exactly.
- **346_MSI2:** 887 nt, contains a 57-nt insertion relative to 328_ref. The target alone does not distinguish it from 324_altUTR; other reference differences are needed.
- **251_Exon1:** 744 nt, has a distinct ending and does not include the target region. L8 HEK's consensus differs only by an extra terminal A; L8 H295R's lacks the terminal A. Their closest-reference read calls strongly favor 251_Exon1 (4,731 and 2,686 reads, respectively).
- **L4 H295R:** its supplied consensus is reverse-complemented relative to the expected references. After orientation correction it differs from 324_altUTR by a 54-nt deletion across the target region and a separate one-base deletion near the 5′ end. The exact placement of gaps within short repeats is not unique. In an augmented comparison, 253 reads preferentially support this observed consensus, versus 32 supporting expected 324_altUTR. These whole-reference calls use a different denominator from the stricter locus-spanning counts.
- **L4 HEK:** no consensus FASTA was supplied in the results directory, and no read meets the stated reference-match criteria. Median read length is 138 nt; only nine reads are at least 700 nt. This sample is inconclusive for the intended constructs, rather than demonstrated target-negative.

## Recommended next steps

1. Confirm which expected construct belongs to each sample and check the actual PCR primer sequences against the expected amplicons.
2. If L4 H295R was intended to retain the insertion, inspect its construct/clone and PCR product: both the supplied consensus and individual reads support a dominant deletion allele. Consider independent confirmation with a fresh PCR and sequencing across both junctions.
3. Check the L4 HEK PCR product size and specificity before repeating sequencing. Its short reads and lack of qualifying alignments do not identify the intended construct.
4. Interpret minority counts cautiously: PCR bias, sample mixtures, chimeras and sequencing artifacts can affect read fractions. This analysis cannot determine their biological source.

## Methods and limits

All 30,022 FASTQ records were analyzed, without quality filtering, trimming, deduplication, or changes to the supplied FASTQs. The exact user reference sequences are saved in `data/SMJYDD_expected_references.fasta`.

Mapping used minimap2 through mappy 2.31 with the `map-ont` preset. Each reference was indexed separately, retaining up to ten hits. This avoids forcing shared reads onto whichever reference happens to receive a primary alignment in a joint index. Whole-reference calls require at least 80 aligned query bases and alignment identity at least 85%. A unique reference preference requires a score advantage of at least 20 over other references. Otherwise the read is ambiguous. Scores are reconstructed from CIGAR with the map-ont match/mismatch and dual-affine gap costs. These calls describe the closest supplied sequence in the aligned portion, not proof of a complete, error-free construct. No reference-covering percentage requirement is imposed; the per-read file includes that percentage.

An explicitly labeled diagnostic reference, `observed_L4_H295R_consensus`, was added in a second comparison. Expected-only and augmented calls are both retained. Its support is descriptive and is not independent validation of a consensus built from the same read data.

For the target assay, reads are mapped to 328_ref and the window covering 1-based positions 714–793 is extracted using the CIGAR, in the reference orientation. This comprises 40 bases on either side of the insertion site after position 753. Insertions in reads are retained. Both outer 30-base flanks must have at least 85% edit-distance identity. The complete extracted window is compared by global edit distance against the absent allele and the two supplied insertion alleles, avoiding sensitivity to equivalent gap placements within repeats. A target-present call allows at most 14 edits to a 134/137-base insertion-containing window; an absent call allows at most eight edits to the 80-base absent window. The winning allele class must be better by at least ten edits. Multiple qualifying hits with different classes yield a conflicting call. Thus these are conservative, alignment-based support counts with explicit thresholds, not validated probabilistic genotypes.

Synthetic checks passed for forward/reverse sequences, both insertion alleles, the absent allele, a single substitution/deletion, and incomplete coverage. All sample totals reconcile, read IDs are unique within samples, and exact-target counts reproduce the original analysis. Raising flank identity to 90% yielded present/absent counts of 76/1,179, 4/549, and 30/234 for L2 HEK, L2 H295R, and L4 H295R. One conflicted L2 HEK read becomes present when a poorer hit fails the stricter check. A smaller-seed mapping sensitivity check (k=11, w=5, minimum chain score 20, minimum minimizer count 2) retained all present calls and changed only one absent call and one other call to uncallable. The main interpretation is stable to these checks.

## Files and reproduction

- `sample_summary.tsv`: counts and read-length summaries, expected-only and augmented calls, and stricter flank checks.
- `read_assignments.tsv`: every read's reference comparison and target-locus classification.
- `alignments.paf`: mappings to each of the four references and the observed consensus, with CIGAR and cs tags. Query names include sample identifiers. Coordinates are PAF-standard 0-based half-open.
- `consensus_comparison.tsv` and `consensus_alignments.txt`: supplied consensus comparisons and readable pairwise alignments.
- `references_with_observed_variant.fasta`: the four expected references plus the explicitly labeled observed consensus.
- `versions.json`, `input_sha256.json`, and `sensitive_mapping_check.json`: software settings, input checksums, and mapping sensitivity results.

Run from the repository root:

```bash
python3 -B analysis/align_smjydd_references.py
```

The script uses project-local dependencies in `tmp/sequence_mapping_deps` (biopython 1.83, mappy 2.31, numpy 1.24.4). On this Apple Silicon Mac, installation required `ARCHFLAGS='-arch arm64'` when building mappy. The original exact-match analysis remains in `outputs/SMJYDD_target_sequence`.

Platform documentation: https://plasmidsaurus.com/faq/amplicon

Aligner documentation: https://github.com/lh3/minimap2 and https://github.com/lh3/minimap2/tree/master/python

# Sample-specific reference alignment

This analysis supersedes the unrestricted reference comparison for interpreting sample identities. Reference panels follow the user's clarification:

- L2: 328_ref and 346_MSI2 only.
- L4: 324_altUTR only.
- L8: 251_Exon1 only.

All 30,022 reads were remapped with minimap2/mappy 2.31, map-ont, against only their permitted reference(s). No observed-consensus reference was added. The original FASTQs were not changed.

## Are both 328 and 346 represented in L2?

Yes: the reads support both the insertion-absent (328-like) and insertion-containing (346-like) versions. Requiring sequence on both sides of the distinguishing region gives:

| Sample | 328-supporting reads | 346-supporting reads | 346 % among resolved reads | Other/conflicting | Not callable |
|---|---:|---:|---:|---:|---:|
| L2 HEK | 1,192 | 74 | 5.85% | 9 | 4,057 |
| L2 H295R | 555 | 4 | 0.72% | 5 | 4,863 |

Each of the 74 and four insertion-containing reads also prefers 346 in the whole-reference comparison. The four H295R reads are 876–887 nt long, align across 98.99–100% of the 346 reference, and have 98.65–99.33% alignment identity. Their IDs are 3163_SMJYDD_2, 3340_SMJYDD_2, 4012_SMJYDD_2, and 4124_SMJYDD_2. This is evidence for a minority insertion-containing read population, not evidence by itself of its biological origin or proof of an error-free 346 construct. In particular, the extracted H295R windows have a net 54-nt insertion versus 57 nt in the supplied 346 reference; small junction differences are allowed in these calls.

The percentages are among resolved, locus-spanning reads, not all reads or original cells/transcripts. Four supporting reads in H295R warrant more caution than 74 in HEK. Shared-region reads are not evidence of two versions. PCR bias and sample/sequencing artifacts cannot be distinguished from genuine biological mixtures by these counts alone.

The HEK count changed from 75 in the earlier comparison to 74 because the target-window comparison now uses only the supplied 346 insertion allele. One read becomes conflicting across its alignments instead of being called insertion-containing.

## L4 and L8

| Sample | Allowed reference | Reads meeting mapping criteria | Reads without qualifying alignment |
|---|---|---:|---:|
| L4 HEK | 324_altUTR | 0 | 4,524 |
| L4 H295R | 324_altUTR | 294 | 3,922 |
| L8 HEK | 251_Exon1 | 4,754 | 722 |
| L8 H295R | 251_Exon1 | 2,695 | 2,352 |

L4 H295R has 30 reads supporting the target-containing window and 237 supporting the deleted window when mapped directly to 324_altUTR; 3,949 reads are not callable at this locus. This describes a deletion relative to the expected L4 sequence, not an assignment to 328. The supplied consensus comparison independently of reference selection is consistent with that deletion (but is derived from the same sequencing data, not an independent experiment).

L4 HEK remains inconclusive under these alignment criteria; no match is not proof that its intended construct is absent. L8 is evaluated against 251 only, and the target-window assay is not applicable to that reference.

## Methods, validation, and files

Whole-reference mapping requires at least 80 aligned query bases and 85% alignment identity. For L2 a unique preference requires a score advantage of at least 20; otherwise the read is ambiguous. Whole-reference preferences can include partially aligned reads and are saved separately from the more conservative locus-spanning counts above.

The locus assay uses 40 bases on either side of the insertion, with at least 85% edit-distance identity in each outer 30-base flank. L2 is evaluated on 328 coordinates against the 328 and 346 windows only. L4 is evaluated directly on 324 coordinates against its intact window or deletion window. Present calls allow at most 14 edits in the full insertion-containing window; absent calls allow at most eight edits in the deletion window; the winning class must improve edit distance by at least ten. Conflicting calls across a read's alignments remain unresolved. At 90% flank identity, L2 present/absent counts are 75/1,179 and 4/549; one previously conflicting HEK read becomes callable when a poorer hit fails that stricter flank check.

Forward/reverse synthetic L4 intact/deletion checks pass. Every read is represented once in the assignment file, all counts reconcile, and every PAF alignment belongs to its sample's permitted panel. Dependencies and input sequences are the same as in the earlier analysis.

- `sample_summary.tsv`: sample-specific mapping and target-window counts.
- `read_assignments.tsv`: read-level calls, compatible references, alignment identity and coverage, and locus metrics.
- `alignments.paf`: alignments only to permitted references, with CIGAR and cs tags.
- `reference_panels.json`: the sample-to-reference mapping.

Reproduce from the project root with `python3 -B analysis/align_smjydd_sample_specific.py`. This uses the shared alignment helpers in `analysis/align_smjydd_references.py` and dependencies in `tmp/sequence_mapping_deps`.

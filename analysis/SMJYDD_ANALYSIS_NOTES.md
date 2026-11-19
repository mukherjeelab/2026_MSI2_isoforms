# Standing analysis settings for SMJYDD

User instruction: exclude every read shorter than 500 nucleotides from all subsequent calculations. Keep reads exactly 500 nt. This applies to means, denominators, mapping summaries, and variant analyses. Older output directories contain historical unfiltered results and must not be used as current summaries.

Use sample-specific references only:
- L2: 328_ref and 346_MSI2.
- L4: 324_altUTR.
- L8: 251_Exon1.

Current entry point: `python3 -B analysis/summarize_smjydd_min500.py`.
Current output directory: `outputs/SMJYDD_min500`.

Distinguish a read mapping to a reference from support for its intact sequence. Report recurrent deletions/other differences separately. All-read percentages use the retained >=500-nt reads; locus/variant-specific percentages must explicitly name their informative-read denominators. Recurring sequence variants are candidates for further validation, not proof of biological origin.

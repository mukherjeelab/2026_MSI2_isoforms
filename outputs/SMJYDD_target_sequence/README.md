# Exact target read counts

Target (57 nt): `TCAGTCCTGAATAGCTACAGTGCTCAACCGAATTTTGGCGCGCCCGCTTCCCCGGCA`

Reverse complement: `TGCCGGGGAAGCGGGCGCGCCAAAATTCGGTTGAGCACTGTAGCTATTCAGGACTGA`

Each FASTQ record is one read. A match requires the full contiguous target or its reverse complement anywhere in the read. Reads matching both orientations count once in matching_reads. No quality filtering, trimming, or deduplication is applied. Nonmatching means no exact full match; it does not establish biological absence, because reads may not span the target or may contain sequencing errors.

Run `python3 analysis/count_fastq_target.py` from the repository root. counts.tsv contains sample summaries; read_classifications.tsv contains every read's classification.

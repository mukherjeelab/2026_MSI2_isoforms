"""Count reads with an exact full target sequence in either orientation.

Run from the repository root: python3 analysis/count_fastq_target.py
No trimming, quality filtering, mismatch allowance, or deduplication is applied.
"""
import csv
from pathlib import Path

TARGET = "TCAGTCCTGAATAGCTACAGTGCTCAACCGAATTTTGGCGCGCCCGCTTCCCCGGCA"
REVERSE = TARGET.translate(str.maketrans("ACGT", "TGCA"))[::-1]
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "SMJYDD_target_sequence"


def records(path):
    """Read these four-line FASTQs, failing on malformed records."""
    with path.open() as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().rstrip("\r\n").upper()
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+") or not seq or len(seq) != len(quality):
                raise ValueError(f"Malformed FASTQ record in {path}: {header.strip()}")
            yield header[1:].strip(), seq


def main():
    paths = sorted((ROOT / "data" / "SMJYDD_fastq").glob("*.fastq"))
    if not paths:
        raise ValueError("No FASTQ files found")
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    with (OUT / "read_classifications.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "read_id", "read_length", "exact_forward", "exact_reverse_complement", "exact_either_orientation"])
        for path in paths:
            total = forward = reverse = either = both = short = 0
            for name, seq in records(path):
                f, r = TARGET in seq, REVERSE in seq
                total += 1
                forward += f
                reverse += r
                either += f or r
                both += f and r
                short += len(seq) < len(TARGET)
                writer.writerow([path.stem, name, len(seq), int(f), int(r), int(f or r)])
            rows.append([path.stem, total, forward, reverse, both, either, total - either, short, f"{100 * either / total:.4f}" if total else "NA"])
    fields = ["sample", "total_reads", "exact_forward", "exact_reverse_complement", "both_orientations", "matching_reads", "nonmatching_reads", "reads_shorter_than_target", "matching_percent"]
    totals = [sum(row[i] for row in rows) for i in range(1, 8)]
    rows.append(["TOTAL", *totals, f"{100 * totals[4] / totals[0]:.4f}" if totals[0] else "NA"])
    with (OUT / "counts.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        writer.writerows(rows)
    (OUT / "README.md").write_text(
        f"# Exact target read counts\n\nTarget ({len(TARGET)} nt): `{TARGET}`\n\n"
        f"Reverse complement: `{REVERSE}`\n\n"
        "Each FASTQ record is one read. A match requires the full contiguous target or its reverse complement anywhere in the read. "
        "Reads matching both orientations count once in matching_reads. No quality filtering, trimming, or deduplication is applied. "
        "Nonmatching means no exact full match; it does not establish biological absence, because reads may not span the target or may contain sequencing errors.\n\n"
        "Run `python3 analysis/count_fastq_target.py` from the repository root. "
        "counts.tsv contains sample summaries; read_classifications.tsv contains every read's classification.\n"
    )
    print("\t".join(fields))
    for row in rows:
        print("\t".join(map(str, row)))


if __name__ == "__main__":
    main()

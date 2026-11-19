#!/usr/bin/env python3
"""Convert featureCounts output into a sample-ordered raw count matrix."""

import argparse
import csv
from pathlib import Path


BAM_SUFFIX = ".mapped.sorted.bam"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_samples(path):
    with open(path, newline="") as handle:
        rows = (
            line for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        )
        reader = csv.DictReader(rows, delimiter="\t")
        if reader.fieldnames is None or "sample" not in reader.fieldnames:
            raise ValueError(f"{path} must contain a sample column")
        samples = [row["sample"].strip() for row in reader]
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("Metadata samples must be present and unique")
    return samples


def sample_from_bam_column(value):
    name = Path(value).name
    if not name.endswith(BAM_SUFFIX):
        raise ValueError(
            f"featureCounts column {value!r} does not end with {BAM_SUFFIX}"
        )
    return name[: -len(BAM_SUFFIX)]


def main():
    args = parse_args()
    samples = read_samples(args.metadata)

    with open(args.counts, newline="") as source:
        rows = (line for line in source if not line.startswith("#"))
        reader = csv.reader(rows, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 7 or header[:6] != [
            "Geneid", "Chr", "Start", "End", "Strand", "Length"
        ]:
            raise ValueError("Unexpected featureCounts header")
        observed = [sample_from_bam_column(value) for value in header[6:]]
        if len(observed) != len(set(observed)) or set(observed) != set(samples):
            raise ValueError(
                f"featureCounts samples {observed} do not match metadata {samples}"
            )
        indices = [6 + observed.index(sample) for sample in samples]

        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="") as destination:
            writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
            writer.writerow(["peak_id"] + samples)
            row_count = 0
            for row in reader:
                if len(row) != len(header):
                    raise ValueError("Inconsistent featureCounts row width")
                selected = [row[index] for index in indices]
                if any(not value.isdigit() for value in selected):
                    raise ValueError(f"Non-integer count for peak {row[0]}")
                writer.writerow([row[0]] + selected)
                row_count += 1

    print(f"Wrote {row_count} peaks and {len(samples)} samples to {output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Combine per-BAM primary mapped-alignment counts in sample-sheet order."""

import argparse
import csv
from pathlib import Path


SUFFIX = ".primary_mapped.txt"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--counts", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def metadata_samples(path):
    with open(path, newline="") as handle:
        rows = (
            line for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        )
        reader = csv.DictReader(rows, delimiter="\t")
        if reader.fieldnames is None or "sample" not in reader.fieldnames:
            raise ValueError(f"{path} must contain a sample column")
        samples = [row["sample"].strip() for row in reader]
    if not samples:
        raise ValueError(f"No samples found in {path}")
    return samples


def main():
    args = parse_args()
    samples = metadata_samples(args.metadata)
    paths = {}
    for value in args.counts:
        path = Path(value)
        if not path.name.endswith(SUFFIX):
            raise ValueError(f"Unexpected library-size filename: {path}")
        sample = path.name[: -len(SUFFIX)]
        if sample in paths:
            raise ValueError(f"Duplicate library-size file for {sample}")
        paths[sample] = path

    if set(paths) != set(samples):
        raise ValueError(
            "Library-size files and metadata differ: "
            f"files={sorted(paths)}, metadata={sorted(samples)}"
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "primary_mapped_alignments"])
        for sample in samples:
            text = paths[sample].read_text().strip()
            try:
                count = int(text)
            except ValueError as error:
                raise ValueError(
                    f"Invalid mapped-alignment count for {sample}: {text!r}"
                ) from error
            if count <= 0:
                raise ValueError(f"Mapped-alignment count for {sample} is {count}")
            writer.writerow([sample, count])

    print(f"Wrote full-library sizes for {len(samples)} samples to {output}")


if __name__ == "__main__":
    main()

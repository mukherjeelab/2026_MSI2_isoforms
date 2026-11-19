#!/usr/bin/env python3
"""Filter and strand-specifically merge a BED file while preserving provenance."""

import argparse
import csv
import re
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple


@dataclass(frozen=True)
class Peak:
    chrom: str
    start: int
    end: int
    name: str
    score: float
    score_text: str
    strand: str
    summit_start: str
    summit_end: str
    gene_id: str
    line_number: int
    original_fields: Tuple[str, ...]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--minimum-score", required=True, type=float)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--members", required=True)
    parser.add_argument("--saf", required=True)
    return parser.parse_args()


def chromosome_key(chrom):
    return [int(token) if token.isdigit() else token.lower()
            for token in re.split(r"(\d+)", chrom)]


def read_peaks(path, minimum_score):
    retained = []
    total = 0
    with open(path, newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            total += 1
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(
                    f"{path}:{line_number} has {len(fields)} columns; at least 6 are required"
                )
            try:
                start = int(fields[1])
                end = int(fields[2])
                score = float(fields[4])
            except ValueError as error:
                raise ValueError(f"Invalid coordinate or score at {path}:{line_number}") from error
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval at {path}:{line_number}")
            strand = fields[5]
            if strand not in {"+", "-", "."}:
                raise ValueError(f"Invalid strand {strand!r} at {path}:{line_number}")
            if score < minimum_score:
                continue
            retained.append(
                Peak(
                    chrom=fields[0],
                    start=start,
                    end=end,
                    name=fields[3] or f"source_peak_{line_number}",
                    score=score,
                    score_text=fields[4],
                    strand=strand,
                    summit_start=fields[6] if len(fields) > 6 else "",
                    summit_end=fields[7] if len(fields) > 7 else "",
                    gene_id=fields[8] if len(fields) > 8 else "",
                    line_number=line_number,
                    original_fields=tuple(fields),
                )
            )
    return total, retained


def merge_peaks(peaks):
    ordered = sorted(
        peaks,
        key=lambda peak: (
            chromosome_key(peak.chrom), peak.strand, peak.start, peak.end, peak.name
        ),
    )
    groups = []
    current = []
    current_end = -1
    current_key = None

    for peak in ordered:
        key = (peak.chrom, peak.strand)
        # BED intervals overlap when the next start is strictly less than the
        # current end. Directly adjacent intervals are not merged.
        if current and key == current_key and peak.start < current_end:
            current.append(peak)
            current_end = max(current_end, peak.end)
        else:
            if current:
                groups.append(current)
            current = [peak]
            current_end = peak.end
            current_key = key
    if current:
        groups.append(current)
    return groups


def format_number(value):
    return f"{value:.15g}"


def main():
    args = parse_args()
    total, retained = read_peaks(args.input, args.minimum_score)
    if not retained:
        raise ValueError(f"No peaks have score >= {args.minimum_score}")
    groups = merge_peaks(retained)

    for path in (args.bed, args.metadata, args.members, args.saf):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    with ExitStack() as stack:
        bed_handle = stack.enter_context(open(args.bed, "w", newline=""))
        metadata_handle = stack.enter_context(open(args.metadata, "w", newline=""))
        members_handle = stack.enter_context(open(args.members, "w", newline=""))
        saf_handle = stack.enter_context(open(args.saf, "w", newline=""))
        bed = csv.writer(bed_handle, delimiter="\t", lineterminator="\n")
        metadata = csv.writer(metadata_handle, delimiter="\t", lineterminator="\n")
        members = csv.writer(members_handle, delimiter="\t", lineterminator="\n")
        saf = csv.writer(saf_handle, delimiter="\t", lineterminator="\n")

        metadata.writerow([
            "peak_id", "chromosome", "start_0based", "end_0based_exclusive",
            "strand", "length", "source_peak_count", "peak_score_max",
            "peak_score_mean", "peak_score_sum", "gene_ids"
        ])
        members.writerow([
            "peak_id", "source_peak_id", "source_score", "chromosome",
            "start_0based", "end_0based_exclusive", "strand", "summit_start",
            "summit_end", "gene_id", "source_line", "original_bed_fields"
        ])
        saf.writerow(["GeneID", "Chr", "Start", "End", "Strand"])

        for index, group in enumerate(groups, start=1):
            peak_id = f"MSI2_peak_{index:06d}"
            chrom = group[0].chrom
            strand = group[0].strand
            start = min(peak.start for peak in group)
            end = max(peak.end for peak in group)
            scores = [peak.score for peak in group]
            score_max = max(scores)
            score_sum = sum(scores)
            score_mean = score_sum / len(scores)
            gene_ids = ";".join(sorted({peak.gene_id for peak in group if peak.gene_id}))

            bed.writerow([chrom, start, end, peak_id, format_number(score_max), strand])
            metadata.writerow([
                peak_id, chrom, start, end, strand, end - start, len(group),
                format_number(score_max), format_number(score_mean),
                format_number(score_sum), gene_ids
            ])
            # SAF uses one-based, inclusive coordinates; BED uses zero-based,
            # half-open coordinates.
            saf.writerow([peak_id, chrom, start + 1, end, strand])

            for peak in group:
                members.writerow([
                    peak_id, peak.name, peak.score_text, peak.chrom, peak.start,
                    peak.end, peak.strand, peak.summit_start, peak.summit_end,
                    peak.gene_id, peak.line_number, "|".join(peak.original_fields)
                ])

    print(f"Input peaks: {total}")
    print(f"Peaks with score >= {args.minimum_score}: {len(retained)}")
    print(f"Merged strand-specific peak groups: {len(groups)}")
    print(f"Peaks absorbed into overlapping groups: {len(retained) - len(groups)}")


if __name__ == "__main__":
    main()

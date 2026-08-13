#!/usr/bin/env python3
"""Build task 1 marker matrix for genomic prediction.

Input depth files are sample-major. This script extracts the marker set that
already passed unsupervised call-rate/MAF filters in the GWAS baseline and
writes a row-major float32 matrix:

    rows = samples
    columns = markers
    value = alternate-read fraction from depth column 6

Missing marker calls are stored as NaN. Downstream R code imputes missing calls
to the marker mean after centering.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from array import array
from pathlib import Path


TRAITS = ("GZZTF", "TCD", "ZZZTF","HYTF","ZGD","SCD","ZLCD","ZXWF")
MISSING = {"", ".", "NA", "NaN", "nan", "N/A", "DIV/0!"}


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except Exception:
        return False


def parse_alt_fraction(depth_field: str) -> float | None:
    if depth_field in MISSING:
        return None
    try:
        counts = [float(x) for x in depth_field.split(",") if x != ""]
    except Exception:
        return None
    if len(counts) < 2:
        return None
    total = sum(counts)
    if total <= 0:
        return None
    return sum(counts[1:]) / total


def read_markers(marker_source: Path, max_markers: int) -> list[dict[str, str]]:
    opener = gzip.open if marker_source.suffix == ".gz" else open
    markers: list[dict[str, str]] = []
    with opener(marker_source, "rt", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            markers.append(row)
    if max_markers and len(markers) > max_markers:
        if max_markers == 1:
            markers = [markers[0]]
        else:
            idxs = sorted(
                set(round(i * (len(markers) - 1) / (max_markers - 1)) for i in range(max_markers))
            )
            markers = [markers[i] for i in idxs]
    return markers


def read_samples(sample_qc: Path) -> list[dict[str, str]]:
    samples: list[dict[str, str]] = []
    with sample_qc.open(errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["depth_empty"] != "0" or row["has_depth"] != "1":
                continue
            if not all(row[trait] not in MISSING and is_number(row[trait]) for trait in TRAITS):
                continue
            samples.append(row)
    return samples


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--marker-source", required=True, type=Path)
    parser.add_argument("--tag", default="task1_gp_full")
    parser.add_argument("--max-markers", type=int, default=0)
    parser.add_argument("--max-lines-per-file", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=10)
    args = parser.parse_args()

    metadata_dir = args.work_dir / "metadata"
    tmp_dir = args.work_dir / "tmp" / "task1_genomic_prediction"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    sample_qc = metadata_dir / "task1_sample_qc.tsv"
    samples = read_samples(sample_qc)
    markers = read_markers(args.marker_source, args.max_markers)

    marker_to_idx = {row["marker"]: i for i, row in enumerate(markers)}
    n_samples = len(samples)
    n_markers = len(markers)
    if n_samples == 0 or n_markers == 0:
        raise SystemExit("No samples or markers selected.")

    matrix_path = tmp_dir / f"{args.tag}.matrix.float32.bin"
    markers_path = metadata_dir / f"{args.tag}.markers.tsv"
    samples_path = metadata_dir / f"{args.tag}.samples.tsv"
    manifest_path = metadata_dir / f"{args.tag}.manifest.tsv"

    with markers_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["marker", "chr", "pos", "ref", "alt", "call_rate", "mean_alt_fraction", "maf_like"])
        for row in markers:
            writer.writerow(
                [
                    row.get("marker", ""),
                    row.get("chr", ""),
                    row.get("pos", ""),
                    row.get("ref", ""),
                    row.get("alt", ""),
                    row.get("call_rate", ""),
                    row.get("mean_alt_fraction", ""),
                    row.get("maf_like", ""),
                ]
            )

    with samples_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "depth_raw_id", "depth_path", *TRAITS])
        for row in samples:
            writer.writerow([row["sample_id"], row["depth_raw_id"], row["depth_path"], *(row[t] for t in TRAITS)])

    matrix_path.unlink(missing_ok=True)
    total_observed = 0
    with matrix_path.open("wb") as out:
        for sample_i, sample in enumerate(samples, 1):
            values = array("f", [math.nan]) * n_markers
            observed = 0
            depth_path = Path(sample["depth_path"])
            with depth_path.open(errors="replace") as handle:
                for line_i, line in enumerate(handle, 1):
                    if args.max_lines_per_file and line_i > args.max_lines_per_file:
                        break
                    cols = line.rstrip("\n\r").split("\t")
                    if len(cols) < 7:
                        continue
                    key = f"{cols[0]}:{cols[1]}:{cols[2]}:{cols[3]}"
                    idx = marker_to_idx.get(key)
                    if idx is None:
                        continue
                    x = parse_alt_fraction(cols[6])
                    if x is None:
                        continue
                    if math.isnan(values[idx]):
                        observed += 1
                    values[idx] = x
            values.tofile(out)
            total_observed += observed
            if sample_i % args.progress_every == 0 or sample_i == n_samples:
                print(
                    f"processed_samples={sample_i}/{n_samples} "
                    f"n_markers={n_markers} observed_calls={total_observed}",
                    flush=True,
                )

    with manifest_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["key", "value"])
        writer.writerow(["tag", args.tag])
        writer.writerow(["n_samples", n_samples])
        writer.writerow(["n_markers", n_markers])
        writer.writerow(["matrix_file", str(matrix_path)])
        writer.writerow(["samples_file", str(samples_path)])
        writer.writerow(["markers_file", str(markers_path)])
        writer.writerow(["marker_source", str(args.marker_source)])
        writer.writerow(["max_markers", args.max_markers])
        writer.writerow(["max_lines_per_file", args.max_lines_per_file])
        writer.writerow(["total_observed_calls", total_observed])

    print(f"Wrote {matrix_path}")
    print(f"Wrote {samples_path}")
    print(f"Wrote {markers_path}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()

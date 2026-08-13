#!/usr/bin/env python3
"""Prepare sample metadata for task 1 G2P baselines.

Task 1:
  X1: 183-sample/*-depth.txt
  Y1: 183-sample/{GZZTF,TCD,ZZZTF}.trait

This script only performs bookkeeping/QC. It does not modify the raw data.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


MISSING = {"", ".", "NA", "NaN", "nan", "N/A", "DIV/0!"}
TRAITS = ("GZZTF", "TCD", "ZZZTF","HYTF","ZGD","SCD","ZLCD","ZXWF")


def norm_sample_id(sample: str) -> str:
    """Normalize the known '-'/'_' naming difference in depth files."""
    return sample.replace("_", "-")


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except Exception:
        return False


def read_trait(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            cols = line.strip().split()
            if not cols:
                continue
            sample = cols[0]
            value = cols[2] if len(cols) > 2 else ""
            values[sample] = value
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    args = parser.parse_args()

    x_dir = args.data_dir / "183-sample"
    metadata_dir = args.work_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    depth_by_norm: dict[str, tuple[str, Path, int, bool]] = {}
    duplicate_norm_ids: dict[str, list[str]] = {}
    for path in sorted(x_dir.glob("*-depth.txt")):
        raw_id = path.name[: -len("-depth.txt")]
        sample_id = norm_sample_id(raw_id)
        size = path.stat().st_size
        is_empty = size == 0
        if sample_id in depth_by_norm:
            duplicate_norm_ids.setdefault(sample_id, []).append(raw_id)
        depth_by_norm[sample_id] = (raw_id, path, size, is_empty)

    traits = {trait: read_trait(x_dir / f"{trait}.trait") for trait in TRAITS}
    all_samples = sorted(set().union(*(set(values) for values in traits.values())))

    sample_qc_path = metadata_dir / "task1_sample_qc.tsv"
    with sample_qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "depth_raw_id",
                "depth_path",
                "depth_size_bytes",
                "depth_empty",
                "has_depth",
                "GZZTF",
                "TCD",
                "ZZZTF",
                "HYTF",
                "ZGD",
                "SCD",
                "ZLCD",
                "ZXWF"
            ]
        )
        for sample in all_samples:
            raw_id, depth_path, size, is_empty = depth_by_norm.get(
                sample, ("", Path(""), 0, True)
            )
            writer.writerow(
                [
                    sample,
                    raw_id,
                    str(depth_path),
                    size,
                    int(is_empty),
                    int(sample in depth_by_norm),
                    traits["GZZTF"].get(sample, ""),
                    traits["TCD"].get(sample, ""),
                    traits["ZZZTF"].get(sample, ""),
                    traits["HYTF"].get(sample, ""),
                    traits["ZGD"].get(sample, ""),
                    traits["SCD"].get(sample, ""),
                    traits["ZLCD"].get(sample, ""),
                    traits["ZXWF"].get(sample, ""),
                ]
            )

    summary_path = metadata_dir / "task1_trait_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "trait",
                "y_samples",
                "numeric_y",
                "missing_y",
                "matched_depth_after_name_normalization",
                "nonempty_depth_and_numeric_y",
            ]
        )
        for trait, values in traits.items():
            y_samples = set(values)
            numeric = {s for s, v in values.items() if v not in MISSING and is_number(v)}
            missing = y_samples - numeric
            matched = {s for s in y_samples if s in depth_by_norm}
            usable = {
                s
                for s in numeric
                if s in depth_by_norm and depth_by_norm[s][2] > 0
            }
            writer.writerow(
                [
                    trait,
                    len(y_samples),
                    len(numeric),
                    len(missing),
                    len(matched),
                    len(usable),
                ]
            )

    if duplicate_norm_ids:
        dup_path = metadata_dir / "task1_duplicate_normalized_depth_ids.tsv"
        with dup_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["normalized_sample_id", "extra_raw_depth_ids"])
            for sample, raw_ids in sorted(duplicate_norm_ids.items()):
                writer.writerow([sample, ",".join(raw_ids)])

    print(f"Wrote {sample_qc_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()

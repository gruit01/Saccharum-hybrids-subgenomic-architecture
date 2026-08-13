#!/usr/bin/env python3
"""Dependency-free GWAS-like baseline for task 1 depth files.

The raw X files are sample-major depth tables, not a ready PLINK matrix. For
each marker we encode the sample value as alternate-read fraction from column 6:

    alt_fraction = sum(alt_depths) / sum(ref_depth + alt_depths)

Then we run per-marker univariate linear regression:

    phenotype ~ alt_fraction

The implementation uses only Python's standard library because the server does
not currently provide numpy/scipy/sklearn/R/GEMMA/PLINK.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import math
import os
import time
from array import array
from pathlib import Path


TRAITS = ("GZZTF", "TCD", "ZZZTF","HYTF","ZGD","SCD","ZLCD","ZXWF")
MISSING = {"", ".", "NA", "NaN", "nan", "N/A", "DIV/0!"}


def norm_sample_id(sample: str) -> str:
    return sample.replace("_", "-")


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except Exception:
        return False


def read_trait(path: Path) -> dict[str, float | None]:
    values: dict[str, float | None] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            cols = line.strip().split()
            if not cols:
                continue
            value = cols[2] if len(cols) > 2 else ""
            values[cols[0]] = float(value) if value not in MISSING and is_number(value) else None
    return values


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


def betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for incomplete beta, Numerical Recipes style."""
    max_iter = 200
    eps = 3.0e-14
    fpmin = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betai(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * betacf(a, b, x) / a
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b


def t_two_sided_pvalue(t_value: float, df: int) -> float:
    if df <= 0 or not math.isfinite(t_value):
        return 1.0
    t_abs = abs(t_value)
    x = df / (df + t_abs * t_abs)
    return max(0.0, min(1.0, betai(0.5 * df, 0.5, x)))


def marker_key(cols: list[str]) -> str:
    return f"{cols[0]}:{cols[1]}:{cols[2]}:{cols[3]}"


def load_inputs(data_dir: Path) -> tuple[list[tuple[str, Path]], dict[str, dict[str, float | None]]]:
    x_dir = data_dir / "183-sample"
    depth_files = []
    for path in sorted(x_dir.glob("*-depth.txt")):
        raw_id = path.name[: -len("-depth.txt")]
        sample_id = norm_sample_id(raw_id)
        if path.stat().st_size == 0:
            continue
        depth_files.append((sample_id, path))
    traits = {trait: read_trait(x_dir / f"{trait}.trait") for trait in TRAITS}
    return depth_files, traits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--min-call-rate", type=float, default=0.80)
    parser.add_argument("--min-maf", type=float, default=0.01)
    parser.add_argument("--top-n", type=int, default=1000)
    parser.add_argument("--max-lines-per-file", type=int, default=0,
                        help="Smoke-test mode: only read this many lines per depth file.")
    parser.add_argument("--progress-every", type=int, default=10)
    args = parser.parse_args()

    result_dir = args.work_dir / "results" / "task1_gwas"
    tmp_dir = args.work_dir / "tmp" / "task1_gwas"
    metadata_dir = args.work_dir / "metadata"
    result_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    depth_files, traits = load_inputs(args.data_dir)
    usable_samples = []
    for sample_id, path in depth_files:
        if any(traits[trait].get(sample_id) is not None for trait in TRAITS):
            usable_samples.append((sample_id, path))

    n_total = len(usable_samples)
    min_n = max(3, int(math.ceil(args.min_call_rate * n_total)))
    print(f"Usable non-empty X files with at least one Y: {n_total}", flush=True)
    print(f"Minimum marker call count: {min_n} ({args.min_call_rate:.3f})", flush=True)

    keys: list[str] = []
    key_to_idx: dict[str, int] = {}
    n = array("I")
    sx = array("d")
    sx2 = array("d")
    tn = {trait: array("I") for trait in TRAITS}
    sy = {trait: array("d") for trait in TRAITS}
    sy2 = {trait: array("d") for trait in TRAITS}
    sxy = {trait: array("d") for trait in TRAITS}

    def add_marker(key: str) -> int:
        idx = key_to_idx.get(key)
        if idx is not None:
            return idx
        idx = len(keys)
        key_to_idx[key] = idx
        keys.append(key)
        n.append(0)
        sx.append(0.0)
        sx2.append(0.0)
        for trait in TRAITS:
            tn[trait].append(0)
            sy[trait].append(0.0)
            sy2[trait].append(0.0)
            sxy[trait].append(0.0)
        return idx

    started = time.time()
    processed_lines = 0
    for sample_i, (sample_id, path) in enumerate(usable_samples, 1):
        sample_y = {trait: traits[trait].get(sample_id) for trait in TRAITS}
        line_count = 0
        with path.open(errors="replace") as handle:
            for line in handle:
                if args.max_lines_per_file and line_count >= args.max_lines_per_file:
                    break
                cols = line.rstrip("\n\r").split("\t")
                if len(cols) < 6:
                    continue
                x = parse_alt_fraction(cols[5])
                if x is None:
                    continue
                idx = add_marker(marker_key(cols))
                n[idx] += 1
                sx[idx] += x
                sx2[idx] += x * x
                for trait, y in sample_y.items():
                    if y is None:
                        continue
                    tn[trait][idx] += 1
                    sy[trait][idx] += y
                    sy2[trait][idx] += y * y
                    sxy[trait][idx] += x * y
                line_count += 1
                processed_lines += 1
        if sample_i % args.progress_every == 0 or sample_i == len(usable_samples):
            elapsed = time.time() - started
            print(
                f"processed_samples={sample_i}/{len(usable_samples)} "
                f"unique_markers={len(keys)} processed_lines={processed_lines} "
                f"elapsed_sec={elapsed:.1f}",
                flush=True,
            )

    run_tag = "smoke" if args.max_lines_per_file else "full"
    summary_rows = []
    for trait in TRAITS:
        all_path = result_dir / f"task1_{trait}_gwas_{run_tag}.all_markers.tsv.gz"
        top_path = result_dir / f"task1_{trait}_gwas_{run_tag}.top{args.top_n}.tsv"
        sig_path = result_dir / f"task1_{trait}_gwas_{run_tag}.significant.tsv"
        top_heap: list[tuple[float, tuple[str, int, float, float, float, float, float]]] = []
        tested = 0
        with gzip.open(all_path, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["marker", "chr", "pos", "ref", "alt", "n", "call_rate", "mean_alt_fraction", "maf_like", "beta", "r2", "t", "p"])
            for idx, key in enumerate(keys):
                nn = tn[trait][idx]
                if nn < min_n:
                    continue
                x_sum = sx[idx]
                x2_sum = sx2[idx]
                y_sum = sy[trait][idx]
                y2_sum = sy2[trait][idx]
                xy_sum = sxy[trait][idx]
                den_x = nn * x2_sum - x_sum * x_sum
                den_y = nn * y2_sum - y_sum * y_sum
                if den_x <= 0.0 or den_y <= 0.0:
                    continue
                mean_x = x_sum / nn
                maf_like = min(mean_x, 1.0 - mean_x)
                if maf_like < args.min_maf:
                    continue
                r_num = nn * xy_sum - x_sum * y_sum
                r2 = (r_num * r_num) / (den_x * den_y)
                if r2 >= 1.0:
                    r2 = 0.999999999999
                beta = r_num / den_x
                t_value = math.copysign(math.sqrt((nn - 2) * r2 / max(1.0 - r2, 1e-300)), r_num)
                p_value = t_two_sided_pvalue(t_value, nn - 2)
                chrom, pos, ref, alt = key.split(":", 3)
                tested += 1
                row = (key, nn, nn / n_total, mean_x, maf_like, beta, r2, t_value, p_value)
                writer.writerow([key, chrom, pos, ref, alt, nn, f"{nn / n_total:.6g}", f"{mean_x:.6g}", f"{maf_like:.6g}", f"{beta:.8g}", f"{r2:.8g}", f"{t_value:.8g}", f"{p_value:.8g}"])
                heap_item = (-p_value, row)
                if len(top_heap) < args.top_n:
                    heapq.heappush(top_heap, heap_item)
                elif heap_item > top_heap[0]:
                    heapq.heapreplace(top_heap, heap_item)

        bonferroni = 0.05 / tested if tested else 0.0
        top_rows = sorted((item[1] for item in top_heap), key=lambda r: r[-1])
        with top_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["marker", "chr", "pos", "ref", "alt", "n", "call_rate", "mean_alt_fraction", "maf_like", "beta", "r2", "t", "p"])
            for row in top_rows:
                key, nn, call_rate, mean_x, maf_like, beta, r2, t_value, p_value = row
                chrom, pos, ref, alt = key.split(":", 3)
                writer.writerow([key, chrom, pos, ref, alt, nn, f"{call_rate:.6g}", f"{mean_x:.6g}", f"{maf_like:.6g}", f"{beta:.8g}", f"{r2:.8g}", f"{t_value:.8g}", f"{p_value:.8g}"])

        sig_count = 0
        with gzip.open(all_path, "rt") as src, sig_path.open("w", newline="") as dst:
            reader = csv.reader(src, delimiter="\t")
            writer = csv.writer(dst, delimiter="\t")
            header = next(reader)
            writer.writerow(header + ["bonferroni_alpha_0.05"])
            for row in reader:
                p_value = float(row[-1])
                if p_value <= bonferroni:
                    sig_count += 1
                    writer.writerow(row + [f"{bonferroni:.8g}"])

        summary_rows.append((trait, tested, bonferroni, sig_count, str(all_path), str(top_path), str(sig_path)))
        print(f"{trait}: tested={tested} bonferroni={bonferroni:.3g} significant={sig_count}", flush=True)

    summary_path = metadata_dir / f"task1_gwas_{run_tag}_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["trait", "tested_markers", "bonferroni_alpha_0.05", "significant_markers", "all_markers_file", "top_snps_file", "significant_file"])
        writer.writerows(summary_rows)
    print(f"Wrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()

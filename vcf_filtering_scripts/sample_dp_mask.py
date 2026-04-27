#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from numbers import Real
from pathlib import Path

try:
    import pysam
except ModuleNotFoundError:
    pysam = None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute sample-specific DP thresholds and mask GT values outside range. "
            "Supports single-pass orchestration plus region-parallel compute/aggregate/mask subcommands."
        )
    )
    subparsers = parser.add_subparsers(dest="command")

    full_parser = subparsers.add_parser("full", help="Legacy two-pass whole-file workflow")
    add_full_args(full_parser)

    stats_parser = subparsers.add_parser(
        "compute-stats",
        help="Scan an input VCF/BCF or one indexed region and write per-sample DP sum/count",
    )
    add_stats_args(stats_parser)

    aggregate_parser = subparsers.add_parser(
        "aggregate-stats",
        help="Combine one or more compute-stats outputs into final sample thresholds",
    )
    add_aggregate_args(aggregate_parser)

    mask_parser = subparsers.add_parser(
        "mask",
        help="Mask GT outside precomputed per-sample DP thresholds, optionally on one region",
    )
    add_mask_args(mask_parser)

    return parser


def add_input_vcf_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input-vcf", required=True, help="Input VCF/BCF path")


def add_threads_arg(parser: argparse.ArgumentParser, help_text: str) -> None:
    parser.add_argument("--threads", type=int, default=1, help=help_text)


def add_region_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--region", help="Optional indexed region such as chr1 or chr1:1-5000000")


def add_full_args(parser: argparse.ArgumentParser) -> None:
    add_input_vcf_arg(parser)
    parser.add_argument("--output-bcf", required=True, help="Output masked BCF path")
    parser.add_argument("--sample-table", required=True, help="Output per-sample DP summary TSV")
    parser.add_argument("--min-factor", type=float, required=True, help="Min DP factor")
    parser.add_argument("--max-factor", type=float, required=True, help="Max DP factor")
    add_threads_arg(parser, "HTSlib reader/writer threads. Default: 1")


def add_stats_args(parser: argparse.ArgumentParser) -> None:
    add_input_vcf_arg(parser)
    parser.add_argument("--output-stats", required=True, help="Output per-sample DP sum/count TSV")
    add_region_arg(parser)
    add_threads_arg(parser, "HTSlib reader threads. Default: 1")


def add_aggregate_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--stats-file",
        action="append",
        default=[],
        help="One compute-stats TSV. Repeat for multiple inputs.",
    )
    parser.add_argument(
        "--stats-list-file",
        help="Optional file containing one compute-stats TSV path per line",
    )
    parser.add_argument("--sample-table", required=True, help="Output per-sample DP summary TSV")
    parser.add_argument("--min-factor", type=float, required=True, help="Min DP factor")
    parser.add_argument("--max-factor", type=float, required=True, help="Max DP factor")


def add_mask_args(parser: argparse.ArgumentParser) -> None:
    add_input_vcf_arg(parser)
    parser.add_argument("--output-bcf", required=True, help="Output masked BCF path")
    parser.add_argument(
        "--thresholds-table",
        required=True,
        help="Sample threshold TSV from full mode or aggregate-stats",
    )
    add_region_arg(parser)
    add_threads_arg(parser, "HTSlib reader/writer threads. Default: 1")


def parse_args() -> argparse.Namespace:
    parser = build_parser()
    argv = None
    import sys

    if len(sys.argv) > 1 and sys.argv[1] not in {"full", "compute-stats", "aggregate-stats", "mask", "-h", "--help"}:
        argv = ["full", *sys.argv[1:]]
    return parser.parse_args(argv)


def parse_dp_value(dp_value) -> float | None:
    if dp_value is None:
        return None
    if isinstance(dp_value, Real):
        if math.isnan(dp_value):
            return None
        return float(dp_value)
    if isinstance(dp_value, (tuple, list)):
        if not dp_value:
            return None
        first_value = dp_value[0]
        if isinstance(first_value, Real) and not math.isnan(first_value):
            return float(first_value)
    return None


def require_pysam() -> None:
    if pysam is None:
        raise ModuleNotFoundError(
            "pysam is required to run this script. Activate the environment that provides pysam."
        )


def iter_records(variant_file: pysam.VariantFile, region: str | None):
    if region:
        yield from variant_file.fetch(region=region)
    else:
        yield from variant_file


def compute_sample_dp_totals(
    input_vcf: str,
    threads: int,
    region: str | None = None,
):
    require_pysam()
    with pysam.VariantFile(input_vcf, threads=threads) as variant_file:
        sample_names = list(variant_file.header.samples)
        sample_count = len(sample_names)
        sums = [0.0] * sample_count
        counts = [0] * sample_count

        for record in iter_records(variant_file, region):
            samples = record.samples
            for i in range(sample_count):
                sample_data = samples[i]
                gt = sample_data.get("GT")
                if gt is None or all(allele is None for allele in gt):
                    continue
                dp = parse_dp_value(sample_data.get("DP"))
                if dp is None:
                    continue
                sums[i] += dp
                counts[i] += 1

    return sample_names, sums, counts


def build_thresholds(sums, counts, min_factor: float, max_factor: float):
    mean_dps = []
    min_dps = []
    max_dps = []
    for total, count in zip(sums, counts):
        mean_dp = total / count if count > 0 else 0.0
        mean_dps.append(mean_dp)
        min_dps.append(mean_dp * min_factor)
        max_dps.append(mean_dp * max_factor)
    return mean_dps, min_dps, max_dps


def write_stats_table(sample_names, sums, counts, output_stats: str) -> None:
    with open(output_stats, "w", encoding="utf-8") as out_handle:
        out_handle.write("sample\tdp_sum\tdp_count\n")
        for sample, dp_sum, dp_count in zip(sample_names, sums, counts):
            out_handle.write(f"{sample}\t{dp_sum:.6f}\t{dp_count}\n")


def read_stats_table(path: str):
    sample_names = []
    sums = []
    counts = []
    with open(path, "r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            sample, dp_sum, dp_count = line.rstrip("\n").split("\t")
            sample_names.append(sample)
            sums.append(float(dp_sum))
            counts.append(int(dp_count))
    return sample_names, sums, counts


def read_path_list(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def aggregate_stats_tables(stats_files: list[str]):
    merged_names = None
    merged_sums = None
    merged_counts = None

    for stats_file in stats_files:
        sample_names, sums, counts = read_stats_table(stats_file)
        if merged_names is None:
            merged_names = sample_names
            merged_sums = sums
            merged_counts = counts
            continue
        if sample_names != merged_names:
            raise ValueError(f"Sample order mismatch in stats file: {stats_file}")
        for i, (dp_sum, dp_count) in enumerate(zip(sums, counts)):
            merged_sums[i] += dp_sum
            merged_counts[i] += dp_count

    if merged_names is None:
        raise ValueError("No stats files provided")
    return merged_names, merged_sums, merged_counts


def write_threshold_table(
    sample_names,
    mean_dps,
    min_dps,
    max_dps,
    sample_table_path: str,
) -> None:
    header = "sample\tmean_dp\tmin_dp\tmax_dp\n"
    with open(sample_table_path, "w", encoding="utf-8") as sample_out:
        sample_out.write(header)
        for sample, mean_dp, min_dp, max_dp in zip(sample_names, mean_dps, min_dps, max_dps):
            row = f"{sample}\t{mean_dp:.6f}\t{min_dp:.6f}\t{max_dp:.6f}\n"
            sample_out.write(row)


def read_thresholds_table(path: str):
    sample_names = []
    min_dps = []
    max_dps = []
    with open(path, "r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            sample, _mean_dp, min_dp, max_dp = line.rstrip("\n").split("\t")
            sample_names.append(sample)
            min_dps.append(float(min_dp))
            max_dps.append(float(max_dp))
    return sample_names, min_dps, max_dps


def mask_variants(
    input_vcf: str,
    output_bcf: str,
    min_dps,
    max_dps,
    threads: int,
    region: str | None = None,
) -> None:
    require_pysam()
    Path(output_bcf).parent.mkdir(parents=True, exist_ok=True)
    with pysam.VariantFile(input_vcf, threads=threads) as in_vcf, pysam.VariantFile(
        output_bcf, "wb", header=in_vcf.header, threads=threads
    ) as out_bcf:
        sample_count = len(min_dps)
        for record in iter_records(in_vcf, region):
            samples = record.samples
            for i in range(sample_count):
                sample_data = samples[i]
                gt = sample_data.get("GT")
                if gt is None or all(allele is None for allele in gt):
                    continue
                dp = parse_dp_value(sample_data.get("DP"))
                if dp is None:
                    continue
                if dp < min_dps[i] or dp > max_dps[i]:
                    ploidy = len(gt) if gt else 2
                    sample_data["GT"] = tuple(None for _ in range(ploidy))
            out_bcf.write(record)


def run_full(args: argparse.Namespace) -> int:
    sample_names, sums, counts = compute_sample_dp_totals(
        args.input_vcf,
        args.threads,
    )
    mean_dps, min_dps, max_dps = build_thresholds(
        sums,
        counts,
        args.min_factor,
        args.max_factor,
    )
    write_threshold_table(
        sample_names,
        mean_dps,
        min_dps,
        max_dps,
        args.sample_table,
    )
    mask_variants(
        args.input_vcf,
        args.output_bcf,
        min_dps,
        max_dps,
        args.threads,
    )
    return 0


def run_compute_stats(args: argparse.Namespace) -> int:
    sample_names, sums, counts = compute_sample_dp_totals(
        args.input_vcf,
        args.threads,
        args.region,
    )
    write_stats_table(sample_names, sums, counts, args.output_stats)
    return 0


def run_aggregate_stats(args: argparse.Namespace) -> int:
    stats_files = list(args.stats_file)
    if args.stats_list_file:
        stats_files.extend(read_path_list(args.stats_list_file))
    if not stats_files:
        raise ValueError("Provide at least one --stats-file or --stats-list-file")
    sample_names, sums, counts = aggregate_stats_tables(stats_files)
    mean_dps, min_dps, max_dps = build_thresholds(
        sums,
        counts,
        args.min_factor,
        args.max_factor,
    )
    write_threshold_table(
        sample_names,
        mean_dps,
        min_dps,
        max_dps,
        args.sample_table,
    )
    return 0


def run_mask(args: argparse.Namespace) -> int:
    sample_names, min_dps, max_dps = read_thresholds_table(args.thresholds_table)
    require_pysam()
    with pysam.VariantFile(args.input_vcf, threads=args.threads) as variant_file:
        header_names = list(variant_file.header.samples)
    if sample_names != header_names:
        raise ValueError("Sample order in thresholds table does not match input VCF/BCF header")
    mask_variants(
        args.input_vcf,
        args.output_bcf,
        min_dps,
        max_dps,
        args.threads,
        args.region,
    )
    return 0


def main() -> int:
    args = parse_args()
    if args.command == "compute-stats":
        return run_compute_stats(args)
    if args.command == "aggregate-stats":
        return run_aggregate_stats(args)
    if args.command == "mask":
        return run_mask(args)
    return run_full(args)


if __name__ == "__main__":
    raise SystemExit(main())

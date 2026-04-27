#!/usr/bin/env python3

from __future__ import annotations

import argparse

try:
    import pysam
except ModuleNotFoundError:
    pysam = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute per-sample missingness from a VCF/BCF file."
    )
    parser.add_argument("--input-vcf", required=True, help="Input VCF/BCF path")
    parser.add_argument("--output-tsv", required=True, help="Output TSV path")
    parser.add_argument("--threads", type=int, default=1, help="HTSlib reader threads")
    return parser.parse_args()


def require_pysam() -> None:
    if pysam is None:
        raise ModuleNotFoundError(
            "pysam is required to run this script. Activate the environment that provides pysam."
        )


def gt_is_missing(gt) -> bool:
    if gt is None:
        return True
    return all(allele is None for allele in gt)


def main() -> int:
    args = parse_args()
    require_pysam()

    with pysam.VariantFile(args.input_vcf, threads=args.threads) as variant_file:
        sample_names = list(variant_file.header.samples)
        total = [0] * len(sample_names)
        missing = [0] * len(sample_names)

        for record in variant_file:
            samples = record.samples
            for idx in range(len(sample_names)):
                gt = samples[idx].get("GT")
                total[idx] += 1
                if gt_is_missing(gt):
                    missing[idx] += 1

    with open(args.output_tsv, "w", encoding="utf-8") as out_handle:
        out_handle.write("sample\tn_genotypes\tn_missing\tmissing_rate\n")
        for sample_name, n_total, n_missing in zip(sample_names, total, missing):
            missing_rate = (n_missing / n_total) if n_total > 0 else 0.0
            out_handle.write(f"{sample_name}\t{n_total}\t{n_missing}\t{missing_rate:.6f}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3

# Step 4 helper for the VCF QC workflow.
# This script performs a two-pass sample-specific DP filter:
# 1. Read the input VCF.GZ and calculate a mean DP for each sample
#    using only non-missing genotypes.
# 2. Re-read the file and mask genotypes whose DP is outside that
#    sample's allowed DP interval.
#
# The output is:
# - a masked VCF.GZ file
# - a per-sample mean DP table

import argparse
import math
from pathlib import Path
from numbers import Real

import pysam


def parse_args() -> argparse.Namespace:
    # Command-line interface used by the shell workflow.
    parser = argparse.ArgumentParser(
        description="Compute sample-specific DP thresholds and mask GT values outside range."
    )
    parser.add_argument("--input-vcf", required=True, help="Input VCF.GZ path")
    parser.add_argument("--output-vcf", required=True, help="Output masked VCF.GZ path")
    parser.add_argument("--sample-table", required=True, help="Output per-sample DP summary TSV")
    parser.add_argument("--min-factor", type=float, required=True, help="Min DP factor")
    parser.add_argument("--max-factor", type=float, required=True, help="Max DP factor")
    return parser.parse_args()


def gt_is_missing(gt) -> bool:
    # pysam stores GT as a tuple such as (0, 1) or (None, None).
    # A genotype is treated as missing if all allele entries are None.
    if gt is None:
        return True
    return not any(allele is not None for allele in gt)


def parse_dp_value(dp_value) -> float | None:
    # FORMAT/DP may be stored as a scalar or, in some cases, as a tuple-like value.
    # This helper normalizes it to a float and returns None when DP is unavailable.
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


def compute_thresholds(
    input_vcf: str,
    min_factor: float,
    max_factor: float,
):
    # First pass over the file:
    # accumulate DP sums and counts for each sample using only callable genotypes.
    with pysam.VariantFile(input_vcf) as variant_file:
        sample_names = list(variant_file.header.samples)
        sample_count = len(sample_names)
        sums = [0.0] * sample_count
        counts = [0] * sample_count

        for record in variant_file:
            for i, sample_name in enumerate(sample_names):
                sample_data = record.samples[sample_name]
                gt = sample_data.get("GT")
                dp = parse_dp_value(sample_data.get("DP"))

                if gt_is_missing(gt) or dp is None:
                    continue

                sums[i] += dp
                counts[i] += 1

    mean_dps = []
    min_dps = []
    max_dps = []

    # Convert the accumulated per-sample means into lower/upper DP thresholds.
    for total, count in zip(sums, counts):
        mean_dp = total / count if count > 0 else 0.0
        mean_dps.append(mean_dp)
        min_dps.append(mean_dp * min_factor)
        max_dps.append(mean_dp * max_factor)

    return sample_names, mean_dps, min_dps, max_dps


def write_sample_tables(
    sample_names,
    mean_dps,
    min_dps,
    max_dps,
    sample_table_path: str,
) -> None:
    # Write the per-sample threshold summary table.
    header = "sample\tmean_dp\tmin_dp\tmax_dp\n"
    with open(sample_table_path, "w", encoding="utf-8") as sample_out:
        sample_out.write(header)
        for sample, mean_dp, min_dp, max_dp in zip(sample_names, mean_dps, min_dps, max_dps):
            row = f"{sample}\t{mean_dp:.6f}\t{min_dp:.6f}\t{max_dp:.6f}\n"
            sample_out.write(row)


def mask_variants(
    input_vcf: str,
    output_vcf: str,
    sample_names,
    min_dps,
    max_dps,
) -> None:
    # Second pass over the file:
    # keep the original record structure but replace GT with missing
    # whenever DP falls outside the sample-specific range.
    with pysam.VariantFile(input_vcf) as in_vcf, pysam.VariantFile(
        output_vcf, "wz", header=in_vcf.header
    ) as out_vcf:
        for record in in_vcf:
            for i, sample_name in enumerate(sample_names):
                sample_data = record.samples[sample_name]
                gt = sample_data.get("GT")
                dp = parse_dp_value(sample_data.get("DP"))

                if gt_is_missing(gt) or dp is None:
                    continue

                if dp < min_dps[i] or dp > max_dps[i]:
                    # Preserve ploidy while masking the genotype call itself.
                    # Other FORMAT fields are kept unchanged.
                    ploidy = len(gt) if gt else 2
                    sample_data["GT"] = tuple(None for _ in range(ploidy))

            out_vcf.write(record)

    # Build a tabix index for downstream use.
    pysam.tabix_index(output_vcf, preset="vcf", force=True)


def main() -> int:
    # Run the two-pass DP masking workflow and write all outputs.
    args = parse_args()

    sample_names, mean_dps, min_dps, max_dps = compute_thresholds(
        args.input_vcf,
        args.min_factor,
        args.max_factor,
    )
    write_sample_tables(
        sample_names,
        mean_dps,
        min_dps,
        max_dps,
        args.sample_table,
    )
    mask_variants(
        args.input_vcf,
        args.output_vcf,
        sample_names,
        min_dps,
        max_dps,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

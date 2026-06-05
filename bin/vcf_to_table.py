#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_info(value):
    if value in ("", "."):
        return {}

    info = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" in item:
            key, raw_value = item.split("=", 1)
            info[key] = raw_value
        else:
            info[item] = "true"
    return info


def format_info(info):
    if not info:
        return ""
    return ";".join(f"{key}={value}" for key, value in info.items())


def parse_sample(format_value, sample_value):
    if format_value in ("", ".") or sample_value in ("", "."):
        return {}

    keys = format_value.split(":")
    values = sample_value.split(":")
    return {key: values[index] if index < len(values) else "" for index, key in enumerate(keys)}


def variant_type(ref, alt):
    if alt in ("", "."):
        return ""
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    if len(ref) < len(alt):
        return "insertion"
    if len(ref) > len(alt):
        return "deletion"
    return "complex"


def allele_specific_value(value, alt_index):
    if value in ("", "."):
        return ""
    values = value.split(",")
    if len(values) == 1:
        return values[0]
    if alt_index < len(values):
        return values[alt_index]
    return ""


def convert_vcf(sample_id, vcf, output):
    fieldnames = [
        "sample_id",
        "vcf_sample",
        "chrom",
        "pos",
        "id",
        "ref",
        "alt",
        "variant_type",
        "qual",
        "filter",
        "depth",
        "allele_depth",
        "genotype",
        "info",
        "format",
        "sample_values",
    ]

    sample_name = ""
    rows = []

    with open_text(vcf) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").split("\t")
                if len(header) > 9:
                    sample_name = header[9]
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 8:
                raise ValueError(f"{vcf}:{line_number}: VCF row has fewer than 8 columns")

            chrom, pos, record_id, ref, alts, qual, filt, info_value = fields[:8]
            format_value = fields[8] if len(fields) > 8 else ""
            sample_value = fields[9] if len(fields) > 9 else ""
            info = parse_info(info_value)
            sample = parse_sample(format_value, sample_value)
            alt_values = alts.split(",") if alts not in ("", ".") else [""]

            for alt_index, alt in enumerate(alt_values):
                rows.append(
                    {
                        "sample_id": sample_id,
                        "vcf_sample": sample_name,
                        "chrom": chrom,
                        "pos": pos,
                        "id": "" if record_id == "." else record_id,
                        "ref": ref,
                        "alt": alt,
                        "variant_type": variant_type(ref, alt),
                        "qual": "" if qual == "." else qual,
                        "filter": "" if filt == "." else filt,
                        "depth": sample.get("DP") or info.get("DP", ""),
                        "allele_depth": allele_specific_value(sample.get("AD", ""), alt_index + 1),
                        "genotype": sample.get("GT", ""),
                        "info": format_info(info),
                        "format": format_value,
                        "sample_values": sample_value,
                    }
                )

    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    convert_vcf(args.sample_id, args.vcf, args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

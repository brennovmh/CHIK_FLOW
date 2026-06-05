#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_single_csv_row(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if len(rows) != 1:
        raise ValueError(f"{path} must contain exactly one data row")
    return rows[0]


def parse_flagstat(path):
    metrics = {
        "reads_total": "",
        "reads_primary": "",
        "reads_mapped": "",
        "reads_mapped_pct": "",
        "reads_properly_paired": "",
        "reads_properly_paired_pct": "",
    }

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            count = line.split(" + ", 1)[0]
            if " in total " in line:
                metrics["reads_total"] = count
            elif " primary" in line and "primary mapped" not in line:
                metrics["reads_primary"] = count
            elif " mapped (" in line and "primary mapped" not in line:
                metrics["reads_mapped"] = count
                metrics["reads_mapped_pct"] = percent_from_flagstat(line)
            elif " properly paired " in line:
                metrics["reads_properly_paired"] = count
                metrics["reads_properly_paired_pct"] = percent_from_flagstat(line)

    return metrics


def percent_from_flagstat(line):
    if "(" not in line or "%" not in line:
        return ""
    return line.split("(", 1)[1].split("%", 1)[0]


def parse_idxstats(path):
    reference_count = 0
    reference_bases = 0
    mapped_reference_reads = 0
    unmapped_reads = 0

    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: idxstats row has fewer than 4 columns")
            name, length, mapped, unmapped = fields[:4]
            try:
                length = int(length)
                mapped = int(mapped)
                unmapped = int(unmapped)
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: idxstats numeric fields are invalid") from error

            if name == "*":
                unmapped_reads += unmapped
                continue
            reference_count += 1
            reference_bases += length
            mapped_reference_reads += mapped

    return {
        "reference_count": reference_count,
        "reference_bases": reference_bases,
        "idxstats_mapped_reads": mapped_reference_reads,
        "idxstats_unmapped_reads": unmapped_reads,
    }


def parse_gene_coverage(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(row)

    if not rows:
        return {
            "features_reported": 0,
            "min_feature_breadth_1x": "",
            "min_feature_breadth_10x": "",
            "lowest_coverage_feature": "",
            "lowest_coverage_feature_type": "",
            "lowest_coverage_feature_breadth_1x": "",
        }

    def as_float(row, key):
        value = row.get(key, "")
        return float(value) if value != "" else 0.0

    lowest = min(rows, key=lambda row: as_float(row, "breadth_1x"))
    return {
        "features_reported": len(rows),
        "min_feature_breadth_1x": f"{min(as_float(row, 'breadth_1x') for row in rows):.6f}",
        "min_feature_breadth_10x": f"{min(as_float(row, 'breadth_10x') for row in rows):.6f}",
        "lowest_coverage_feature": lowest.get("feature_name", ""),
        "lowest_coverage_feature_type": lowest.get("feature_type", ""),
        "lowest_coverage_feature_breadth_1x": f"{as_float(lowest, 'breadth_1x'):.6f}",
    }


def parse_consensus(path):
    bases = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            bases.append(line.strip().upper())

    sequence = "".join(bases)
    length = len(sequence)
    n_count = sequence.count("N")
    acgt_count = sum(sequence.count(base) for base in "ACGT")
    gc_count = sequence.count("G") + sequence.count("C")

    return {
        "consensus_length": length,
        "consensus_n_bases": n_count,
        "consensus_n_fraction": f"{n_count / length:.6f}" if length else "",
        "consensus_gc_fraction": f"{gc_count / acgt_count:.6f}" if acgt_count else "",
    }


def parse_low_coverage_bed(path):
    intervals = 0
    bases = 0
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: BED row has fewer than 3 columns")
            start = int(fields[1])
            end = int(fields[2])
            intervals += 1
            bases += max(0, end - start)
    return {
        "low_coverage_intervals": intervals,
        "low_coverage_bases": bases,
    }


def count_vcf_records(path):
    count = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            if line.strip():
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--flagstat", required=True)
    parser.add_argument("--idxstats", required=True)
    parser.add_argument("--coverage-summary", required=True)
    parser.add_argument("--gene-coverage", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--low-coverage-bed", required=True)
    parser.add_argument("--variants-vcf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    row = {"sample_id": args.sample_id}
    row.update(parse_flagstat(args.flagstat))
    row.update(parse_idxstats(args.idxstats))

    coverage = read_single_csv_row(args.coverage_summary)
    for key in (
        "total_bases",
        "covered_bases_1x",
        "covered_bases_10x",
        "mean_depth",
        "breadth_1x",
        "breadth_10x",
    ):
        row[f"genome_{key}"] = coverage.get(key, "")

    row.update(parse_gene_coverage(args.gene_coverage))
    row.update(parse_consensus(args.consensus))
    row.update(parse_low_coverage_bed(args.low_coverage_bed))
    row["variant_records"] = count_vcf_records(args.variants_vcf)

    fieldnames = [
        "sample_id",
        "reads_total",
        "reads_primary",
        "reads_mapped",
        "reads_mapped_pct",
        "reads_properly_paired",
        "reads_properly_paired_pct",
        "reference_count",
        "reference_bases",
        "idxstats_mapped_reads",
        "idxstats_unmapped_reads",
        "genome_total_bases",
        "genome_covered_bases_1x",
        "genome_covered_bases_10x",
        "genome_mean_depth",
        "genome_breadth_1x",
        "genome_breadth_10x",
        "features_reported",
        "min_feature_breadth_1x",
        "min_feature_breadth_10x",
        "lowest_coverage_feature",
        "lowest_coverage_feature_type",
        "lowest_coverage_feature_breadth_1x",
        "consensus_length",
        "consensus_n_bases",
        "consensus_n_fraction",
        "consensus_gc_fraction",
        "low_coverage_intervals",
        "low_coverage_bases",
        "variant_records",
    ]

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

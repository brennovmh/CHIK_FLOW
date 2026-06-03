#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys


REQUIRED_COLUMNS = ["sample", "fastq_1", "fastq_2"]


def clean_sample_name(sample):
    sample = sample.strip()
    sample = re.sub(r"\s+", "_", sample)
    sample = re.sub(r"[^A-Za-z0-9_.-]", "_", sample)
    return sample


def validate_fastq(path, row_number, column):
    if not path:
        return ""
    if not path.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")):
        raise ValueError(f"Row {row_number}: {column} is not a FASTQ file: {path}")
    if not os.path.exists(path):
        raise ValueError(f"Row {row_number}: {column} does not exist: {path}")
    return os.path.abspath(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != REQUIRED_COLUMNS:
            raise ValueError(
                "Samplesheet must contain exactly these columns: "
                + ",".join(REQUIRED_COLUMNS)
            )

        for index, row in enumerate(reader, start=2):
            sample = clean_sample_name(row["sample"])
            if not sample:
                raise ValueError(f"Row {index}: sample is empty")

            fastq_1 = validate_fastq(row["fastq_1"].strip(), index, "fastq_1")
            fastq_2 = validate_fastq(row["fastq_2"].strip(), index, "fastq_2")
            single_end = "true" if not fastq_2 else "false"

            rows.append(
                {
                    "sample": sample,
                    "fastq_1": fastq_1,
                    "fastq_2": fastq_2,
                    "single_end": single_end,
                }
            )

    if not rows:
        raise ValueError("Samplesheet has no samples")

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", "fastq_1", "fastq_2", "single_end"]
        )
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

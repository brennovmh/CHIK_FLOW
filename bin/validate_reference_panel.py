#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys


VALID_SEQUENCE = re.compile(r"^[ACGTNacgtn.-]+$")


def parse_fasta(path):
    records = []
    current_id = None
    current_description = ""
    current_sequence = []

    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    records.append(
                        (current_id, current_description, "".join(current_sequence))
                    )
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Line {line_number}: FASTA header is empty")
                parts = header.split(maxsplit=1)
                current_id = parts[0]
                current_description = parts[1] if len(parts) > 1 else ""
                current_sequence = []
            else:
                if not current_id:
                    raise ValueError(
                        f"Line {line_number}: sequence found before first FASTA header"
                    )
                if not VALID_SEQUENCE.match(line):
                    raise ValueError(
                        f"Line {line_number}: FASTA sequence contains unsupported bases"
                    )
                current_sequence.append(line.upper())

    if current_id:
        records.append((current_id, current_description, "".join(current_sequence)))

    if not records:
        raise ValueError("Reference FASTA has no records")

    seen = set()
    for record_id, _, sequence in records:
        if record_id in seen:
            raise ValueError(f"Duplicate FASTA record id: {record_id}")
        seen.add(record_id)
        if not sequence:
            raise ValueError(f"FASTA record has no sequence: {record_id}")

    return records


def parse_gff_seqids(path, sequence_lengths):
    seqids = set()
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Line {line_number}: GFF row must have 9 columns")
            seqid = fields[0]
            if seqid not in sequence_lengths:
                seqids.add(seqid)
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as error:
                raise ValueError(
                    f"Line {line_number}: GFF start/end must be integers"
                ) from error
            if start < 1 or end < start or end > sequence_lengths[seqid]:
                raise ValueError(
                    f"Line {line_number}: GFF coordinates exceed FASTA sequence length"
                )
            seqids.add(seqid)
    return seqids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gff")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    records = parse_fasta(args.fasta)
    fasta_ids = {record_id for record_id, _, _ in records}
    sequence_lengths = {record_id: len(sequence) for record_id, _, sequence in records}

    gff_seqids = set()
    if args.gff:
        gff_seqids = parse_gff_seqids(args.gff, sequence_lengths)
        missing = gff_seqids - fasta_ids
        if missing:
            raise ValueError(
                "GFF contains seqids not present in FASTA: " + ",".join(sorted(missing))
            )

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["reference_id", "description", "length", "has_annotation"],
        )
        writer.writeheader()
        for record_id, description, sequence in records:
            writer.writerow(
                {
                    "reference_id": record_id,
                    "description": description,
                    "length": len(sequence),
                    "has_annotation": str(record_id in gff_seqids).lower(),
                }
            )

    if args.gff and not gff_seqids:
        print("WARNING: GFF has no feature rows", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

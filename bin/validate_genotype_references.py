#!/usr/bin/env python3

import argparse
import csv
import re
import sys


VALID_SEQUENCE = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn.-]+$")
VALID_SOURCES = {"wild", "vaccine"}
SOURCE_ALIASES = {
    "wild": "wild",
    "wild-type": "wild",
    "wild_type": "wild",
    "wildtype": "wild",
    "field": "wild",
    "clinical": "wild",
    "vaccine": "vaccine",
    "vaccinal": "vaccine",
    "vaccine-strain": "vaccine",
    "vaccine_strain": "vaccine",
    "live-attenuated": "vaccine",
    "live_attenuated": "vaccine",
}


def normalize_source(value):
    normalized = value.strip().lower().replace(" ", "_")
    return SOURCE_ALIASES.get(normalized)


def parse_metadata(header):
    metadata = {}
    for item in header.replace(";", "|").split("|")[1:]:
        if "=" in item:
            key, value = item.split("=", 1)
            metadata[key.strip().lower()] = value.strip()
    return metadata


def parse_fasta(path):
    records = []
    current_header = None
    current_sequence = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_sequence)))
                current_header = line[1:].strip()
                current_sequence = []
                if not current_header:
                    raise ValueError(f"Line {line_number}: FASTA header is empty")
            else:
                if current_header is None:
                    raise ValueError(
                        f"Line {line_number}: sequence found before first FASTA header"
                    )
                if not VALID_SEQUENCE.match(line):
                    raise ValueError(
                        f"Line {line_number}: FASTA sequence contains unsupported bases"
                    )
                current_sequence.append(line.upper())
    if current_header is not None:
        records.append((current_header, "".join(current_sequence)))
    if not records:
        raise ValueError("Genotype reference FASTA has no records")
    return records


def record_id(header):
    return header.split()[0].split("|")[0]


def validate(records):
    rows = []
    seen_ids = set()
    seen_sources = set()

    for header, sequence in records:
        ref_id = record_id(header)
        if ref_id in seen_ids:
            raise ValueError(f"Duplicate FASTA record id: {ref_id}")
        seen_ids.add(ref_id)
        if not sequence:
            raise ValueError(f"FASTA record has no sequence: {ref_id}")

        metadata = parse_metadata(header)
        missing = [
            key
            for key in ("genotype", "lineage", "source")
            if not metadata.get(key)
        ]
        if missing:
            raise ValueError(
                f"{ref_id}: missing required header metadata: {','.join(missing)}"
            )

        source = normalize_source(metadata["source"])
        if source not in VALID_SOURCES:
            raise ValueError(
                f"{ref_id}: source must be one of wild,vaccine "
                f"(observed: {metadata['source']})"
            )
        seen_sources.add(source)

        rows.append(
            {
                "reference_id": ref_id,
                "genotype": metadata["genotype"],
                "lineage": metadata["lineage"],
                "source": source,
                "length": len(sequence),
            }
        )

    if len(seen_sources) == 1:
        print(
            "WARNING: genotype reference FASTA contains only "
            f"{next(iter(seen_sources))} source records",
            file=sys.stderr,
        )

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = validate(parse_fasta(args.fasta))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["reference_id", "genotype", "lineage", "source", "length"],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

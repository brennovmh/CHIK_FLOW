#!/usr/bin/env python3

import argparse
import csv
import sys


FIELDNAMES = [
    "sample_id",
    "best_reference",
    "genotype",
    "lineage",
    "source",
    "identity",
    "distance",
    "compared_bases",
    "ambiguous_bases",
    "status",
    "note",
]


KNOWN_REFERENCE_LABELS = {
    "NC_004162.2": ("ECSA", "S27-African-prototype", "wild"),
    "AF369024": ("ECSA", "S27-African-prototype", "wild"),
}


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


def parse_fasta(path):
    records = []
    current = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records.append((current, "".join(chunks).upper()))
                current = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        records.append((current, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    return records


def header_id(header):
    return header.split()[0].split("|")[0]


def header_metadata(header):
    metadata = {}
    for item in header.replace(";", "|").split("|")[1:]:
        if "=" in item:
            key, value = item.split("=", 1)
            metadata[key.strip().lower()] = value.strip()
    accession = header_id(header)
    if accession in KNOWN_REFERENCE_LABELS:
        genotype, lineage, source = KNOWN_REFERENCE_LABELS[accession]
        metadata.setdefault("genotype", genotype)
        metadata.setdefault("lineage", lineage)
        metadata.setdefault("source", source)
    metadata["source"] = normalize_source(metadata.get("source", "unknown"))
    return metadata


def normalize_source(value):
    if not value:
        return "unknown"
    normalized = value.strip().lower().replace(" ", "_")
    return SOURCE_ALIASES.get(normalized, "unknown")


def compare_sequences(query, reference):
    compared = 0
    mismatches = 0
    ambiguous = 0
    for query_base, reference_base in zip(query.upper(), reference.upper()):
        if query_base in "N-" or reference_base in "N-":
            ambiguous += 1
            continue
        if query_base not in "ACGT" or reference_base not in "ACGT":
            ambiguous += 1
            continue
        compared += 1
        if query_base != reference_base:
            mismatches += 1
    if compared == 0:
        return None
    distance = mismatches / compared
    return {
        "identity": 1 - distance,
        "distance": distance,
        "compared_bases": compared,
        "ambiguous_bases": ambiguous,
    }


def assign(sample_id, consensus, references):
    best = None
    for ref_header, ref_sequence in references:
        metrics = compare_sequences(consensus, ref_sequence)
        if metrics is None:
            continue
        candidate = {
            "reference_header": ref_header,
            "reference_id": header_id(ref_header),
            **header_metadata(ref_header),
            **metrics,
        }
        if best is None or candidate["distance"] < best["distance"]:
            best = candidate

    if best is None:
        return {
            "sample_id": sample_id,
            "best_reference": "",
            "genotype": "unassigned",
            "lineage": "unassigned",
            "source": "unknown",
            "identity": "",
            "distance": "",
            "compared_bases": 0,
            "ambiguous_bases": len(consensus),
            "status": "failed",
            "note": "No comparable A/C/G/T bases against genotype references",
        }

    genotype = best.get("genotype", "unclassified")
    lineage = best.get("lineage", "unclassified")
    source = best.get("source", "unknown")
    status = "assigned" if genotype != "unclassified" else "nearest_reference_only"
    note = ""
    if len(references) == 1:
        note = "Only one genotype reference was available; assignment is nearest-reference only"
    if source == "unknown":
        note = "; ".join(
            item for item in [note, "Nearest reference has no wild/vaccine source label"] if item
        )

    return {
        "sample_id": sample_id,
        "best_reference": best["reference_id"],
        "genotype": genotype,
        "lineage": lineage,
        "source": source,
        "identity": f"{best['identity']:.6f}",
        "distance": f"{best['distance']:.6f}",
        "compared_bases": best["compared_bases"],
        "ambiguous_bases": best["ambiguous_bases"],
        "status": status,
        "note": note,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--references", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    consensus_records = parse_fasta(args.consensus)
    if len(consensus_records) != 1:
        raise ValueError(f"{args.consensus}: expected exactly one consensus record")

    references = parse_fasta(args.references)
    result = assign(args.sample_id, consensus_records[0][1], references)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerow(result)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

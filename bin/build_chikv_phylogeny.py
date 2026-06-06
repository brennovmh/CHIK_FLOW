#!/usr/bin/env python3

import argparse
import csv
import re
import sys


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


KNOWN_REFERENCE_LABELS = {
    "NC_004162.2": ("ECSA", "S27-African-prototype", "wild"),
    "AF369024": ("ECSA", "S27-African-prototype", "wild"),
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
    return records


def header_id(header):
    return header.split()[0].split("|")[0]


def normalize_source(value):
    if not value:
        return "unknown"
    normalized = value.strip().lower().replace(" ", "_")
    return SOURCE_ALIASES.get(normalized, "unknown")


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
    metadata.setdefault("genotype", "unclassified")
    metadata.setdefault("lineage", "unclassified")
    metadata["source"] = normalize_source(metadata.get("source", "unknown"))
    return metadata


def safe_tree_label(name, source):
    label = f"{name}__source-{source}"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label)


def distance(left, right):
    compared = 0
    mismatches = 0
    for left_base, right_base in zip(left, right):
        if left_base not in "ACGT" or right_base not in "ACGT":
            continue
        compared += 1
        if left_base != right_base:
            mismatches += 1
    if compared == 0:
        return 1.0, compared
    return mismatches / compared, compared


def nearest_reference(sequence, references):
    best = None
    for reference in references:
        current_distance, compared = distance(sequence, reference["sequence"])
        candidate = {
            "nearest_reference": reference["id"],
            "nearest_distance": current_distance,
            "compared_bases": compared,
            "genotype": reference["genotype"],
            "lineage": reference["lineage"],
            "source": reference["source"],
        }
        if best is None or candidate["nearest_distance"] < best["nearest_distance"]:
            best = candidate
    if best is None:
        return {
            "nearest_reference": "",
            "nearest_distance": "",
            "compared_bases": 0,
            "genotype": "unclassified",
            "lineage": "unclassified",
            "source": "unknown",
        }
    best["nearest_distance"] = f"{best['nearest_distance']:.6f}"
    return best


def write_alignment(records, output):
    with open(output, "w") as handle:
        for record in records:
            handle.write(f">{record['tree_label']}\n")
            sequence = record["sequence"]
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")


def write_distance_matrix(records, output):
    names = [record["tree_label"] for record in records]
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", *names])
        for left in records:
            row = [left["tree_label"]]
            for right in records:
                row.append(f"{distance(left['sequence'], right['sequence'])[0]:.6f}")
            writer.writerow(row)


def write_metadata(records, output):
    fieldnames = [
        "tree_label",
        "record_id",
        "role",
        "genotype",
        "lineage",
        "source",
        "nearest_reference",
        "nearest_distance",
        "compared_bases",
    ]
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow({field: record.get(field, "") for field in fieldnames})


def upgma(records):
    if len(records) == 1:
        return f"{records[0]['tree_label']}:0.000000;"

    clusters = {
        record["tree_label"]: {
            "members": [record["tree_label"]],
            "height": 0.0,
            "newick": record["tree_label"],
        }
        for record in records
    }
    sequences = {record["tree_label"]: record["sequence"] for record in records}

    def cluster_distance(left, right):
        values = []
        for left_member in clusters[left]["members"]:
            for right_member in clusters[right]["members"]:
                values.append(distance(sequences[left_member], sequences[right_member])[0])
        return sum(values) / len(values)

    counter = 0
    while len(clusters) > 1:
        names = sorted(clusters)
        best = None
        for index, left in enumerate(names):
            for right in names[index + 1:]:
                current = cluster_distance(left, right)
                if best is None or current < best[0]:
                    best = (current, left, right)
        current_distance, left, right = best
        new_height = current_distance / 2
        left_branch = max(new_height - clusters[left]["height"], 0)
        right_branch = max(new_height - clusters[right]["height"], 0)
        counter += 1
        new_name = f"node{counter}"
        clusters[new_name] = {
            "members": clusters[left]["members"] + clusters[right]["members"],
            "height": new_height,
            "newick": (
                f"({clusters[left]['newick']}:{left_branch:.6f},"
                f"{clusters[right]['newick']}:{right_branch:.6f})"
            ),
        }
        del clusters[left]
        del clusters[right]

    return next(iter(clusters.values()))["newick"] + ";"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus", nargs="+", required=True)
    parser.add_argument("--reference")
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--distances", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--tree", required=True)
    args = parser.parse_args()

    records = []
    references = []
    if args.reference:
        for header, sequence in parse_fasta(args.reference):
            metadata = header_metadata(header)
            record_id = header_id(header)
            reference = {
                "record_id": record_id,
                "id": record_id,
                "role": "reference",
                "genotype": metadata["genotype"],
                "lineage": metadata["lineage"],
                "source": metadata["source"],
                "nearest_reference": record_id,
                "nearest_distance": "0.000000",
                "compared_bases": len(sequence),
                "sequence": sequence,
            }
            references.append(reference)
            records.append(reference)
    for path in args.consensus:
        for header, sequence in parse_fasta(path):
            record_id = header_id(header)
            inferred = nearest_reference(sequence, references)
            records.append(
                {
                    "record_id": record_id,
                    "role": "sample",
                    "sequence": sequence,
                    **inferred,
                }
            )

    deduplicated = []
    seen = set()
    for record in records:
        base_label = safe_tree_label(record["record_id"], record["source"])
        unique_name = base_label
        suffix = 1
        while unique_name in seen:
            suffix += 1
            unique_name = f"{base_label}_{suffix}"
        seen.add(unique_name)
        record["tree_label"] = unique_name
        deduplicated.append(record)

    if not deduplicated:
        raise ValueError("No FASTA records were provided")

    write_alignment(deduplicated, args.alignment)
    write_distance_matrix(deduplicated, args.distances)
    write_metadata(deduplicated, args.metadata)
    with open(args.tree, "w") as handle:
        handle.write(upgma(deduplicated) + "\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

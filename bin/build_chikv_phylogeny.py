#!/usr/bin/env python3

import argparse
import csv
import sys


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
                    records.append((current.split()[0], "".join(chunks).upper()))
                current = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        records.append((current.split()[0], "".join(chunks).upper()))
    return records


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


def write_alignment(records, output):
    with open(output, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")


def write_distance_matrix(records, output):
    names = [name for name, sequence in records]
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", *names])
        for left_name, left_sequence in records:
            row = [left_name]
            for right_name, right_sequence in records:
                row.append(f"{distance(left_sequence, right_sequence)[0]:.6f}")
            writer.writerow(row)


def upgma(records):
    if len(records) == 1:
        return f"{records[0][0]}:0.000000;"

    clusters = {
        name: {
            "members": [name],
            "height": 0.0,
            "newick": name,
        }
        for name, sequence in records
    }
    sequences = {name: sequence for name, sequence in records}

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
    parser.add_argument("--tree", required=True)
    args = parser.parse_args()

    records = []
    if args.reference:
        records.extend(parse_fasta(args.reference))
    for path in args.consensus:
        records.extend(parse_fasta(path))

    deduplicated = []
    seen = set()
    for name, sequence in records:
        unique_name = name
        suffix = 1
        while unique_name in seen:
            suffix += 1
            unique_name = f"{name}_{suffix}"
        seen.add(unique_name)
        deduplicated.append((unique_name, sequence))

    if not deduplicated:
        raise ValueError("No FASTA records were provided")

    write_alignment(deduplicated, args.alignment)
    write_distance_matrix(deduplicated, args.distances)
    with open(args.tree, "w") as handle:
        handle.write(upgma(deduplicated) + "\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

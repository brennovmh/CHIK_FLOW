#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_fasta(path):
    records = []
    current = None
    chunks = []
    with open_text(path) as handle:
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


def parse_fastq_sequences(paths, max_reads):
    sequences = []
    for path in paths:
        with open_text(path) as handle:
            for line_number, line in enumerate(handle):
                if line_number % 4 == 1:
                    sequences.append(line.strip().upper())
                    if max_reads and len(sequences) >= max_reads:
                        return sequences
    if not sequences:
        raise ValueError("No FASTQ reads found")
    return sequences


def kmers(sequence, kmer_size):
    for index in range(0, len(sequence) - kmer_size + 1):
        kmer = sequence[index:index + kmer_size]
        if set(kmer) <= {"A", "C", "G", "T"}:
            yield kmer


def reference_kmers(sequence, kmer_size):
    return set(kmers(sequence, kmer_size))


def score_reference(reads, ref_kmers, kmer_size):
    total = 0
    matches = 0
    matched_reads = 0
    for read in reads:
        read_total = 0
        read_matches = 0
        for kmer in kmers(read, kmer_size):
            read_total += 1
            if kmer in ref_kmers:
                read_matches += 1
        total += read_total
        matches += read_matches
        if read_matches:
            matched_reads += 1
    identity = matches / total if total else 0.0
    return matches, total, matched_reads, identity


def write_fasta(header, sequence, output):
    with open(output, "w") as handle:
        handle.write(f">{header}\n")
        for index in range(0, len(sequence), 80):
            handle.write(sequence[index:index + 80] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--reads", nargs="+", required=True)
    parser.add_argument("--references", required=True)
    parser.add_argument("--selected-fasta", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--kmer-size", type=int, default=31)
    parser.add_argument("--max-reads", type=int, default=10000)
    args = parser.parse_args()

    references = parse_fasta(args.references)
    reads = parse_fastq_sequences(args.reads, args.max_reads)
    rows = []
    best = None

    for header, sequence in references:
        ref_id = header.split()[0].split("|")[0]
        matches, total, matched_reads, identity = score_reference(
            reads,
            reference_kmers(sequence, args.kmer_size),
            args.kmer_size,
        )
        row = {
            "sample_id": args.sample_id,
            "reference_id": ref_id,
            "score": matches,
            "total_kmers": total,
            "matched_reads": matched_reads,
            "read_count": len(reads),
            "kmer_identity": f"{identity:.6f}",
            "selected": "false",
        }
        rows.append(row)
        candidate = (matches, matched_reads, identity, ref_id, header, sequence)
        if best is None or candidate[:4] > best[:4]:
            best = candidate

    if best is None or best[0] == 0:
        raise ValueError("No reference k-mer matches found for sample")

    for row in rows:
        if row["reference_id"] == best[3]:
            row["selected"] = "true"

    write_fasta(best[4], best[5], args.selected_fasta)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "reference_id",
                "score",
                "total_kmers",
                "matched_reads",
                "read_count",
                "kmer_identity",
                "selected",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

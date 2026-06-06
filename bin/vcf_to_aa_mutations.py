#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys
from urllib.parse import unquote


GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


FIELDNAMES = [
    "sample_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "feature_id",
    "feature_name",
    "product",
    "protein_id",
    "cds_position",
    "codon_position",
    "aa_position",
    "ref_codon",
    "alt_codon",
    "ref_aa",
    "alt_aa",
    "mutation",
    "effect",
]


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_attributes(value):
    attributes = {}
    for item in value.split(";"):
        if not item or "=" not in item:
            continue
        key, raw_value = item.split("=", 1)
        attributes[key] = unquote(raw_value)
    return attributes


def feature_name(attributes, fallback):
    for key in ("gene", "Name", "locus_tag", "product", "ID"):
        if attributes.get(key):
            return attributes[key]
    return fallback


def parse_gff_cds(path):
    features = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: GFF row must have 9 columns")

            seqid, source, feature_type, start, end, score, strand, phase, attrs = fields
            if feature_type != "CDS":
                continue

            start = int(start)
            end = int(end)
            if strand not in ("+", "-"):
                continue

            attributes = parse_attributes(attrs)
            fallback = f"{seqid}:{start}-{end}"
            features.append(
                {
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "phase": 0 if phase == "." else int(phase),
                    "id": attributes.get("ID", fallback),
                    "name": feature_name(attributes, fallback),
                    "product": attributes.get("product", ""),
                    "protein_id": attributes.get("protein_id", ""),
                }
            )

    if not features:
        raise ValueError(f"{path}: no CDS features found")
    return features


def parse_fasta(path):
    records = {}
    current = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    records[current] = "".join(chunks).upper()
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current:
        records[current] = "".join(chunks).upper()
    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    return records


def reverse_complement(sequence):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1].upper()


def translate(codon):
    if len(codon) != 3 or any(base not in "ACGT" for base in codon):
        return "X"
    return GENETIC_CODE.get(codon, "X")


def parse_vcf(path):
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_number}: VCF row has fewer than 8 columns")
            chrom, pos, record_id, ref, alts = fields[:5]
            for alt in alts.split(","):
                if alt not in ("", "."):
                    yield {
                        "chrom": chrom,
                        "pos": int(pos),
                        "ref": ref.upper(),
                        "alt": alt.upper(),
                    }


def overlaps(feature, chrom, start, end):
    return feature["seqid"] == chrom and start <= feature["end"] and end >= feature["start"]


def coding_offset(feature, genomic_pos):
    if feature["strand"] == "+":
        return genomic_pos - feature["start"] - feature["phase"]
    return feature["end"] - genomic_pos - feature["phase"]


def codon_bounds(feature, codon_start_offset):
    if feature["strand"] == "+":
        start = feature["start"] + feature["phase"] + codon_start_offset
        return start, start + 2
    end = feature["end"] - feature["phase"] - codon_start_offset
    return end - 2, end


def classify_indel(ref, alt):
    delta = len(alt) - len(ref)
    if delta % 3:
        return "frameshift"
    return "inframe_indel"


def annotate_substitution(sample_id, variant, feature, reference):
    offset = coding_offset(feature, variant["pos"])
    if offset < 0:
        return None

    codon_start_offset = (offset // 3) * 3
    touched_offsets = [
        coding_offset(feature, variant["pos"] + index)
        for index in range(len(variant["ref"]))
        if overlaps(feature, variant["chrom"], variant["pos"] + index, variant["pos"] + index)
    ]
    if any((item // 3) * 3 != codon_start_offset for item in touched_offsets):
        return row(
            sample_id,
            variant,
            feature,
            cds_position=offset + 1,
            codon_position=(offset % 3) + 1,
            aa_position=(codon_start_offset // 3) + 1,
            mutation=f"p.complex_substitution@{(codon_start_offset // 3) + 1}",
            effect="complex_substitution",
        )

    codon_start, codon_end = codon_bounds(feature, codon_start_offset)
    if codon_start < feature["start"] or codon_end > feature["end"]:
        return None

    ref_codon = reference[feature["seqid"]][codon_start - 1:codon_end]
    if feature["strand"] == "-":
        ref_codon = reverse_complement(ref_codon)

    alt_codon = list(ref_codon)
    for index, base in enumerate(variant["alt"]):
        genomic_pos = variant["pos"] + index
        if not overlaps(feature, variant["chrom"], genomic_pos, genomic_pos):
            continue
        alt_offset = coding_offset(feature, genomic_pos)
        codon_index = alt_offset % 3
        alt_codon[codon_index] = reverse_complement(base)[0] if feature["strand"] == "-" else base
    alt_codon = "".join(alt_codon)

    ref_aa = translate(ref_codon)
    alt_aa = translate(alt_codon)
    aa_position = (codon_start_offset // 3) + 1
    effect = "synonymous" if ref_aa == alt_aa else "missense"
    if alt_aa == "*":
        effect = "stop_gained"
    elif ref_aa == "*" and alt_aa != "*":
        effect = "stop_lost"

    return row(
        sample_id,
        variant,
        feature,
        cds_position=offset + 1,
        codon_position=(offset % 3) + 1,
        aa_position=aa_position,
        ref_codon=ref_codon,
        alt_codon=alt_codon,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        mutation=f"{ref_aa}{aa_position}{alt_aa}",
        effect=effect,
    )


def annotate_indel(sample_id, variant, feature):
    offset = coding_offset(feature, variant["pos"])
    if offset < 0:
        return None
    aa_position = (offset // 3) + 1
    effect = classify_indel(variant["ref"], variant["alt"])
    return row(
        sample_id,
        variant,
        feature,
        cds_position=offset + 1,
        codon_position=(offset % 3) + 1,
        aa_position=aa_position,
        mutation=f"p.{effect}@{aa_position}",
        effect=effect,
    )


def row(
    sample_id,
    variant,
    feature,
    cds_position="",
    codon_position="",
    aa_position="",
    ref_codon="",
    alt_codon="",
    ref_aa="",
    alt_aa="",
    mutation="",
    effect="",
):
    return {
        "sample_id": sample_id,
        "chrom": variant["chrom"],
        "pos": variant["pos"],
        "ref": variant["ref"],
        "alt": variant["alt"],
        "feature_id": feature["id"],
        "feature_name": feature["name"],
        "product": feature["product"],
        "protein_id": feature["protein_id"],
        "cds_position": cds_position,
        "codon_position": codon_position,
        "aa_position": aa_position,
        "ref_codon": ref_codon,
        "alt_codon": alt_codon,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "mutation": mutation,
        "effect": effect,
    }


def annotate_variant(sample_id, variant, features, reference):
    variant_end = variant["pos"] + len(variant["ref"]) - 1
    rows = []
    for feature in features:
        if not overlaps(feature, variant["chrom"], variant["pos"], variant_end):
            continue
        if len(variant["ref"]) == len(variant["alt"]):
            annotated = annotate_substitution(sample_id, variant, feature, reference)
        else:
            annotated = annotate_indel(sample_id, variant, feature)
        if annotated:
            rows.append(annotated)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--reference-gff", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    reference = parse_fasta(args.reference_fasta)
    features = parse_gff_cds(args.reference_gff)
    rows = []
    for variant in parse_vcf(args.vcf):
        rows.extend(annotate_variant(args.sample_id, variant, features, reference))

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

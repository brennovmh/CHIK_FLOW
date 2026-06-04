#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict
from urllib.parse import unquote


DEFAULT_FEATURE_TYPES = ("gene", "CDS")


def parse_attributes(value):
    attributes = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" not in item:
            continue
        key, raw_value = item.split("=", 1)
        attributes[key] = unquote(raw_value)
    return attributes


def feature_name(attributes, fallback):
    for key in ("gene", "Name", "locus_tag", "product", "ID"):
        if attributes.get(key):
            return attributes[key]
    return fallback


def parse_gff(path, feature_types):
    features = []
    wanted_types = set(feature_types)

    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Line {line_number}: GFF row must have 9 columns")

            seqid, source, feature_type, start, end, score, strand, phase, attrs = fields
            if feature_type not in wanted_types:
                continue

            try:
                start = int(start)
                end = int(end)
            except ValueError as error:
                raise ValueError(
                    f"Line {line_number}: GFF start/end must be integers"
                ) from error

            if start < 1 or end < start:
                raise ValueError(f"Line {line_number}: invalid GFF interval")

            attributes = parse_attributes(attrs)
            fallback = f"{seqid}:{start}-{end}"
            features.append(
                {
                    "seqid": seqid,
                    "source": source,
                    "feature_type": feature_type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "id": attributes.get("ID", fallback),
                    "name": feature_name(attributes, fallback),
                    "product": attributes.get("product", ""),
                }
            )

    if not features:
        raise ValueError(
            "GFF has no matching feature rows for: " + ",".join(feature_types)
        )

    return features


def parse_depth(path):
    depths = defaultdict(dict)
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"Line {line_number}: depth row must have 3 columns")
            seqid, position, depth = fields[:3]
            try:
                position = int(position)
                depth = int(depth)
            except ValueError as error:
                raise ValueError(
                    f"Line {line_number}: depth position/depth must be integers"
                ) from error
            depths[seqid][position] = depth
    return depths


def summarize_feature(feature, depths):
    seq_depth = depths.get(feature["seqid"], {})
    length = feature["end"] - feature["start"] + 1
    depth_sum = 0
    covered_1x = 0
    covered_10x = 0
    min_depth = None
    max_depth = 0

    for position in range(feature["start"], feature["end"] + 1):
        depth = seq_depth.get(position, 0)
        depth_sum += depth
        if depth >= 1:
            covered_1x += 1
        if depth >= 10:
            covered_10x += 1
        min_depth = depth if min_depth is None else min(min_depth, depth)
        max_depth = max(max_depth, depth)

    mean_depth = depth_sum / length if length else 0
    return {
        "length": length,
        "covered_bases_1x": covered_1x,
        "covered_bases_10x": covered_10x,
        "mean_depth": mean_depth,
        "min_depth": min_depth or 0,
        "max_depth": max_depth,
        "breadth_1x": covered_1x / length if length else 0,
        "breadth_10x": covered_10x / length if length else 0,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--depth", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--feature-types",
        default=",".join(DEFAULT_FEATURE_TYPES),
        help="Comma-separated GFF feature types to summarize.",
    )
    args = parser.parse_args()

    feature_types = [item.strip() for item in args.feature_types.split(",") if item.strip()]
    if not feature_types:
        raise ValueError("--feature-types must include at least one feature type")

    features = parse_gff(args.gff, feature_types)
    depths = parse_depth(args.depth)

    fieldnames = [
        "sample_id",
        "seqid",
        "feature_type",
        "feature_id",
        "feature_name",
        "product",
        "start",
        "end",
        "strand",
        "length",
        "covered_bases_1x",
        "covered_bases_10x",
        "mean_depth",
        "min_depth",
        "max_depth",
        "breadth_1x",
        "breadth_10x",
    ]

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for feature in features:
            summary = summarize_feature(feature, depths)
            writer.writerow(
                {
                    "sample_id": args.sample_id,
                    "seqid": feature["seqid"],
                    "feature_type": feature["feature_type"],
                    "feature_id": feature["id"],
                    "feature_name": feature["name"],
                    "product": feature["product"],
                    "start": feature["start"],
                    "end": feature["end"],
                    "strand": feature["strand"],
                    "length": summary["length"],
                    "covered_bases_1x": summary["covered_bases_1x"],
                    "covered_bases_10x": summary["covered_bases_10x"],
                    "mean_depth": f"{summary['mean_depth']:.4f}",
                    "min_depth": summary["min_depth"],
                    "max_depth": summary["max_depth"],
                    "breadth_1x": f"{summary['breadth_1x']:.6f}",
                    "breadth_10x": f"{summary['breadth_10x']:.6f}",
                }
            )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

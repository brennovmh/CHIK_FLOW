#!/usr/bin/env python3

import argparse
import csv
import html
import sys
from pathlib import Path


SOURCE_COLORS = {
    "wild": "#1b9e77",
    "vaccine": "#d95f02",
    "unknown": "#7570b3",
}


def read_csv(path):
    if not path or not Path(path).exists():
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def html_table(rows):
    if not rows:
        return "<p>No records available.</p>"
    columns = list(rows[0].keys())
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(str(row.get(column, '')))}</td>" for column in columns)
        body.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def to_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_text(path):
    if not path or not Path(path).exists():
        return ""
    return Path(path).read_text().strip()


def metadata_by_label(rows):
    return {row.get("tree_label", ""): row for row in rows}


def source_summary(genotype_rows):
    counts = {"wild": 0, "vaccine": 0, "unknown": 0}
    for row in genotype_rows:
        source = row.get("source") or "unknown"
        if source not in counts:
            source = "unknown"
        counts[source] += 1
    return counts


def alert_rows(sample_rows, genotype_rows):
    alerts = []
    for row in genotype_rows:
        status = row.get("status", "")
        source = row.get("source", "unknown")
        if status in {"low_confidence", "ambiguous", "failed"} or source == "unknown":
            alerts.append(
                {
                    "sample_id": row.get("sample_id", ""),
                    "type": "genotyping",
                    "severity": "high" if status == "failed" else "medium",
                    "message": f"status={status or 'missing'}; source={source}",
                }
            )
    for row in sample_rows:
        breadth = to_float(row.get("genome_breadth_1x", ""))
        n_fraction = to_float(row.get("consensus_n_fraction", ""))
        if breadth < 0.8:
            alerts.append(
                {
                    "sample_id": row.get("sample_id", ""),
                    "type": "coverage",
                    "severity": "medium",
                    "message": f"genome breadth at 1x is {breadth:.3f}",
                }
            )
        if n_fraction > 0.2:
            alerts.append(
                {
                    "sample_id": row.get("sample_id", ""),
                    "type": "consensus",
                    "severity": "medium",
                    "message": f"consensus N fraction is {n_fraction:.3f}",
                }
            )
    return alerts


def summary_cards(sample_rows, genotype_rows, alerts):
    counts = source_summary(genotype_rows)
    cards = [
        ("Samples", str(len(sample_rows))),
        ("Wild-like", str(counts["wild"])),
        ("Vaccine-like", str(counts["vaccine"])),
        ("Unknown", str(counts["unknown"])),
        ("Alerts", str(len(alerts))),
    ]
    return (
        '<div class="cards">'
        + "".join(
            f'<div class="card"><div class="card-label">{html.escape(label)}</div>'
            f'<div class="card-value">{html.escape(value)}</div></div>'
            for label, value in cards
        )
        + "</div>"
    )


def parse_newick_label(text, index):
    start = index
    while index < len(text) and text[index] not in ",():;":
        index += 1
    return text[start:index], index


def parse_newick_length(text, index):
    if index >= len(text) or text[index] != ":":
        return 0.0, index
    index += 1
    start = index
    while index < len(text) and text[index] not in ",();":
        index += 1
    try:
        return float(text[start:index] or 0), index
    except ValueError:
        return 0.0, index


def parse_newick_node(text, index=0):
    if index < len(text) and text[index] == "(":
        index += 1
        children = []
        while index < len(text):
            child, index = parse_newick_node(text, index)
            children.append(child)
            if index < len(text) and text[index] == ",":
                index += 1
                continue
            if index < len(text) and text[index] == ")":
                index += 1
                break
        label, index = parse_newick_label(text, index)
        length, index = parse_newick_length(text, index)
        return {"label": label, "length": length, "children": children}, index

    label, index = parse_newick_label(text, index)
    length, index = parse_newick_length(text, index)
    return {"label": label, "length": length, "children": []}, index


def parse_newick(text):
    if not text:
        return None
    node, _ = parse_newick_node(text.strip().rstrip(";"))
    return node


def leaves(node):
    if not node:
        return []
    if not node["children"]:
        return [node]
    found = []
    for child in node["children"]:
        found.extend(leaves(child))
    return found


def assign_tree_coordinates(node, depth, y_positions, max_depth):
    current_depth = depth + node.get("length", 0.0)
    max_depth[0] = max(max_depth[0], current_depth)
    if not node["children"]:
        node["x_depth"] = current_depth
        node["y_index"] = y_positions[node["label"]]
        return
    for child in node["children"]:
        assign_tree_coordinates(child, current_depth, y_positions, max_depth)
    node["x_depth"] = current_depth
    node["y_index"] = sum(child["y_index"] for child in node["children"]) / len(node["children"])


def draw_tree_edges(node, scale_x, scale_y, x_offset, y_offset, elements):
    x = x_offset + node["x_depth"] * scale_x
    y = y_offset + node["y_index"] * scale_y
    for child in node["children"]:
        child_x = x_offset + child["x_depth"] * scale_x
        child_y = y_offset + child["y_index"] * scale_y
        elements.append(
            f'<line x1="{x:.1f}" y1="{child_y:.1f}" x2="{child_x:.1f}" y2="{child_y:.1f}" stroke="#52606d" stroke-width="1.5"/>'
        )
        elements.append(
            f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x:.1f}" y2="{child_y:.1f}" stroke="#52606d" stroke-width="1.5"/>'
        )
        draw_tree_edges(child, scale_x, scale_y, x_offset, y_offset, elements)


def tree_figure(metadata_rows, tree, css_class="tree-figure"):
    if not metadata_rows:
        return "<p>No phylogeny metadata available.</p>"
    root = parse_newick(tree)
    if not root:
        return "<p>No tree available.</p>"

    leaf_nodes = leaves(root)
    rows_by_label = metadata_by_label(metadata_rows)
    ordered_labels = [leaf["label"] for leaf in leaf_nodes]
    width = 980
    row_height = 34
    top = 54
    left = 28
    tree_width = 400
    label_x = left + tree_width + 24
    height = top + row_height * max(len(ordered_labels), 1) + 56
    y_positions = {label: index for index, label in enumerate(ordered_labels)}
    max_depth = [0.0]
    assign_tree_coordinates(root, 0.0, y_positions, max_depth)
    scale_x = tree_width / max(max_depth[0], 0.000001)
    scale_y = row_height

    elements = [
        f'<svg class="{css_class}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Phylogeny with wild vaccine labels">',
        '<rect x="0" y="0" width="980" height="100%" fill="#ffffff"/>',
        '<text x="32" y="24" class="svg-title">Phylogeny source labels</text>',
    ]

    for source, color in SOURCE_COLORS.items():
        legend_x = 620 + list(SOURCE_COLORS).index(source) * 110
        elements.append(f'<circle cx="{legend_x}" cy="20" r="6" fill="{color}"/>')
        elements.append(f'<text x="{legend_x + 12}" y="24" class="svg-label">{html.escape(source)}</text>')

    draw_tree_edges(root, scale_x, scale_y, left, top, elements)

    for label in ordered_labels:
        row = rows_by_label.get(label, {"tree_label": label, "source": "unknown"})
        y = top + y_positions[label] * row_height
        source = row.get("source") or "unknown"
        color = SOURCE_COLORS.get(source, SOURCE_COLORS["unknown"])
        genotype = row.get("genotype") or "unclassified"
        lineage = row.get("lineage") or "unclassified"
        role = row.get("role") or "record"
        nearest = row.get("nearest_reference") or ""

        elements.extend(
            [
                f'<circle cx="{label_x - 12}" cy="{y}" r="5" fill="{color}"/>',
                f'<text x="{label_x}" y="{y + 4}" class="svg-label">{html.escape(label)}</text>',
                f'<text x="660" y="{y + 4}" class="svg-meta">{html.escape(role)} | {html.escape(source)} | {html.escape(genotype)} | {html.escape(lineage)}</text>',
            ]
        )
        if nearest:
            elements.append(
                f'<text x="880" y="{y + 4}" class="svg-meta">nearest={html.escape(nearest)}</text>'
            )

    if tree:
        compact_tree = tree if len(tree) <= 150 else tree[:147] + "..."
        elements.append(f'<text x="32" y="{height - 18}" class="svg-meta">{html.escape(compact_tree)}</text>')
    elements.append("</svg>")
    return "".join(elements)


def coverage_figure(gene_rows):
    if not gene_rows:
        return "<p>No gene coverage records available.</p>"
    rows = sorted(gene_rows, key=lambda row: (row.get("sample_id", ""), row.get("start", "0"), row.get("feature_name", "")))
    width = 980
    row_height = 24
    left = 260
    top = 36
    bar_width = 620
    height = top + row_height * len(rows) + 30
    elements = [
        f'<svg class="coverage-figure" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Gene coverage by sample">',
        '<rect x="0" y="0" width="980" height="100%" fill="#ffffff"/>',
        '<text x="24" y="22" class="svg-title">Gene coverage breadth at 1x</text>',
    ]
    for index, row in enumerate(rows):
        y = top + index * row_height
        breadth = max(0.0, min(to_float(row.get("breadth_1x", "")), 1.0))
        color = "#1b9e77" if breadth >= 0.9 else "#e6ab02" if breadth >= 0.5 else "#d95f02"
        label = f"{row.get('sample_id', '')} | {row.get('feature_name', '')}"
        elements.extend(
            [
                f'<text x="24" y="{y + 11}" class="svg-label">{html.escape(label[:38])}</text>',
                f'<rect x="{left}" y="{y}" width="{bar_width}" height="14" fill="#edf2f7"/>',
                f'<rect x="{left}" y="{y}" width="{bar_width * breadth:.1f}" height="14" fill="{color}"/>',
                f'<text x="{left + bar_width + 12}" y="{y + 11}" class="svg-meta">{breadth:.3f}</text>',
            ]
        )
    elements.append("</svg>")
    return "".join(elements)


def write_html(output, sample_rows, genotype_rows, tree, phylogeny_rows, gene_rows, phylogeny_svg):
    alerts = alert_rows(sample_rows, genotype_rows)
    phylogeny = tree_figure(phylogeny_rows, tree)
    if phylogeny_svg:
        Path(phylogeny_svg).write_text(tree_figure(phylogeny_rows, tree, css_class="tree-figure-export"))
    css = """
    body { font-family: Arial, sans-serif; margin: 32px; color: #1f2933; }
    h1, h2 { color: #102a43; }
    table { border-collapse: collapse; width: 100%; margin: 16px 0 28px; font-size: 13px; }
    th, td { border: 1px solid #d9e2ec; padding: 7px 9px; text-align: left; }
    th { background: #f0f4f8; }
    code, pre { background: #f0f4f8; padding: 10px; display: block; overflow-x: auto; }
    .meta { color: #52606d; }
    .cards { display: flex; flex-wrap: wrap; gap: 12px; margin: 16px 0 28px; }
    .card { border: 1px solid #d9e2ec; padding: 12px 14px; min-width: 128px; background: #fbfcfd; }
    .card-label { color: #52606d; font-size: 12px; }
    .card-value { color: #102a43; font-size: 24px; font-weight: 700; margin-top: 4px; }
    .tree-figure, .coverage-figure { width: 100%; max-width: 980px; border: 1px solid #d9e2ec; margin: 12px 0 22px; }
    .svg-title { font: 700 16px Arial, sans-serif; fill: #102a43; }
    .svg-label { font: 12px Arial, sans-serif; fill: #102a43; }
    .svg-meta { font: 11px Arial, sans-serif; fill: #52606d; }
    """
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>CHIK-FLOW Report</title>
  <style>{css}</style>
</head>
<body>
  <h1>CHIK-FLOW Report</h1>
  <p class="meta">Automated batch report generated from pipeline CSV outputs.</p>
  <h2>Batch Overview</h2>
  {summary_cards(sample_rows, genotype_rows, alerts)}
  <h2>Alerts</h2>
  {html_table(alerts)}
  <h2>Sample Summary</h2>
  {html_table(sample_rows)}
  <h2>Genotyping</h2>
  {html_table(genotype_rows)}
  <h2>Gene Coverage</h2>
  {coverage_figure(gene_rows)}
  <h2>Phylogeny</h2>
  {phylogeny}
  <h2>Phylogeny Metadata</h2>
  {html_table(phylogeny_rows)}
  <h2>Newick</h2>
  <pre>{html.escape(tree or 'No tree available.')}</pre>
</body>
</html>
"""
    Path(output).write_text(content)


def pdf_escape(value):
    return value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def write_pdf(output, lines):
    text_lines = []
    for line in lines:
        while len(line) > 92:
            text_lines.append(line[:92])
            line = line[92:]
        text_lines.append(line)

    commands = ["BT", "/F1 10 Tf", "50 790 Td"]
    for index, line in enumerate(text_lines[:70]):
        if index:
            commands.append("0 -14 Td")
        commands.append(f"({pdf_escape(line)}) Tj")
    commands.append("ET")
    stream = "\n".join(commands).encode()

    objects = []
    objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    objects.append(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
    objects.append(b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>")
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    objects.append(b"<< /Length " + str(len(stream)).encode() + b" >>\nstream\n" + stream + b"\nendstream")

    content = bytearray(b"%PDF-1.4\n")
    offsets = [0]
    for number, obj in enumerate(objects, start=1):
        offsets.append(len(content))
        content.extend(f"{number} 0 obj\n".encode())
        content.extend(obj)
        content.extend(b"\nendobj\n")
    xref = len(content)
    content.extend(f"xref\n0 {len(objects) + 1}\n".encode())
    content.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        content.extend(f"{offset:010d} 00000 n \n".encode())
    content.extend(f"trailer << /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n".encode())
    Path(output).write_bytes(content)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-summary", required=True)
    parser.add_argument("--genotypes", nargs="*", default=[])
    parser.add_argument("--tree")
    parser.add_argument("--phylogeny-metadata")
    parser.add_argument("--gene-coverages", nargs="*", default=[])
    parser.add_argument("--html", required=True)
    parser.add_argument("--pdf", required=True)
    parser.add_argument("--phylogeny-svg")
    args = parser.parse_args()

    sample_rows = read_csv(args.sample_summary)
    genotype_rows = []
    for path in args.genotypes:
        genotype_rows.extend(read_csv(path))
    tree = read_text(args.tree)
    phylogeny_rows = read_csv(args.phylogeny_metadata)
    gene_rows = []
    for path in args.gene_coverages:
        gene_rows.extend(read_csv(path))

    alerts = alert_rows(sample_rows, genotype_rows)
    write_html(args.html, sample_rows, genotype_rows, tree, phylogeny_rows, gene_rows, args.phylogeny_svg)

    lines = ["CHIK-FLOW Report", "", "Batch Overview"]
    counts = source_summary(genotype_rows)
    lines.append(
        f"samples={len(sample_rows)}, wild={counts['wild']}, vaccine={counts['vaccine']}, unknown={counts['unknown']}, alerts={len(alerts)}"
    )
    lines.extend(["", "Alerts"])
    for row in alerts:
        lines.append(", ".join(f"{key}={value}" for key, value in row.items()))
    lines.extend(["", "Sample Summary"])
    for row in sample_rows:
        lines.append(", ".join(f"{key}={value}" for key, value in row.items()))
    lines.extend(["", "Genotyping"])
    for row in genotype_rows:
        lines.append(", ".join(f"{key}={value}" for key, value in row.items()))
    lines.extend(["", "Phylogeny"])
    for row in phylogeny_rows:
        lines.append(
            "{tree_label}: source={source}, genotype={genotype}, lineage={lineage}".format(
                tree_label=row.get("tree_label", ""),
                source=row.get("source", ""),
                genotype=row.get("genotype", ""),
                lineage=row.get("lineage", ""),
            )
        )
    lines.extend(["", "Newick", tree or "No tree available."])
    write_pdf(args.pdf, lines)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

#!/usr/bin/env python3

import argparse
import base64
import csv
import html
import sys
from datetime import UTC, datetime
from pathlib import Path


DEFAULT_LOGO_URL = "https://github.com/user-attachments/assets/c115b8e6-ce38-4b45-b1c3-e0b82c9da7e5"

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
        return '<p class="empty-state">No records available.</p>'
    columns = list(rows[0].keys())
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(str(row.get(column, '')))}</td>" for column in columns)
        body.append(f"<tr>{cells}</tr>")
    return (
        '<div class="table-wrap">'
        f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"
        "</div>"
    )


def to_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def clamp(value, minimum=0.0, maximum=1.0):
    return max(minimum, min(value, maximum))


def format_percent(value):
    return f"{clamp(value) * 100:.1f}%"


def format_pdf_percent(value):
    return f"{clamp(to_float(value)) * 100:.1f}%"


def coverage_color(value):
    value = clamp(value)
    if value >= 0.95:
        return "#057a55"
    if value >= 0.85:
        return "#2f9e44"
    if value >= 0.70:
        return "#f59f00"
    if value >= 0.50:
        return "#f08c00"
    return "#c92a2a"


def coverage_text_color(value):
    return "#ffffff" if clamp(value) < 0.7 or clamp(value) >= 0.95 else "#102a43"


def feature_sort_key(row):
    feature_type_rank = {
        "gene": 0,
        "CDS": 1,
        "protein_domain": 2,
    }
    priority = {
        "E2": -2,
        "E1": -1,
    }
    name = row.get("feature_name", "")
    return (
        row.get("seqid", ""),
        priority.get(name, 0),
        to_float(row.get("start", "")),
        to_float(row.get("end", "")),
        feature_type_rank.get(row.get("feature_type", ""), 9),
        name,
        row.get("feature_id", ""),
    )


def distinct_features(gene_rows):
    features = {}
    for row in gene_rows:
        name = row.get("feature_name") or row.get("feature_id") or row.get("product") or "feature"
        key = (row.get("seqid", ""), row.get("start", ""), row.get("end", ""), name)
        features.setdefault(key, row)
    return sorted(features.values(), key=feature_sort_key)


def display_feature_name(row):
    name = row.get("feature_name") or row.get("feature_id") or row.get("product") or "feature"
    if name in {"E1", "E2", "E3", "Capsid"}:
        return name
    if name == "CHIKVgp1":
        return "NS polyprotein"
    if name == "CHIKVgp2":
        return "Structural polyprotein"
    if name == "NP_690588.1":
        return "NS CDS"
    if name == "NP_690589.2":
        return "Structural CDS"
    return name


def sample_ids(sample_rows, gene_rows):
    ids = [row.get("sample_id", "") for row in sample_rows if row.get("sample_id", "")]
    for row in gene_rows:
        sample_id = row.get("sample_id", "")
        if sample_id and sample_id not in ids:
            ids.append(sample_id)
    return ids


def read_text(path):
    if not path or not Path(path).exists():
        return ""
    return Path(path).read_text().strip()


def image_data_uri(path, fallback_url=DEFAULT_LOGO_URL):
    if not path:
        return fallback_url
    image_path = Path(path)
    if not image_path.exists():
        return fallback_url
    suffix = image_path.suffix.lower()
    mime = "image/png" if suffix == ".png" else "image/jpeg" if suffix in {".jpg", ".jpeg"} else "image/svg+xml"
    encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def metadata_by_label(rows):
    return {row.get("tree_label", ""): row for row in rows}


def text_width(text, char_width=7, minimum=0, maximum=None):
    width = len(str(text or "")) * char_width
    width = max(minimum, width)
    if maximum is not None:
        width = min(maximum, width)
    return width


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


def report_findings(sample_rows, genotype_rows, alerts):
    findings = []
    if not sample_rows:
        findings.append("No sample summary records were available for interpretation.")
        return findings

    low_genome = [
        row.get("sample_id", "")
        for row in sample_rows
        if to_float(row.get("genome_breadth_1x", "")) < 0.8
    ]
    high_n = [
        row.get("sample_id", "")
        for row in sample_rows
        if to_float(row.get("consensus_n_fraction", "")) > 0.2
    ]
    assigned = [
        row.get("sample_id", "")
        for row in genotype_rows
        if row.get("status", "") not in {"low_confidence", "ambiguous", "failed"}
        and row.get("source", "unknown") != "unknown"
    ]
    uncertain = sorted({row.get("sample_id", "") for row in alerts if row.get("type") in {"coverage", "consensus", "genotyping"}})

    findings.append(f"{len(sample_rows)} sample(s) were included in this batch report.")
    findings.append(f"{len(assigned)} sample(s) have reportable genotype/source assignments.")
    if low_genome:
        findings.append("Low genome breadth was detected in: " + ", ".join(filter(None, low_genome)) + ".")
    else:
        findings.append("All samples met the default genome breadth screening threshold at 1x.")
    if high_n:
        findings.append("Consensus ambiguity above threshold was detected in: " + ", ".join(filter(None, high_n)) + ".")
    if uncertain:
        findings.append("Review recommended for: " + ", ".join(filter(None, uncertain)) + ".")
    elif not alerts:
        findings.append("No automated quality or genotyping alerts were raised.")
    return findings


def findings_html(findings):
    return "<ul class=\"findings\">" + "".join(f"<li>{html.escape(item)}</li>" for item in findings) + "</ul>"


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
    longest_label = max((text_width(label) for label in ordered_labels), default=0)
    width = max(1080, 520 + longest_label + 560)
    row_height = 34
    top = 54
    left = 32
    tree_width = 360
    label_x = left + tree_width + 24
    meta_x = label_x + longest_label + 34
    nearest_x = meta_x + 300
    height = top + row_height * max(len(ordered_labels), 1) + 56
    y_positions = {label: index for index, label in enumerate(ordered_labels)}
    max_depth = [0.0]
    assign_tree_coordinates(root, 0.0, y_positions, max_depth)
    scale_x = tree_width / max(max_depth[0], 0.000001)
    scale_y = row_height

    elements = [
        f'<svg class="{css_class}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Phylogeny with wild vaccine labels">',
        f'<rect x="0" y="0" width="{width}" height="100%" fill="#ffffff"/>',
        '<text x="32" y="24" class="svg-title">Phylogeny source labels</text>',
    ]

    for source, color in SOURCE_COLORS.items():
        legend_x = width - 360 + list(SOURCE_COLORS).index(source) * 110
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
                f'<text x="{meta_x}" y="{y + 4}" class="svg-meta">{html.escape(role)} | {html.escape(source)} | {html.escape(genotype)} | {html.escape(lineage)}</text>',
            ]
        )
        if nearest:
            elements.append(
                f'<text x="{nearest_x}" y="{y + 4}" class="svg-meta">nearest={html.escape(nearest)}</text>'
            )

    if tree:
        compact_tree = tree if len(tree) <= 150 else tree[:147] + "..."
        elements.append(f'<text x="32" y="{height - 18}" class="svg-meta">{html.escape(compact_tree)}</text>')
    elements.append("</svg>")
    return "".join(elements)


def genome_coverage_figure(sample_rows):
    if not sample_rows:
        return "<p>No sample coverage records available.</p>"

    rows = sorted(sample_rows, key=lambda row: row.get("sample_id", ""))
    label_width = max((text_width(row.get("sample_id", "")) for row in rows), default=0)
    left = max(230, min(460, label_width + 44))
    bar_width = 590
    right_space = 255
    width = max(980, left + bar_width + right_space)
    row_height = 32
    top = 52
    height = top + row_height * len(rows) + 42
    elements = [
        f'<svg class="coverage-figure" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Genome coverage breadth by sample">',
        f'<rect x="0" y="0" width="{width}" height="100%" rx="0" fill="#ffffff"/>',
        '<text x="24" y="25" class="svg-title">Genome breadth at 1x</text>',
        f'<text x="{left + bar_width + 14}" y="25" class="svg-meta">1x / 10x / mean</text>',
    ]

    for tick in range(0, 5):
        x = left + bar_width * tick / 4
        elements.append(f'<line x1="{x:.1f}" y1="38" x2="{x:.1f}" y2="{height - 24}" stroke="#e6ebf1" stroke-width="1"/>')

    for index, row in enumerate(rows):
        y = top + index * row_height
        breadth_1x = clamp(to_float(row.get("genome_breadth_1x", "")))
        breadth_10x = clamp(to_float(row.get("genome_breadth_10x", "")))
        mean_depth = to_float(row.get("genome_mean_depth", ""))
        fill = coverage_color(breadth_1x)
        elements.extend(
            [
                f'<text x="24" y="{y + 15}" class="svg-label strong">{html.escape(row.get("sample_id", ""))}</text>',
                f'<rect x="{left}" y="{y}" width="{bar_width}" height="18" rx="4" fill="#edf2f7"/>',
                f'<rect x="{left}" y="{y}" width="{bar_width * breadth_1x:.1f}" height="18" rx="4" fill="{fill}"/>',
                f'<rect x="{left}" y="{y + 21}" width="{bar_width * breadth_10x:.1f}" height="4" rx="2" fill="#2b8a9f"/>',
                f'<text x="{left + bar_width + 14}" y="{y + 14}" class="svg-label">{format_percent(breadth_1x)}</text>',
                f'<text x="{left + bar_width + 76}" y="{y + 14}" class="svg-meta">{format_percent(breadth_10x)} / {mean_depth:.0f}x</text>',
            ]
        )

    elements.append('<text x="24" y="{0}" class="svg-meta">Thin teal line shows 10x breadth.</text>'.format(height - 10))
    elements.append("</svg>")
    return "".join(elements)


def gene_coverage_heatmap(gene_rows):
    if not gene_rows:
        return "<p>No gene coverage records available.</p>"

    features = distinct_features(gene_rows)
    samples = sample_ids([], gene_rows)
    values = {}
    for row in gene_rows:
        name = row.get("feature_name") or row.get("feature_id") or row.get("product") or "feature"
        key = (row.get("seqid", ""), row.get("start", ""), row.get("end", ""), name)
        values[(row.get("sample_id", ""), key)] = clamp(to_float(row.get("breadth_10x", "")))

    cell_width = 82
    row_height = 32
    label_width = max((text_width(sample_id) for sample_id in samples), default=0)
    left = max(190, min(460, label_width + 44))
    top = 104
    width = max(980, left + cell_width * len(features) + 42)
    height = top + row_height * len(samples) + 60
    elements = [
        f'<svg class="coverage-figure coverage-heatmap" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Gene coverage breadth heatmap">',
        f'<rect x="0" y="0" width="{width}" height="100%" fill="#ffffff"/>',
        '<text x="24" y="25" class="svg-title">Gene coverage breadth at >=10x</text>',
        '<text x="24" y="49" class="svg-meta">E2 and E1 are shown first, followed by other annotated genomic features.</text>',
    ]

    legend = [(">=95%", "#057a55"), ("85-95%", "#2f9e44"), ("70-85%", "#f59f00"), ("50-70%", "#f08c00"), ("<50%", "#c92a2a")]
    for index, (label, color) in enumerate(legend):
        x = 24 + index * 78
        elements.extend(
            [
                f'<rect x="{x}" y="64" width="14" height="14" rx="3" fill="{color}"/>',
                f'<text x="{x + 19}" y="76" class="svg-meta">{html.escape(label)}</text>',
            ]
        )

    for col, row in enumerate(features):
        x = left + col * cell_width + cell_width / 2
        label = display_feature_name(row)[:18]
        elements.append(
            f'<text x="{x:.1f}" y="94" class="svg-meta rotated" transform="rotate(-35 {x:.1f} 94)">{html.escape(label)}</text>'
        )

    for row_index, sample_id in enumerate(samples):
        y = top + row_index * row_height
        elements.append(f'<text x="24" y="{y + 20}" class="svg-label strong">{html.escape(sample_id)}</text>')
        for col, feature in enumerate(features):
            name = feature.get("feature_name") or feature.get("feature_id") or feature.get("product") or "feature"
            key = (feature.get("seqid", ""), feature.get("start", ""), feature.get("end", ""), name)
            value = values.get((sample_id, key))
            x = left + col * cell_width
            if value is None:
                fill = "#f1f3f5"
                label = "NA"
                text_color = "#6b7280"
            else:
                fill = coverage_color(value)
                label = f"{value * 100:.0f}"
                text_color = coverage_text_color(value)
            elements.extend(
                [
                    f'<rect x="{x}" y="{y}" width="{cell_width - 4}" height="24" rx="4" fill="{fill}"/>',
                    f'<text x="{x + (cell_width - 4) / 2:.1f}" y="{y + 16}" text-anchor="middle" class="svg-cell" fill="{text_color}">{html.escape(label)}</text>',
                ]
            )

    elements.append("</svg>")
    return "".join(elements)


def coverage_section(sample_rows, gene_rows):
    return genome_coverage_figure(sample_rows) + gene_coverage_heatmap(gene_rows)


def write_html(output, sample_rows, genotype_rows, tree, phylogeny_rows, gene_rows, phylogeny_svg, logo=None):
    alerts = alert_rows(sample_rows, genotype_rows)
    phylogeny = tree_figure(phylogeny_rows, tree)
    findings = report_findings(sample_rows, genotype_rows, alerts)
    generated_at = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    logo_src = image_data_uri(logo)
    if phylogeny_svg:
        Path(phylogeny_svg).write_text(tree_figure(phylogeny_rows, tree, css_class="tree-figure-export"))
    css = """
    :root { color-scheme: light; --ink: #152331; --muted: #5f6f7f; --line: #d8e0e8; --panel: #f7f9fb; --brand: #0f766e; --accent: #ba3a20; --header: #e7f1f2; }
    * { box-sizing: border-box; }
    body { font-family: Arial, Helvetica, sans-serif; margin: 0; color: var(--ink); background: #eef3f7; }
    main { max-width: 1180px; margin: 0 auto; padding: 24px 28px 40px; background: #ffffff; min-height: 100vh; }
    .report-header { display: grid; grid-template-columns: 112px 1fr; gap: 22px; align-items: center; padding: 20px 0 24px; border-bottom: 4px solid var(--brand); }
    .report-logo { width: 108px; height: 108px; object-fit: contain; }
    .eyebrow { color: var(--brand); font-size: 12px; font-weight: 700; letter-spacing: 0; text-transform: uppercase; margin: 0 0 6px; }
    h1 { color: var(--ink); font-size: 34px; line-height: 1.05; margin: 0; }
    h2 { color: var(--ink); font-size: 20px; margin: 32px 0 10px; padding-bottom: 8px; border-bottom: 1px solid var(--line); }
    h3 { color: var(--ink); font-size: 15px; margin: 20px 0 8px; }
    .table-wrap { width: 100%; overflow-x: auto; margin: 14px 0 26px; border: 1px solid var(--line); background: #ffffff; }
    table { border-collapse: collapse; width: max-content; min-width: 100%; margin: 0; font-size: 12px; }
    th, td { border: 1px solid var(--line); padding: 7px 9px; text-align: left; vertical-align: top; max-width: 280px; overflow-wrap: anywhere; }
    th { background: var(--header); color: #243746; font-weight: 700; position: sticky; top: 0; }
    tr:nth-child(even) td { background: #fbfcfd; }
    code, pre { background: #f1f5f9; border: 1px solid var(--line); padding: 12px; display: block; overflow-x: auto; }
    .meta { color: var(--muted); }
    .report-meta { margin: 8px 0 0; color: var(--muted); font-size: 13px; }
    .section-note { color: var(--muted); font-size: 13px; margin: 0 0 10px; }
    .findings { margin: 14px 0 26px; padding-left: 22px; line-height: 1.45; }
    .findings li { margin: 6px 0; }
    .empty-state { margin: 12px 0 24px; padding: 12px 14px; border: 1px solid var(--line); background: var(--panel); color: var(--muted); }
    .cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(145px, 1fr)); gap: 12px; margin: 16px 0 28px; }
    .card { border: 1px solid var(--line); border-left: 5px solid var(--brand); padding: 13px 14px; min-width: 128px; background: var(--panel); }
    .card-label { color: var(--muted); font-size: 12px; font-weight: 700; text-transform: uppercase; }
    .card-value { color: var(--ink); font-size: 26px; font-weight: 700; margin-top: 4px; }
    .figure-wrap { overflow-x: auto; margin: 12px 0 24px; border: 1px solid var(--line); background: #ffffff; }
    .tree-figure, .coverage-figure { width: max-content; min-width: 100%; display: block; }
    .svg-title { font: 700 16px Arial, sans-serif; fill: #102a43; }
    .svg-label { font: 12px Arial, sans-serif; fill: #102a43; }
    .svg-label.strong { font-weight: 700; }
    .svg-meta { font: 11px Arial, sans-serif; fill: #52606d; }
    .svg-cell { font: 700 11px Arial, sans-serif; }
    .rotated { dominant-baseline: middle; }
    @media print {
      body { background: #ffffff; }
      main { max-width: none; padding: 16px; }
      .table-wrap { overflow: visible; border: 0; }
      table { width: 100%; font-size: 9px; }
      th, td { padding: 4px 5px; max-width: 150px; }
      .figure-wrap { overflow: visible; break-inside: avoid; }
      .coverage-figure, .tree-figure { width: 100%; min-width: 0; }
    }
    """
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>ChikScan Report</title>
  <style>{css}</style>
</head>
<body>
  <main>
    <header class="report-header">
      <img class="report-logo" src="{logo_src}" alt="ChikScan logo">
      <div>
        <p class="eyebrow">Chikungunya sequencing surveillance</p>
        <h1>ChikScan Report</h1>
        <p class="report-meta">Automated batch report generated from pipeline CSV outputs. Generated {generated_at}.</p>
      </div>
    </header>

    <h2>Batch Overview</h2>
    {summary_cards(sample_rows, genotype_rows, alerts)}

    <h2>Automated Interpretation</h2>
    {findings_html(findings)}

    <h2>Alerts</h2>
    {html_table(alerts)}

    <h2>Coverage</h2>
    <p class="section-note">Genome and gene-level breadth summarize the fraction of positions covered at the requested thresholds.</p>
    <div class="figure-wrap">{coverage_section(sample_rows, gene_rows)}</div>

    <h2>Sample Summary</h2>
    {html_table(sample_rows)}

    <h2>Genotyping</h2>
    {html_table(genotype_rows)}

    <h2>Phylogeny</h2>
    <div class="figure-wrap">{phylogeny}</div>

    <h2>Phylogeny Metadata</h2>
    {html_table(phylogeny_rows)}

    <h2>Newick</h2>
    <pre>{html.escape(tree or 'No tree available.')}</pre>
  </main>
</body>
</html>
"""
    Path(output).write_text(content)


def pdf_escape(value):
    return value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def pdf_rgb(color):
    color = color.lstrip("#")
    red = int(color[0:2], 16) / 255
    green = int(color[2:4], 16) / 255
    blue = int(color[4:6], 16) / 255
    return f"{red:.3f} {green:.3f} {blue:.3f}"


def pdf_text(x, y, text, size=10, font="F1", color="#152331"):
    safe = pdf_escape(str(text))
    return f"BT /{font} {size} Tf {pdf_rgb(color)} rg {x:.1f} {y:.1f} Td ({safe}) Tj ET"


def pdf_rect(x, y, width, height, color):
    return f"q {pdf_rgb(color)} rg {x:.1f} {y:.1f} {width:.1f} {height:.1f} re f Q"


def pdf_line(x1, y1, x2, y2, color="#d8e0e8", width=1):
    return f"q {pdf_rgb(color)} RG {width:.1f} w {x1:.1f} {y1:.1f} m {x2:.1f} {y2:.1f} l S Q"


def wrap_text(text, max_chars):
    words = str(text).split()
    if not words:
        return [""]
    lines = []
    current = words[0]
    for word in words[1:]:
        if len(current) + len(word) + 1 > max_chars:
            lines.append(current)
            current = word
        else:
            current += " " + word
    lines.append(current)
    return lines


def pdf_short_sample_id(sample_id, limit=34):
    sample_id = str(sample_id or "")
    if len(sample_id) <= limit:
        return sample_id
    return sample_id[: limit - 3] + "..."


def genotype_by_sample(genotype_rows):
    return {row.get("sample_id", ""): row for row in genotype_rows}


def key_gene_rows(gene_rows):
    wanted = {"E2", "E1", "CHIKVgp1", "CHIKVgp2"}
    rows = [row for row in gene_rows if row.get("feature_name") in wanted]
    return sorted(rows, key=lambda row: (row.get("sample_id", ""), feature_sort_key(row)))


def build_pdf_pages(sample_rows, genotype_rows, gene_rows, alerts, findings, generated_at):
    pages = []
    current = []
    y = 744

    def add_page():
        nonlocal current, y
        current = [
            pdf_rect(0, 0, 612, 792, "#ffffff"),
            pdf_rect(0, 760, 612, 32, "#0f766e"),
            pdf_text(40, 771, "ChikScan", 15, "F2", "#ffffff"),
        ]
        y = 730
        pages.append(current)

    def ensure(space):
        nonlocal y
        if y - space < 48:
            add_page()

    def heading(text):
        nonlocal y
        y -= 14
        ensure(38)
        current.append(pdf_text(40, y, text, 15, "F2"))
        current.append(pdf_line(40, y - 8, 572, y - 8))
        y -= 28

    def paragraph(text, size=10, color="#152331", max_chars=86, leading=13):
        nonlocal y
        for line in wrap_text(text, max_chars):
            ensure(leading + 2)
            current.append(pdf_text(48, y, line, size, "F1", color))
            y -= leading

    add_page()
    current.append(pdf_text(40, y, "ChikScan Batch Report", 24, "F2"))
    current.append(pdf_text(40, y - 22, "Chikungunya sequencing surveillance", 11, "F1", "#5f6f7f"))
    current.append(pdf_text(40, y - 38, f"Generated {generated_at}", 9, "F1", "#5f6f7f"))
    y -= 72

    counts = source_summary(genotype_rows)
    cards = [
        ("Samples", len(sample_rows), "#0f766e"),
        ("Wild-like", counts["wild"], SOURCE_COLORS["wild"]),
        ("Vaccine-like", counts["vaccine"], SOURCE_COLORS["vaccine"]),
        ("Unknown", counts["unknown"], SOURCE_COLORS["unknown"]),
        ("Alerts", len(alerts), "#ba3a20" if alerts else "#52606d"),
    ]
    card_width = 96
    for index, (label, value, color) in enumerate(cards):
        x = 40 + index * 106
        current.append(pdf_rect(x, y - 48, card_width, 48, "#f7f9fb"))
        current.append(pdf_rect(x, y - 48, 5, 48, color))
        current.append(pdf_text(x + 12, y - 18, label, 8, "F2", "#5f6f7f"))
        current.append(pdf_text(x + 12, y - 39, value, 18, "F2"))
    y -= 78

    heading("Executive Interpretation")
    for finding in findings[:5]:
        paragraph("- " + finding, max_chars=84)
    if len(findings) > 5:
        paragraph(f"- {len(findings) - 5} additional interpretation note(s) omitted from this compact PDF.", color="#5f6f7f")

    heading("Sample Results")
    genotype_lookup = genotype_by_sample(genotype_rows)
    current.append(pdf_rect(40, y - 17, 532, 19, "#e7f1f2"))
    headers = [("Sample", 44), ("Source", 238), ("Genotype", 306), ("Genome 1x", 382), ("Genome 10x", 454), ("Mean", 526)]
    for label, x in headers:
        current.append(pdf_text(x, y - 11, label, 8, "F2", "#243746"))
    y -= 27
    for row in sorted(sample_rows, key=lambda item: item.get("sample_id", ""))[:18]:
        ensure(22)
        genotype = genotype_lookup.get(row.get("sample_id", ""), {})
        current.append(pdf_line(40, y + 6, 572, y + 6, "#eef2f6", 0.8))
        current.append(pdf_text(44, y - 6, pdf_short_sample_id(row.get("sample_id", "")), 8))
        current.append(pdf_text(238, y - 6, genotype.get("source", "unknown"), 8))
        current.append(pdf_text(306, y - 6, genotype.get("genotype", "unclassified"), 8))
        current.append(pdf_text(382, y - 6, format_pdf_percent(row.get("genome_breadth_1x", "")), 8))
        current.append(pdf_text(454, y - 6, format_pdf_percent(row.get("genome_breadth_10x", "")), 8))
        current.append(pdf_text(526, y - 6, f"{to_float(row.get('genome_mean_depth', '')):.0f}x", 8))
        y -= 20
    if len(sample_rows) > 18:
        paragraph(f"{len(sample_rows) - 18} additional sample(s) omitted from compact PDF. See HTML/CSV outputs for full tables.", color="#5f6f7f")

    heading("E1/E2 Coverage At >=10x")
    key_rows = key_gene_rows(gene_rows)
    if key_rows:
        current.append(pdf_rect(40, y - 17, 532, 19, "#e7f1f2"))
        headers = [("Sample", 44), ("Feature", 250), ("10x breadth", 336), ("Mean depth", 430), ("Min depth", 516)]
        for label, x in headers:
            current.append(pdf_text(x, y - 11, label, 8, "F2", "#243746"))
        y -= 27
        for row in key_rows[:24]:
            ensure(22)
            current.append(pdf_line(40, y + 6, 572, y + 6, "#eef2f6", 0.8))
            current.append(pdf_text(44, y - 6, pdf_short_sample_id(row.get("sample_id", "")), 8))
            current.append(pdf_text(250, y - 6, display_feature_name(row), 8))
            current.append(pdf_text(336, y - 6, format_pdf_percent(row.get("breadth_10x", "")), 8))
            current.append(pdf_text(430, y - 6, f"{to_float(row.get('mean_depth', '')):.0f}x", 8))
            current.append(pdf_text(516, y - 6, f"{to_float(row.get('min_depth', '')):.0f}x", 8))
            y -= 20
    else:
        paragraph("No gene/domain coverage rows were available for PDF summary.", color="#5f6f7f")

    if alerts:
        heading("Alerts")
        for alert in alerts[:10]:
            paragraph(
                "{sample_id}: {type} / {severity} - {message}".format(**alert),
                color="#8a1c0a",
                max_chars=82,
            )
        if len(alerts) > 10:
            paragraph(f"{len(alerts) - 10} additional alert(s) omitted from compact PDF.", color="#5f6f7f")

    heading("Notes")
    paragraph("This compact PDF intentionally omits raw Newick strings, full metadata tables, and exhaustive CSV fields. Use the HTML report and CSV outputs for audit-level detail.", color="#5f6f7f")
    return pages


def write_pdf(output, sample_rows, genotype_rows, gene_rows, alerts, findings, generated_at):
    page_streams = ["\n".join(page).encode() for page in build_pdf_pages(sample_rows, genotype_rows, gene_rows, alerts, findings, generated_at)]

    objects = []
    objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    page_object_numbers = []
    content_object_numbers = []
    next_object = 5
    for _page in page_streams:
        page_object_numbers.append(next_object)
        content_object_numbers.append(next_object + 1)
        next_object += 2
    kids = " ".join(f"{number} 0 R" for number in page_object_numbers).encode()
    objects.append(b"<< /Type /Pages /Kids [" + kids + b"] /Count " + str(len(page_streams)).encode() + b" >>")
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>")
    for page_number, stream in enumerate(page_streams, start=1):
        content_number = content_object_numbers[page_number - 1]
        objects.append(
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 3 0 R /F2 4 0 R >> >> /Contents {content_number} 0 R >>".encode()
        )
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
    parser.add_argument("--logo")
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
    findings = report_findings(sample_rows, genotype_rows, alerts)
    write_html(args.html, sample_rows, genotype_rows, tree, phylogeny_rows, gene_rows, args.phylogeny_svg, args.logo)
    generated_at = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    write_pdf(args.pdf, sample_rows, genotype_rows, gene_rows, alerts, findings, generated_at)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

#!/usr/bin/env python3

import argparse
import csv
import html
import sys
from pathlib import Path


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


def read_text(path):
    if not path or not Path(path).exists():
        return ""
    return Path(path).read_text().strip()


def write_html(output, sample_rows, genotype_rows, tree):
    css = """
    body { font-family: Arial, sans-serif; margin: 32px; color: #1f2933; }
    h1, h2 { color: #102a43; }
    table { border-collapse: collapse; width: 100%; margin: 16px 0 28px; font-size: 13px; }
    th, td { border: 1px solid #d9e2ec; padding: 7px 9px; text-align: left; }
    th { background: #f0f4f8; }
    code, pre { background: #f0f4f8; padding: 10px; display: block; overflow-x: auto; }
    .meta { color: #52606d; }
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
  <h2>Sample Summary</h2>
  {html_table(sample_rows)}
  <h2>Genotyping</h2>
  {html_table(genotype_rows)}
  <h2>Phylogeny</h2>
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
    parser.add_argument("--html", required=True)
    parser.add_argument("--pdf", required=True)
    args = parser.parse_args()

    sample_rows = read_csv(args.sample_summary)
    genotype_rows = []
    for path in args.genotypes:
        genotype_rows.extend(read_csv(path))
    tree = read_text(args.tree)

    write_html(args.html, sample_rows, genotype_rows, tree)

    lines = ["CHIK-FLOW Report", "", "Sample Summary"]
    for row in sample_rows:
        lines.append(", ".join(f"{key}={value}" for key, value in row.items()))
    lines.extend(["", "Genotyping"])
    for row in genotype_rows:
        lines.append(", ".join(f"{key}={value}" for key, value in row.items()))
    lines.extend(["", "Phylogeny", tree or "No tree available."])
    write_pdf(args.pdf, lines)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

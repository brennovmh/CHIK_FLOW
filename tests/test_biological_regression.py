#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def run(command, expect_success=True):
    result = subprocess.run(
        command,
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if expect_success and result.returncode != 0:
        raise AssertionError(
            f"Command failed: {' '.join(command)}\nSTDOUT={result.stdout}\nSTDERR={result.stderr}"
        )
    if not expect_success and result.returncode == 0:
        raise AssertionError(f"Command unexpectedly succeeded: {' '.join(command)}")
    return result


def read_one_csv(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise AssertionError(f"{path}: expected one row, observed {len(rows)}")
    return rows[0]


def write(path, content):
    path.write_text(content)
    return path


def test_wild_and_vaccine_assignment(tmpdir):
    refs = write(
        tmpdir / "refs.fasta",
        ">wild_ref|genotype=ECSA|lineage=wild_test|source=wild\n"
        "ACGTACGTACGT\n"
        ">vaccine_ref|genotype=Asian|lineage=vaccine_test|source=vaccine\n"
        "TTTTACGTACGT\n",
    )
    wild_consensus = write(tmpdir / "wild.fasta", ">wild_sample\nACGTACGTACGT\n")
    vaccine_consensus = write(tmpdir / "vaccine.fasta", ">vaccine_sample\nTTTTACGTACGT\n")

    run([
        "python3",
        "bin/validate_genotype_references.py",
        "--fasta",
        str(refs),
        "--output",
        str(tmpdir / "validated.csv"),
    ])

    run([
        "python3",
        "bin/assign_chikv_genotype.py",
        "--sample-id",
        "wild_sample",
        "--consensus",
        str(wild_consensus),
        "--references",
        str(refs),
        "--output",
        str(tmpdir / "wild.csv"),
    ])
    wild = read_one_csv(tmpdir / "wild.csv")
    assert wild["source"] == "wild"
    assert wild["best_reference"] == "wild_ref"

    run([
        "python3",
        "bin/assign_chikv_genotype.py",
        "--sample-id",
        "vaccine_sample",
        "--consensus",
        str(vaccine_consensus),
        "--references",
        str(refs),
        "--output",
        str(tmpdir / "vaccine.csv"),
    ])
    vaccine = read_one_csv(tmpdir / "vaccine.csv")
    assert vaccine["source"] == "vaccine"
    assert vaccine["best_reference"] == "vaccine_ref"


def test_invalid_headers_fail(tmpdir):
    invalid_source = write(
        tmpdir / "invalid_source.fasta",
        ">bad|genotype=ECSA|lineage=x|source=laboratory\nACGT\n",
    )
    missing_source = write(
        tmpdir / "missing_source.fasta",
        ">bad|genotype=ECSA|lineage=x\nACGT\n",
    )

    run(
        [
            "python3",
            "bin/validate_genotype_references.py",
            "--fasta",
            str(invalid_source),
            "--output",
            str(tmpdir / "invalid.csv"),
        ],
        expect_success=False,
    )
    run(
        [
            "python3",
            "bin/validate_genotype_references.py",
            "--fasta",
            str(missing_source),
            "--output",
            str(tmpdir / "missing.csv"),
        ],
        expect_success=False,
    )


def test_no_comparable_bases(tmpdir):
    refs = write(
        tmpdir / "refs.fasta",
        ">wild_ref|genotype=ECSA|lineage=wild_test|source=wild\nACGTACGT\n",
    )
    consensus = write(tmpdir / "ambiguous.fasta", ">ambiguous\nNNNNNNNN\n")
    run([
        "python3",
        "bin/assign_chikv_genotype.py",
        "--sample-id",
        "ambiguous",
        "--consensus",
        str(consensus),
        "--references",
        str(refs),
        "--output",
        str(tmpdir / "ambiguous.csv"),
    ])
    row = read_one_csv(tmpdir / "ambiguous.csv")
    assert row["status"] == "failed"
    assert row["source"] == "unknown"
    assert row["compared_bases"] == "0"


def main():
    with tempfile.TemporaryDirectory() as directory:
        tmpdir = Path(directory)
        test_wild_and_vaccine_assignment(tmpdir)
        test_invalid_headers_fail(tmpdir)
        test_no_comparable_bases(tmpdir)
    print("Biological regression tests passed")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

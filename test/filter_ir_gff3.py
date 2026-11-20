#!/usr/bin/env python3
import argparse
from collections import Counter
from pathlib import Path

CHROMOSOMES = {
    "CM061116.1": "galeata",
    "CM125067.1": "obtusa",
    "CM029730.1": "pulicaria",
    "NC_060020.1": "pulex",
    "CM058171.1": "magna",
    "CM023190.1": "sinensis",
    "NC_081333.1": "carinata",
}


def main():
    parser = argparse.ArgumentParser(
        description="Extract IR annotations on the seven test chromosomes."
    )
    parser.add_argument("input", type=Path, help="Combined input GFF3")
    parser.add_argument("output", type=Path, help="Filtered output GFF3")
    args = parser.parse_args()

    if args.input.resolve() == args.output.resolve():
        parser.error("Input and output must be different files.")

    counts = Counter()
    records = []

    with args.input.open() as source:
        for line_number, line in enumerate(source, 1):
            if line.startswith("##FASTA"):
                break

            if line.startswith("##gff-version"):
                continue

            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) != 4:
                    parser.error(
                        f"Malformed sequence-region directive at line {line_number}."
                    )
                if fields[1] in CHROMOSOMES:
                    records.append(line)
                continue

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                continue

            if fields[0] in CHROMOSOMES:
                records.append(line)
                counts[fields[0]] += 1

    missing = [chrom for chrom in CHROMOSOMES if counts[chrom] == 0]
    if missing:
        parser.error(
            "No feature records found for: "
            + ", ".join(missing)
            + ". Check the chromosome identifiers in GFF3 column 1."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as target:
        target.write("##gff-version 3\n")
        for record in records:
            target.write(record.rstrip("\r\n") + "\n")

    print("chromosome\tspecies\tfeature_records")
    for chrom, species in CHROMOSOMES.items():
        print(f"{chrom}\t{species}\t{counts[chrom]}")


if __name__ == "__main__":
    main()


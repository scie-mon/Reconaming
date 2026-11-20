#!/usr/bin/env python3
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Replace the GFF3 source column with IRnotator."
    )
    parser.add_argument("input", type=Path, help="Input GFF3")
    parser.add_argument("output", type=Path, help="Output GFF3")
    args = parser.parse_args()

    if args.input.resolve() == args.output.resolve():
        parser.error("Input and output must be different files.")

    updated = 0
    skipped = 0
    in_fasta = False

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.input.open() as source, args.output.open("w") as target:
        for line in source:
            if in_fasta:
                target.write(line)
                continue

            if line.startswith("##FASTA"):
                in_fasta = True
                target.write(line)
                continue

            if line.startswith("#") or not line.strip():
                target.write(line)
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                skipped += 1
                continue

            fields[1] = "IRnotator"
            target.write("\t".join(fields) + "\n")
            updated += 1

    print(f"Updated feature records: {updated}")
    print(f"Skipped malformed lines: {skipped}")


if __name__ == "__main__":
    main()


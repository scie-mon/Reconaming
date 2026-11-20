#!/usr/bin/env python3
import argparse
from pathlib import Path


def update_cds(line, line_number):
    fields = line.rstrip("\r\n").split("\t")

    if len(fields) != 9 or fields[2] != "CDS":
        return line, False

    attributes = fields[8].split(";")
    name = None

    for attribute in attributes:
        key, separator, value = attribute.partition("=")
        if separator and key == "Name":
            name = value
            break

    if not name:
        raise ValueError(
            f"Line {line_number}: CDS feature has no nonempty Name attribute."
        )

    # Remove the exact suffix, not any combination of its characters.
    transcript = name.removesuffix(" CDS")

    if not transcript:
        raise ValueError(
            f"Line {line_number}: Name produces an empty transcript value."
        )

    updated = []
    transcript_written = False

    for attribute in attributes:
        if not attribute:
            continue

        key, separator, _ = attribute.partition("=")

        if separator and key == "transcript":
            if not transcript_written:
                updated.append(f"transcript={transcript}")
                transcript_written = True
        else:
            updated.append(attribute)

    if not transcript_written:
        updated.append(f"transcript={transcript}")

    fields[8] = ";".join(updated)
    return "\t".join(fields) + "\n", True


def main():
    parser = argparse.ArgumentParser(
        description="Set each CDS transcript attribute from its Name attribute."
    )
    parser.add_argument("input", type=Path, help="Input GFF3")
    parser.add_argument("output", type=Path, help="Updated GFF3")
    args = parser.parse_args()

    if args.input.resolve() == args.output.resolve():
        parser.error("Input and output must be different files.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".part")
    updated_count = 0
    in_fasta = False

    try:
        with args.input.open() as source, temporary.open("w") as target:
            for line_number, line in enumerate(source, 1):
                if line.startswith("##FASTA"):
                    in_fasta = True

                if in_fasta or line.startswith("#") or not line.strip():
                    target.write(line)
                    continue

                result, updated = update_cds(line, line_number)
                target.write(result)
                updated_count += int(updated)

        temporary.replace(args.output)
    finally:
        temporary.unlink(missing_ok=True)

    print(f"Updated CDS records: {updated_count}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}")


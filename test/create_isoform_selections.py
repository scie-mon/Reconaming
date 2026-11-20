#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


INITIAL = (
    "NC_060020.1",  # pulex
    "CM058171.1",   # magna
    "CM023190.1",   # sinensis
    "NC_081333.1",  # carinata
)

ADDITIONAL = (
    "CM061116.1",   # galeata
    "CM125067.1",   # obtusa
    "CM029730.1",   # pulicaria
)

CHROMOSOMES = set(INITIAL + ADDITIONAL)


def attributes(text):
    result = {}
    for item in text.split(";"):
        key, separator, value = item.partition("=")
        if separator:
            result[unquote(key)] = unquote(value)
    return result


def read_lookup(path):
    identifiers = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    raise ValueError("Reference FASTA contains an empty header.")
                identifiers.add(header.split()[0])

    if not identifiers:
        raise ValueError("No FASTA headers found in the reference file.")

    return identifiers


def read_gff(path):
    transcripts = {}
    aliases = defaultdict(set)

    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[0] not in CHROMOSOMES:
                continue

            chrom, _, feature = fields[:3]
            attr = attributes(fields[8])

            if feature.lower() in {"mrna", "transcript"}:
                transcript_id = attr.get("ID", "")
                # Split before decoding, preserving encoded commas in IDs.
                raw_attr = dict(
                    item.split("=", 1)
                    for item in fields[8].split(";")
                    if "=" in item
                )
                parents = [
                    unquote(value)
                    for value in raw_attr.get("Parent", "").split(",")
                    if value
                ]

                if not transcript_id or len(parents) != 1:
                    raise ValueError(
                        f"Line {number}: transcript requires an ID "
                        "and exactly one gene Parent."
                    )

                entry = (chrom, parents[0])
                if transcript_id in transcripts:
                    raise ValueError(
                        f"Duplicate transcript ID: {transcript_id}"
                    )
                transcripts[transcript_id] = entry

            elif feature == "CDS":
                short_id = attr.get("transcript", "")
                if not short_id:
                    short_id = attr.get("Name", "").removesuffix(" CDS")
                if not short_id:
                    raise ValueError(
                        f"Line {number}: CDS has no transcript or Name value."
                    )

                raw_attr = dict(
                    item.split("=", 1)
                    for item in fields[8].split(";")
                    if "=" in item
                )
                for parent in raw_attr.get("Parent", "").split(","):
                    if parent:
                        aliases[unquote(parent)].add(short_id)

    by_gene = defaultdict(list)
    for transcript_id, short_ids in aliases.items():
        if transcript_id not in transcripts:
            raise ValueError(
                f"CDS Parent has no corresponding transcript: {transcript_id}"
            )
        if len(short_ids) != 1:
            raise ValueError(
                f"Inconsistent CDS transcript values for {transcript_id}: "
                f"{sorted(short_ids)}"
            )

        chrom, gene_id = transcripts[transcript_id]
        by_gene[(chrom, gene_id)].append(
            (transcript_id, next(iter(short_ids)))
        )

    return by_gene


def write_table(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "transcript_id"])
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Build four- and seven-species isoform selection tables."
    )
    parser.add_argument(
        "--gff", type=Path,
        default=Path("test/seven_species_IRs_transcript.gff3"),
    )
    parser.add_argument(
        "--reference", type=Path,
        default=Path(
            "/home/primeline/Data/maker/metagenome/"
            "14_phylogeny/IR_iGluR.faa"
        ),
    )
    parser.add_argument("--outdir", type=Path, default=Path("test"))
    args = parser.parse_args()

    lookup = read_lookup(args.reference)
    by_gene = read_gff(args.gff)
    rows_by_chrom = defaultdict(list)
    unresolved = []
    ambiguous = []
    seen_genes = set()

    for (chrom, gene_id), candidates in sorted(by_gene.items()):
        if gene_id in seen_genes:
            raise ValueError(f"Gene ID is not globally unique: {gene_id}")
        seen_genes.add(gene_id)

        if len(candidates) == 1:
            continue

        matches = [
            transcript_id
            for transcript_id, short_id in candidates
            if short_id in lookup
        ]

        if len(matches) > 1:
            ambiguous.append(f"{gene_id}: {', '.join(matches)}")
            continue

        selected = matches[0] if matches else ""
        rows_by_chrom[chrom].append((gene_id, selected))

        if not matches:
            unresolved.append(gene_id)

    if ambiguous:
        raise ValueError(
            "Multiple lookup matches for the same gene; "
            "no tables were written:\n" + "\n".join(ambiguous)
        )

    rows4 = [
        row for chrom in INITIAL for row in rows_by_chrom[chrom]
    ]
    rows7 = rows4 + [
        row for chrom in ADDITIONAL for row in rows_by_chrom[chrom]
    ]

    args.outdir.mkdir(parents=True, exist_ok=True)

    lookup_path = args.outdir / "selected_protein_ids.txt"
    lookup_path.write_text(
        "".join(identifier + "\n" for identifier in sorted(lookup))
    )

    path4 = args.outdir / "isoform_selection_4_species.tsv"
    path7 = args.outdir / "isoform_selection_7_species.tsv"
    write_table(path4, rows4)
    write_table(path7, rows7)

    print(f"Lookup: {lookup_path} ({len(lookup)} identifiers)")
    for path, rows in ((path4, rows4), (path7, rows7)):
        selected = sum(bool(transcript) for _, transcript in rows)
        print(
            f"{path}: {len(rows)} multi-isoform genes; "
            f"{selected} selected; {len(rows) - selected} unresolved"
        )

    if unresolved:
        print("\nNo lookup match; selection left empty:")
        for gene_id in unresolved:
            print(f"  {gene_id}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}")


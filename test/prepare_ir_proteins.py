#!/usr/bin/env python3
import argparse
import base64
import csv
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Seq import Seq


SPECIES = {
    "NC_060020.1": "pulex",
    "CM058171.1": "magna",
    "CM023190.1": "sinensis",
    "NC_081333.1": "carinata",
    "CM061116.1": "galeata",
    "CM125067.1": "obtusa",
    "CM029730.1": "pulicaria",
}
INITIAL = set(list(SPECIES)[:4])


def parse_attributes(text):
    return {
        unquote(key): unquote(value)
        for item in text.split(";")
        if "=" in item
        for key, value in [item.split("=", 1)]
    }


def parents(text):
    raw = dict(
        item.split("=", 1)
        for item in text.split(";")
        if "=" in item
    )
    return [
        unquote(value)
        for value in raw.get("Parent", "").split(",")
        if value
    ]


def canonical_id(species, transcript_id):
    payload = f"{species}\x1f{transcript_id}".encode()
    return "seq_" + base64.urlsafe_b64encode(payload).decode().rstrip("=")


def read_gff(path):
    transcripts = {}
    cds = defaultdict(list)

    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[0] not in SPECIES:
                continue

            chrom, _, feature, start, end, _, strand, phase, text = fields
            attr = parse_attributes(text)
            parent_ids = parents(text)

            if feature.lower() in {"mrna", "transcript"}:
                tid = attr.get("ID", "")
                if not tid or len(parent_ids) != 1:
                    raise ValueError(
                        f"Line {number}: transcript needs ID and one gene Parent."
                    )
                if tid in transcripts:
                    raise ValueError(f"Duplicate transcript ID: {tid}")
                transcripts[tid] = {
                    "chrom": chrom,
                    "gene_id": parent_ids[0],
                    "attributes": attr,
                }

            elif feature == "CDS":
                if not parent_ids:
                    raise ValueError(f"Line {number}: CDS has no Parent.")
                for tid in parent_ids:
                    cds[tid].append({
                        "chrom": chrom,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                        "attributes": attr,
                    })

    return transcripts, cds


def clean_protein(text):
    return "".join(text.split()).upper().removesuffix("*")


def protein_for(transcript, segments, genome):
    corrected = set()
    for attr in [transcript["attributes"]] + [
        segment["attributes"] for segment in segments
    ]:
        value = attr.get("translation", "")
        if value.strip():
            corrected.add(clean_protein(value))

    if len(corrected) > 1:
        raise ValueError("Conflicting translation attributes.")

    if corrected:
        protein = next(iter(corrected))
        method = "translation_attribute"
    else:
        strands = {segment["strand"] for segment in segments}
        contigs = {segment["chrom"] for segment in segments}
        if len(strands) != 1 or contigs != {transcript["chrom"]}:
            raise ValueError("Inconsistent CDS chromosome or strand.")

        strand = next(iter(strands))
        if strand not in {"+", "-"}:
            raise ValueError("Invalid CDS strand.")

        ordered = sorted(
            segments, key=lambda row: row["start"], reverse=(strand == "-")
        )
        fragments = []
        for segment in ordered:
            start, end = segment["start"], segment["end"]
            if not 1 <= start <= end <= len(genome):
                raise ValueError("CDS interval outside chromosome.")
            fragment = genome[start - 1:end]
            if strand == "-":
                fragment = fragment.reverse_complement()
            fragments.append(str(fragment))

        coding = "".join(fragments)
        if len(coding) % 3:
            raise ValueError("CDS length is not divisible by three.")

        protein = str(Seq(coding).translate(table=1)).removesuffix("*")
        method = "genomic_CDS"

    if not protein:
        raise ValueError("Empty protein sequence.")
    if "*" in protein:
        raise ValueError("Internal stop codon in protein.")
    if any(letter not in "ACDEFGHIKLMNPQRSTVWYBXZJUO" for letter in protein):
        raise ValueError("Unexpected character in protein.")

    return protein, method


def read_selections(path):
    selections = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["gene_id", "transcript_id"]:
            raise ValueError(f"{path}: expected gene_id and transcript_id columns.")

        for row in reader:
            gene = row["gene_id"].strip()
            tid = row["transcript_id"].strip()
            if not gene:
                continue
            if gene in selections:
                raise ValueError(f"{path}: duplicate gene {gene}")
            if not tid:
                raise ValueError(f"{path}: unresolved selection for {gene}")
            selections[gene] = tid
    return selections


def write_tsv(path, header, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_fasta(path, records):
    with path.open("w") as handle:
        for identifier, protein in records:
            handle.write(f">{identifier}\n")
            for offset in range(0, len(protein), 80):
                handle.write(protein[offset:offset + 80] + "\n")


def representatives(chromosomes, selections, transcripts, cds, proteins):
    groups = defaultdict(list)
    for tid in cds:
        transcript = transcripts[tid]
        if transcript["chrom"] in chromosomes:
            groups[transcript["gene_id"]].append(tid)

    unknown = set(selections) - set(groups)
    if unknown:
        raise ValueError(f"Selections contain unknown genes: {sorted(unknown)}")

    chosen = []
    for gene, tids in sorted(groups.items()):
        if len(tids) == 1:
            tid = tids[0]
        else:
            tid = selections.get(gene)
            if tid not in tids:
                raise ValueError(f"No valid representative selected for {gene}")

        if tid not in proteins:
            raise ValueError(
                f"Representative {tid} for {gene} could not be translated. "
                "See protein_extraction_report.tsv."
            )
        chosen.append(tid)
    return chosen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gff", type=Path,
        default=Path("test/seven_species_IRs_transcript.gff3"),
    )
    parser.add_argument(
        "--chromosomes", type=Path,
        default=Path(".test-data/chromosomes"),
    )
    parser.add_argument(
        "--selection4", type=Path,
        default=Path("test/isoform_selection_4_species.tsv"),
    )
    parser.add_argument(
        "--selection7", type=Path,
        default=Path("test/isoform_selection_7_species.tsv"),
    )
    parser.add_argument(
        "--outdir", type=Path, default=Path("test/proteins"),
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    transcripts, cds = read_gff(args.gff)
    missing = set(cds) - set(transcripts)
    if missing:
        raise ValueError(f"Missing transcript features: {sorted(missing)}")

    selection4 = read_selections(args.selection4)
    selection7 = read_selections(args.selection7)
    for gene, tid in selection4.items():
        if selection7.get(gene) != tid:
            raise ValueError(f"Four-/seven-species selection differs for {gene}")

    proteins = {}
    metadata = {}
    report = []
    short_ids = set()

    for chrom, species in SPECIES.items():
        path = args.chromosomes / f"{chrom}.fna"
        record = SeqIO.read(path, "fasta")
        if record.id != chrom:
            raise ValueError(f"{path}: expected FASTA identifier {chrom}")

        for tid in sorted(cds):
            transcript = transcripts[tid]
            if transcript["chrom"] != chrom:
                continue

            aliases = {
                segment["attributes"].get("transcript", "")
                or segment["attributes"].get("Name", "").removesuffix(" CDS")
                for segment in cds[tid]
            }
            if len(aliases) != 1 or not next(iter(aliases)):
                raise ValueError(f"Inconsistent/missing short transcript ID: {tid}")

            short_id = next(iter(aliases))
            if short_id in short_ids:
                raise ValueError(f"Duplicate short transcript ID: {short_id}")
            short_ids.add(short_id)

            sequence_id = short_id
            metadata[tid] = (
                sequence_id, short_id, transcript["gene_id"], species
            )

            try:
                protein, method = protein_for(
                    transcript, cds[tid], record.seq
                )
            except ValueError as error:
                report.append([
                    tid, short_id, species, "excluded", "", str(error)
                ])
                continue

            proteins[tid] = protein
            report.append([
                tid, short_id, species, "accepted", method, ""
            ])

    write_tsv(
        args.outdir / "protein_extraction_report.tsv",
        ["transcript_id", "short_id", "species_id",
         "status", "translation_source", "reason"],
        report,
    )

    write_fasta(
        args.outdir / "all_isoforms.faa",
        [(metadata[tid][1], proteins[tid]) for tid in sorted(proteins)],
    )
    write_tsv(
        args.outdir / "protein_identifiers.tsv",
        ["sequence_id", "short_id", "gene_id", "species_id", "transcript_id"],
        [list(metadata[tid]) + [tid] for tid in sorted(proteins)],
    )

    for size, chromosomes, selection in [
        (4, INITIAL, selection4),
        (7, set(SPECIES), selection7),
    ]:
        chosen = representatives(
            chromosomes, selection, transcripts, cds, proteins
        )
        write_fasta(
            args.outdir / f"representatives_{size}_species.faa",
            [(metadata[tid][0], proteins[tid]) for tid in chosen],
        )
        write_tsv(
            args.outdir / f"id_species_{size}_species.tsv",
            ["sequence_id", "species_id"],
            [(metadata[tid][0], metadata[tid][3]) for tid in chosen],
        )
        print(f"{size} species: {len(chosen)} representatives")

    accepted = len(proteins)
    corrected = sum(row[4] == "translation_attribute" for row in report)
    print(f"Candidate isoforms translated: {accepted}")
    print(f"Using corrected translation attributes: {corrected}")
    print(f"Excluded candidate isoforms: {len(cds) - accepted}")
    print(f"Output directory: {args.outdir}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}")


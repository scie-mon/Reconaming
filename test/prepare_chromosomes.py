#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import os
import shutil
import urllib.request
from pathlib import Path

FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

# species, chromosome, assembly directory, initial four-species set, role
DATASETS = [
    (
        "pulex", "NC_060020.1",
        "GCF/021/134/715/GCF_021134715.1_ASM2113471v1",
        True, "focal",
    ),
    (
        "magna", "CM058171.1",
        "GCA/030/254/905/GCA_030254905.1_UOB_LRV0_1",
        True, "focal",
    ),
    (
        "sinensis", "CM023190.1",
        "GCA/013/167/095/GCA_013167095.2_Dsi",
        True, "focal",
    ),
    (
        "carinata", "NC_081333.1",
        "GCF/022/539/665/GCF_022539665.2_CSIRO_AGI_Dcar_HiC_V3",
        True, "focal",
    ),
    (
        "galeata", "CM061116.1",
        "GCA/030/770/115/GCA_030770115.1_ASM3077011v1",
        False, "focal",
    ),
    (
        "obtusa", "CM125067.1",
        "GCA/052/054/635/GCA_052054635.1_ASM5205463v1",
        False, "focal",
    ),
    (
        "pulicaria", "CM029730.1",
        "GCA/017/493/165/GCA_017493165.1_DPuli_LK16",
        False, "focal",
    ),
    (
        "Simocephalus", "CM082768.1",
        "GCA/041/003/635/GCA_041003635.1_ASM4100363v1",
        False, "outgroup",
    ),
]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, destination):
    if destination.exists():
        print(f"  Reusing assembly: {destination}")
        return

    temporary = destination.with_name(destination.name + ".part")
    request = urllib.request.Request(
        url, headers={"User-Agent": "Reconaming-test-data-preparation/1.0"}
    )

    print(f"  Downloading: {url}")
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            with temporary.open("wb") as target:
                shutil.copyfileobj(response, target, length=1024 * 1024)

        if temporary.stat().st_size == 0:
            raise ValueError("Downloaded file is empty.")

        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def validate_chromosome(path, accession):
    headers = []
    length = 0

    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line[1:].split()[0])
            elif line.strip():
                if not headers:
                    raise ValueError(f"{path}: sequence before FASTA header.")
                length += len("".join(line.split()))

    if headers != [accession]:
        raise ValueError(
            f"{path}: expected one record named {accession}; found {headers}."
        )
    if length == 0:
        raise ValueError(f"{path}: empty chromosome sequence.")

    return length


def extract_chromosome(assembly, destination, accession):
    temporary = destination.with_name(destination.name + ".part")
    matches = 0
    selected = False

    print(f"  Extracting: {accession}")
    try:
        with gzip.open(assembly, "rt") as source:
            with temporary.open("w") as target:
                for line in source:
                    if line.startswith(">"):
                        identifier = line[1:].split()[0]
                        selected = identifier == accession
                        if selected:
                            matches += 1

                    if selected:
                        target.write(line)

        # Scan the entire assembly so duplicate records and gzip errors
        # are detected before publishing the extracted chromosome.
        if matches != 1:
            raise ValueError(
                f"{assembly}: expected exactly one {accession} record; "
                f"found {matches}."
            )

        length = validate_chromosome(temporary, accession)
        os.replace(temporary, destination)
        return length
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(
        description="Download assemblies and extract the eight test chromosomes."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(".test-data"),
        help="Destination directory (default: .test-data).",
    )
    args = parser.parse_args()

    downloads = args.data_dir / "downloads"
    chromosomes = args.data_dir / "chromosomes"
    downloads.mkdir(parents=True, exist_ok=True)
    chromosomes.mkdir(parents=True, exist_ok=True)

    rows = []

    for species, accession, assembly_directory, initial, role in DATASETS:
        assembly_name = assembly_directory.rsplit("/", 1)[1]
        filename = f"{assembly_name}_genomic.fna.gz"
        url = f"{FTP_ROOT}/{assembly_directory}/{filename}"
        assembly_path = downloads / filename
        chromosome_path = chromosomes / f"{accession}.fna"

        print(f"\n{species}: {accession}")

        if chromosome_path.exists():
            length = validate_chromosome(chromosome_path, accession)
            print(f"  Reusing validated chromosome: {chromosome_path}")
        else:
            download(url, assembly_path)
            length = extract_chromosome(
                assembly_path, chromosome_path, accession
            )

        checksum = sha256(chromosome_path)
        print(f"  Length: {length:,} bp")
        print(f"  SHA-256: {checksum}")

        rows.append({
            "species_id": species,
            "role": role,
            "initial_focal_set": "yes" if initial else "no",
            "chromosome_accession": accession,
            "assembly_accession": assembly_name.split("_", 2)[0]
                                  + "_"
                                  + assembly_name.split("_", 2)[1],
            "source_url": url,
            "chromosome_file": str(chromosome_path.resolve()),
            "sequence_length": length,
            "chromosome_sha256": checksum,
        })

    manifest = args.data_dir / "chromosome_manifest.tsv"
    temporary = manifest.with_name(manifest.name + ".part")
    try:
        with temporary.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, manifest)
    finally:
        temporary.unlink(missing_ok=True)

    print(f"\nPrepared {len(rows)} chromosomes.")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, EOFError) as error:
        raise SystemExit(f"ERROR: {error}")

#!/bin/bash -euo pipefail
python3 prepare_proteins.py --genome CM058173.1.fa --gff merged.rev0.gff3 --species-id 'Species_1' --gene-id-attribute 'ID' --transcript-id-attribute 'ID' --protein-id-attribute '' --proteins isoform_proteins.faa --manifest isoform_manifest.tsv --report protein_preparation_report.tsv

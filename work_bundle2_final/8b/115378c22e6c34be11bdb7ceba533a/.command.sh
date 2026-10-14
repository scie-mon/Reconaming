#!/bin/bash -euo pipefail
shopt -s nullglob
proteins=(proteins/*/*.faa)
manifests=(manifests/*/*.tsv)
(( ${#proteins[@]} > 0 )) || { echo 'No isoform protein FASTA files received.' >&2; exit 1; }
(( ${#manifests[@]} > 0 )) || { echo 'No isoform manifest TSV files received.' >&2; exit 1; }

cat "${proteins[@]}" > isoform_proteins.faa
head -n 1 "${manifests[0]}" > isoform_manifest.tsv
for manifest in "${manifests[@]}"; do
    tail -n +2 "$manifest"
done >> isoform_manifest.tsv

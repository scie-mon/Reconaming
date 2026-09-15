process MERGE_ISOFORM_OUTPUTS {
    tag 'merge isoform outputs'

    input:
    path protein_files, stageAs: 'proteins/??/*'
    path manifest_files, stageAs: 'manifests/??/*'

    output:
    path 'isoform_proteins.faa', emit: proteins
    path 'isoform_manifest.tsv', emit: manifest

    script:
    """
    shopt -s nullglob
    proteins=(proteins/*/*.faa)
    manifests=(manifests/*/*.tsv)
    (( \${#proteins[@]} > 0 )) || { echo 'No isoform protein FASTA files received.' >&2; exit 1; }
    (( \${#manifests[@]} > 0 )) || { echo 'No isoform manifest TSV files received.' >&2; exit 1; }

    cat "\${proteins[@]}" > isoform_proteins.faa
    head -n 1 "\${manifests[0]}" > isoform_manifest.tsv
    for manifest in "\${manifests[@]}"; do
        tail -n +2 "\$manifest"
    done >> isoform_manifest.tsv
    """
}

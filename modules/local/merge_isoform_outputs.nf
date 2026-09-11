process MERGE_ISOFORM_OUTPUTS {
    tag 'merge focal-genome isoforms'
    input:
    path fastas
    path manifests
    output:
    path 'all_isoform_proteins.faa', emit: proteins
    path 'all_isoform_manifest.tsv', emit: manifest
    script:
    """
    python3 ${projectDir}/bin/merge_isoform_outputs.py \\
      --fastas ${fastas} --manifests ${manifests} \\
      --proteins all_isoform_proteins.faa --manifest all_isoform_manifest.tsv
    """
}

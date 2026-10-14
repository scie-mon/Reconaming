#!/bin/bash -euo pipefail
select_species_tree_inputs.py --representative-manifest representative_manifest.tsv --sequence-inputs test_species_tree_inputs.tsv --busco-mode auto --output selected_species_inputs.tsv

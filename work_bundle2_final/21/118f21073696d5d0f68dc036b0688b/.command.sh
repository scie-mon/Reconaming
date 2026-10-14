#!/bin/bash -euo pipefail
assess_busco_completeness.py --busco-dirs Species_2 Species_3 Root_species Species_1 --warning-threshold 90 --output busco_completeness.tsv

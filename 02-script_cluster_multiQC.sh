#!/bin/bash

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


# chargement des modules nécessaires
module load multiqc/1.13 


# répertoire de base (le répertoire depuis lequel vous lancez le script)
base_dir="$PWD"

echo "=============================================================="
echo "Running MultiQC to generate a summary report"
echo "=============================================================="
mkdir -p "${base_dir}/multiQC_reads"
multiqc -o "${base_dir}/multiQC_reads" "${base_dir}/reads_qc" 
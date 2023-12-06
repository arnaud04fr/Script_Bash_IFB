#!/bin/bash

#SBATCH --mem=96G
#SBATCH --cpus-per-task=8

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


# chargement des modules nécessaires
module load star/2.7.9a


# répertoire de base et les données (le répertoire depuis lequel vous lancez le script)
base_dir="$PWD"
# répertoire contenant les fichiers du génome de référence
# (séquence et annotations)
genome_dir="${base_dir}/genome"
# chemin et nom du fichier contenant le génome de référence.
genome_file="${genome_dir}/GRCh38.primary_assembly.genome.fa"
# chemin et nom du fichier contenant les annotations
annotation_file="${genome_dir}/gencode.v44.basic.annotation.gtf"


# On indexe le génome qu'une seule fois.
echo "=============================================================="
echo "Indexer le génome de référence"
echo "=============================================================="
mkdir -p "${base_dir}/genome_index"
srun STAR --runThreadN "${SLURM_CPUS_PER_TASK}" \
--runMode genomeGenerate \
--genomeDir "${base_dir}/genome_index" \
--genomeFastaFiles "${genome_file}" \
--sjdbGTFfile "${annotation_file}" \
--sjdbOverhang 100 \
--genomeSAindexNbases 14
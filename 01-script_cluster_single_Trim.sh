#!/bin/bash

#SBATCH --partition=long
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user= your_email@mail.com
#SBATCH --array=0-7

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


# chargement des modules nécessaires
module load fastqc/0.11.9
module load trimmomatic/0.39 
module load star/2.7.10b
module load samtools/1.15.1
module load htseq/0.13.5
module load cufflinks/2.2.1

# répertoire de base (le répertoire depuis lequel vous lancez le script)
base_dir="$PWD"
# répertoire contenant les données
data_dir="/path_to_folder"
# répertoire contenant les fichiers du génome de référence
# (séquence et annotations)
genome_dir="${data_dir}/genome"
# chemin et nom du fichier contenant le génome de référence.
genome_file="${genome_dir}/GRCh38.primary_assembly.genome.fa"
# chemin et nom du fichier contenant les annotations
annotation_file="${genome_dir}/gencode.v44.basic.annotation.gtf"
# répertoire contenant les fichiers .fastq.gz
fastq_dir="${base_dir}/reads" #à adapter
fastq_files=("${fastq_dir}"/*.fastq.gz)
# Output directory for filtered FastQ files
filtered_fastq_dir="${base_dir}/filtered_reads"

# extraction de l'identifiant de l'échantillon
# à partir du nom de fichier :
sample=$(basename -s .fastq.gz "${fastq_files[$SLURM_ARRAY_TASK_ID]}")


echo "=============================================================="
echo "Contrôler la qualité : échantillon ${sample}"
echo "=============================================================="
mkdir -p "${base_dir}/reads_qc"
srun fastqc "${fastq_dir}/${sample}.fastq.gz" --outdir "${base_dir}/reads_qc"


echo "=============================================================="
echo "Filtrage des séquences d'adaptatrices:  échantillons ${sample}"
echo "=============================================================="
mkdir -p "${base_dir}/filtered_reads"
# Input and output files
input_fastq="${fastq_dir}/${sample}.fastq.gz"
output_fastq="${filtered_fastq_dir}/${sample}_filtered.fastq.gz"

# Run trimmomatic to filter out all adapter sequences
trimmomatic SE -phred33  "${input_fastq}" "${output_fastq}" ILLUMINACLIP:"${base_dir}/TruSeq3-SE.fa":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25

echo "=============================================================="
echo "Contrôler la qualité après filtration: échantillon ${sample}"
echo "=============================================================="

mkdir -p "${base_dir}/filtered_reads_qc"
srun fastqc "${base_dir}/filtered_reads/${sample}_filtered.fastq.gz" --outdir "${base_dir}/filtered_reads_qc"


echo "=============================================================="
echo "Aligner les reads sur le génome de référence : échantillon ${sample}"
echo "=============================================================="

mkdir -p "${base_dir}/reads_map"
srun STAR --runThreadN "${SLURM_CPUS_PER_TASK}" \
--runMode alignReads \
--genomeDir "${data_dir}/genome_index" \
--sjdbGTFfile "${annotation_file}" \
--readFilesCommand zcat \
--readFilesIn "${filtered_fastq_dir}/${sample}_filtered.fastq.gz" \
--outFilterType BySJout \
--alignIntronMin 10 \
--alignIntronMax 3000 \
--outFileNamePrefix "${base_dir}/reads_map/${sample}_" \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMtype BAM Unsorted


echo "=============================================================="
echo "Trier les reads alignés : échantillon ${sample}"
echo "=============================================================="
srun samtools sort "${base_dir}/reads_map/${sample}_Aligned.out.bam" \
-o "${base_dir}/reads_map/${sample}_Aligned.sorted.out.bam"


echo "=============================================================="
echo "Indexer les reads alignés : échantillon ${sample}"
echo "=============================================================="
srun samtools index "${base_dir}/reads_map/${sample}_Aligned.sorted.out.bam"


echo "=============================================================="
echo "Compter les reads : échantillon ${sample}"
echo "=============================================================="
mkdir -p "${base_dir}/counts/${sample}"
srun htseq-count --order=pos --stranded=reverse \
--mode=intersection-nonempty \
"${base_dir}/reads_map/${sample}_Aligned.sorted.out.bam" \
"${annotation_file}" > "${base_dir}/counts/${sample}/count_${sample}.txt"


echo "=============================================================="
echo "Compter les transcrits : échantillon ${sample}"
echo "=============================================================="
srun cuffquant --num-threads "${SLURM_CPUS_PER_TASK}" \
--library-type=fr-firststrand "${annotation_file}" \
"${base_dir}/reads_map/${sample}_Aligned.sorted.out.bam" \
--output-dir "${base_dir}/counts/${sample}"

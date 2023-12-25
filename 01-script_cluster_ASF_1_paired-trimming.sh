#!/bin/bash

#SBATCH --partition=long
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@mail.com
#SBATCH --array=0-5

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
data_dir="/shared/projects/202304_duo/amailleux/DU_project_WATIPF"
# répertoire contenant les fichiers du génome de référence
# (séquence et annotations)
genome_dir="${data_dir}/genome"
# chemin et nom du fichier contenant le génome de référence.
genome_file="${genome_dir}/GRCh38.primary_assembly.genome.fa"
# chemin et nom du fichier contenant les annotations
annotation_file="${genome_dir}/gencode.v44.basic.annotation.gtf"
fastq_dir="${base_dir}/reads"
# Output directory for filtered FastQ files
filtered_fastq_dir="${base_dir}/filtered_reads"

# liste de tous les fichiers _R1 dans un tableau
fastq_files=(${fastq_dir}/*_1.fastq.gz) #à adpater en fonction des fichiers.
# extraction de l'identifiant de l'échantillon
sample=$(basename -s _1.fastq.gz "${fastq_files[$SLURM_ARRAY_TASK_ID]}")
# Nom du fichier _R2 correspondant
fastq_file_r2="${fastq_dir}/${sample}_2.fastq.gz" #à adpater en fonction des fichiers avant filtrage.
filtered_fastq_file_r2="${fastq_dir}/${sample}_2_filtered.fastq.gz" #à adpater en fonction des fichiers après filtrage.


echo "=============================================================="
echo "Contrôler la qualité : échantillon ${sample}"
echo "=============================================================="
mkdir -p "${base_dir}/reads_qc"
srun fastqc "${fastq_dir}/${sample}_1.fastq.gz" --outdir "${base_dir}/reads_qc"
srun fastqc "${fastq_file_r2}" --outdir "${base_dir}/reads_qc"

echo "=============================================================="
echo "Filtrage des séquences adaptatrices et de mauvaise qualité:  échantillons ${sample}"
echo "=============================================================="
mkdir -p "${base_dir}/filtered_reads"
# Input and output files for R1
input_fastq="${fastq_dir}/${sample}_1.fastq.gz"
output_fastq="${filtered_fastq_dir}/${sample}_1_filtered.fastq.gz"
echo "vérification valeur variable -R1"
echo $input_fastq
echo $output_fastq
# Run trimmomatic to filter out all adapter sequences (need to get TrueSeq file in base_dir
trimmomatic SE -phred33  "${input_fastq}" "${output_fastq}" ILLUMINACLIP:"${base_dir}/TruSeq3-SE.fa":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25

# Input and output files for R2
input_fastq="${fastq_dir}/${sample}_2.fastq.gz"
output_fastq="${filtered_fastq_dir}/${sample}_2_filtered.fastq.gz"
echo "vérification valeur variable -R2"
echo $input_fastq
echo $output_fastq
# Run trimmomatic to filter out all adapter sequences
trimmomatic SE -phred33  "${input_fastq}" "${output_fastq}" ILLUMINACLIP:"${base_dir}/TruSeq3-SE.fa":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25

echo "=============================================================="
echo "Contrôler la qualité après filtration: échantillon ${sample}"
echo "=============================================================="

mkdir -p "${base_dir}/filtered_reads_qc"
srun fastqc "${base_dir}/filtered_reads/${sample}_1_filtered.fastq.gz" --outdir "${base_dir}/filtered_reads_qc"
srun fastqc "${base_dir}/filtered_reads/${sample}_2_filtered.fastq.gz" --outdir "${base_dir}/filtered_reads_qc"


echo "=============================================================="
echo "Aligner les reads sur le génome de référence : échantillon ${sample}"
echo "=============================================================="
mkdir -p "${base_dir}/reads_map"
srun STAR --runThreadN "${SLURM_CPUS_PER_TASK}" \
--runMode alignReads \
--genomeDir "${data_dir}/genome_index" \
--sjdbGTFfile "${annotation_file}" \
--readFilesCommand zcat \
--readFilesIn "${filtered_fastq_dir}/${sample}_1_filtered.fastq.gz" "${filtered_fastq_file_r2}" \
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

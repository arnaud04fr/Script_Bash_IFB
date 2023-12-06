#!/bin/bash

#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=youremail@mail.com

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


# Specify the input FASTQ.gz files
input_files=("XX-1.fastq.gz" "XX-2.fastq.gz" "XX-3.fastq.gz")

# Specify the output file (compressed FASTQ.gz)
output_file="cXX.fastq.gz"

echo "########################################"
echo "###    Concatenation $input_files   ###"   
echo "########################################"

# Concatenate and compress the input files
for file in "${input_files[@]}"; do
  zcat "$file" >> "temp.fastq"
done

# Compress the concatenated file and remove the temporary file
gzip -c "temp.fastq" > "$output_file"
rm "temp.fastq"

echo "Concatenation and compression complete. Output saved to $output_file."

###############################################################################################

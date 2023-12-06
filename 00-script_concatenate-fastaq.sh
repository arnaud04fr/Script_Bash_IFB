#!/bin/bash

#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arnaud04@gmail.com

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


# Specify the input FASTQ.gz files
input_files=("SRR18918146.fastq.gz" "SRR18918147.fastq.gz" "SRR18918148.fastq.gz" "SRR18918149.fastq.gz")

# Specify the output file (compressed FASTQ.gz)
output_file="cSRR18918146.fastq.gz"

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

# Specify the input FASTQ.gz files
input_files=("SRR18918150.fastq.gz" "SRR18918151.fastq.gz" "SRR18918152.fastq.gz" "SRR18918153.fastq.gz")

echo "########################################"
echo "###    Concatenation $input_files   ###"   
echo "########################################"

# Specify the output file (compressed FASTQ.gz)
output_file="cSRR18918150.fastq.gz"

# Concatenate and compress the input files
for file in "${input_files[@]}"; do
  zcat "$file" >> "temp.fastq"
done

# Compress the concatenated file and remove the temporary file
gzip -c "temp.fastq" > "$output_file"
rm "temp.fastq"

echo "Concatenation and compression complete. Output saved to $output_file."

###############################################################################################

# Specify the input FASTQ.gz files
input_files=("SRR18918154.fastq.gz" "SRR18918155.fastq.gz" "SRR18918156.fastq.gz" "SRR18918157.fastq.gz")

echo "########################################"
echo "###    Concatenation $input_files   ###"   
echo "########################################"

# Specify the output file (compressed FASTQ.gz)
output_file="cSRR18918154.fastq.gz"

# Concatenate and compress the input files
for file in "${input_files[@]}"; do
  zcat "$file" >> "temp.fastq"
done

# Compress the concatenated file and remove the temporary file
gzip -c "temp.fastq" > "$output_file"
rm "temp.fastq"

echo "Concatenation and compression complete. Output saved to $output_file."

###############################################################################################

# Specify the input FASTQ.gz files
input_files=("SRR18918158.fastq.gz" "SRR18918159.fastq.gz" "SRR18918160.fastq.gz" "SRR18918161.fastq.gz")


echo "########################################"
echo "###    Concatenation $input_files   ###"   
echo "########################################"

# Specify the output file (compressed FASTQ.gz)
output_file="cSRR18918158.fastq.gz"

# Concatenate and compress the input files
for file in "${input_files[@]}"; do
  zcat "$file" >> "temp.fastq"
done

# Compress the concatenated file and remove the temporary file
gzip -c "temp.fastq" > "$output_file"
rm "temp.fastq"

echo "Concatenation and compression complete. Output saved to $output_file."

###############################################################################################

# Specify the input FASTQ.gz files
input_files=("SRR18918162.fastq.gz" "SRR18918163.fastq.gz" "SRR18918164.fastq.gz" "SRR18918165.fastq.gz")


echo "########################################"
echo "###    Concatenation $input_files   ###"   
echo "########################################"

# Specify the output file (compressed FASTQ.gz)
output_file="cSRR18918162.fastq.gz"

# Concatenate and compress the input files
for file in "${input_files[@]}"; do
  zcat "$file" >> "temp.fastq"
done

# Compress the concatenated file and remove the temporary file
gzip -c "temp.fastq" > "$output_file"
rm "temp.fastq"

echo "Concatenation and compression complete. Output saved to $output_file."

###############################################################################################

# Specify the input FASTQ.gz files
input_files=("SRR18918166.fastq.gz" "SRR18918167.fastq.gz" "SRR18918168.fastq.gz" "SRR18918169.fastq.gz")


echo "########################################"
echo "###    Concatenation $input_files   ###"   
echo "########################################"

# Specify the output file (compressed FASTQ.gz)
output_file="cSRR18918166.fastq.gz"

# Concatenate and compress the input files
for file in "${input_files[@]}"; do
  zcat "$file" >> "temp.fastq"
done

# Compress the concatenated file and remove the temporary file
gzip -c "temp.fastq" > "$output_file"
rm "temp.fastq"

echo "Concatenation and compression complete. Output saved to $output_file."


###############################################################################################

# Specify the input FASTQ.gz files
input_files=("SRR18918170.fastq.gz" "SRR18918171.fastq.gz" "SRR18918172.fastq.gz" "SRR18918173.fastq.gz")


echo "########################################"
echo "###    Concatenation $input_files   ###"   
echo "########################################"

# Specify the output file (compressed FASTQ.gz)
output_file="cSRR18918170.fastq.gz"

# Concatenate and compress the input files
for file in "${input_files[@]}"; do
  zcat "$file" >> "temp.fastq"
done

# Compress the concatenated file and remove the temporary file
gzip -c "temp.fastq" > "$output_file"
rm "temp.fastq"

echo "Concatenation and compression complete. Output saved to $output_file."
#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status
# Step 1: Align the data Uncomment and modify the following command to align your data
./dragen-os -r /home/tigem/h.poddar/short_reads_pipelines/ref_genome/GRCh38/hash/ -1 
~/mnt/nas-gamba/AshkenazimTrio/P/WGS_150bp/140528_D00360_0018_AH8VC6ADXX_P/140528_D00360_0018_AH8VC6ADXX_R1_fastq.gz -2 
~/mnt/nas-gamba/AshkenazimTrio/P/WGS_150bp/140528_D00360_0018_AH8VC6ADXX_P/140528_D00360_0018_AH8VC6ADXX_R2_fastq.gz
# Step 2: Collate, fixmate, sort, and mark duplicates Uncomment and modify the following command to process your BAM file
samtools collate -@ 32 -O -u sample.bam | samtools fixmate -@ 32 -m -u - - | samtools sort -@ 32 -u - | samtools markdup -f stats.txt -@ 32 -S 
-r - HG002.bam
# Index the final BAM file Uncomment and modify the following command to index your BAM file
samtools index -@ 32 HG002.bam
# Define base directory
BASE="/home/tigem/h.poddar/pipeline_test"
# Define input and output directories
INPUT_DIR="${BASE}/input" OUTPUT_DIR="${BASE}/output"
# Create local directory structure
mkdir -p "${INPUT_DIR}" mkdir -p "${OUTPUT_DIR}"
# Move BAM and BAI files to input directory
mv HG002.bam HG002.bam.bai "${INPUT_DIR}/"
# Create intermediate results directory
INTERMEDIATE_DIRECTORY="${OUTPUT_DIR}/intermediate_results_dir" mkdir -p "${INTERMEDIATE_DIRECTORY}"
# Run DeepVariant
BIN_VERSION="1.6.1" udocker run \ -v "${INPUT_DIR}":"/input" \ -v "${OUTPUT_DIR}":"/output" \ google/deepvariant:"${BIN_VERSION}" \ 
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \ --ref=/home/tigem/h.poddar/short_reads_pipelines/ref_genome/GRCh38/GRCh38.primary_assembly.genome.fa \ 
    --reads=/input/HG002.bam \ --output_vcf=/output/output.vcf.gz \ --output_gvcf=/output/output.vcf.gz \ 
    --intermediate_results_dir=/output/intermediate_results_dir \ --num_shards=20
echo "Pipeline completed successfully!"

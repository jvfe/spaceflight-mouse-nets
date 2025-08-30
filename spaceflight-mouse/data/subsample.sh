#!/bin/bash
set -e

echo "--- Setting up downsampled data directory ---"
mkdir -p sra_data_downsampled

# Seed for reproducibility
SEED=1024
# Fraction of reads to keep
FRACTION=0.1

echo "--- Starting downsampling of FASTQ files ---"

# Loop through all single-end FASTQ files (e.g., *.fastq.gz)
for fastq_file in raw_data/*.fastq.gz; do
  
  # Get the base name of the file, removing the .fastq.gz suffix
  base=$(basename "$fastq_file" .fastq.gz)
  
  echo "Downsampling $base..."
  
  # Downsample the single-end FASTQ file using seqtk
  # The output is compressed with gzip and saved to the new directory
  seqtk sample -s $SEED "$fastq_file" $FRACTION | gzip > sra_data_downsampled/${base}.fastq.gz

done

echo "--- Downsampling complete. Files are in sra_data_downsampled/ ---"

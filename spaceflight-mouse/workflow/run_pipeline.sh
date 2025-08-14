#!/usr/bin/bash
set -e

# Run the pipeline with adjusted paths
nextflow run jvfe/simplernaseq \
    --input data/samplesheet.csv \
    --fasta data/references/chr22.fa \
    --gtf data/references/chr22.gtf \
    --outdir results \
    -profile docker,small_genome \
    -resume
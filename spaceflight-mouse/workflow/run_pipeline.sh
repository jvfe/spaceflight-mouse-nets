#!/usr/bin/bash
set -e

# Run the pipeline with adjusted paths
nextflow run jvfe/simplernaseq \
    --input data/samplesheet.csv \
    --fasta data/references/Mus_musculus.GRCm39.dna.chromosome.19.fa \
    --gtf data/references/Mus_musculus.GRCm39.112.chr19.gtf \
    --outdir results \
    -profile docker,small_genome \
    -resume
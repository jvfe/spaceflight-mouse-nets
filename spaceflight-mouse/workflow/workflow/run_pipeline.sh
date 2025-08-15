#!/usr/bin/bash
set -e

# Run the pipeline with adjusted paths
nextflow run workflow/simplernaseq \
    --input data/samplesheet.csv \
    --fasta data/references/Mus_musculus.GRCm39.dna.chrs_17-19.fa \
    --gtf data/references/Mus_musculus.GRCm39.112.chrs_17-19.gtf \
    --outdir results \
    -profile docker,large_genome \
    -resume

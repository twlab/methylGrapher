#!/bin/bash


# Generate example data
python3 make_example_test.py


input_gfa="./gg.gfa"
work_dir=$(pwd)
genome_prefix="methylgrapher.genome"


# 1. Prepare genome
# It only needs to perform once for each genome.
python3 ../src/main.py PrepareGenome -gfa $input_gfa -lp ./lambda.fa -prefix $genome_prefix


# 2. Run MethylGrapher with one-line command
python3 ../src/main.py simple -fq1 ./R1.fastq -fq2 ./R2.fastq -index_prefixv $genome_prefix -work_dir $work_dir -directional Y


# (OPTIONAL) sorting output
sort -k1n -k2n graph.methyl -o graph.methyl





# Optional
# Clean-up
#rm *.gaf
#rm C2T.R*
#rm G2A.R*




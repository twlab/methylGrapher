#!/bin/bash


# Generate example data
python3 make_example.py


input_gfa="./gg.gfa"
work_dir=$(pwd)
genome_prefix="methylgrapher.genome"


methylGrapher="methylGrapher"

# My own debugging use
# methylGrapher="python3 ../../src/main.py"

# 1. Prepare genome
# It only needs to perform once for each genome.
$methylGrapher PrepareGenome -gfa $input_gfa -lp ./lambda.fa -prefix $genome_prefix


# 2. Alignment
$methylGrapher Align -t 3 -fq1 ./R1.fastq -fq2 ./R2.fastq -index_prefix $genome_prefix -work_dir $work_dir -directional Y

# 3. Methylation Extraction
# -minimum_mapq should be used for real data
$methylGrapher MethylCall -t 3 -index_prefix $genome_prefix -work_dir $work_dir -cg_only Y -genotyping_cytosine Y # -minimum_mapq 10

# 4. Merge CpG
$methylGrapher MergeCpG -index_prefix $genome_prefix -work_dir $work_dir

# 4. Simulate bisulfite conversion rate
$methylGrapher ConversionRate -index_prefix $genome_prefix -work_dir $work_dir


# (OPTIONAL) sorting output
# sort -k1n -k2n graph.methyl -o graph.methyl


# Optional
# Clean-up
#rm *.gaf
#rm C2T.R*
#rm G2A.R*




#!/bin/zsh


# Generate example data
python3 make_example.py


input_gfa="./gg.gfa"
work_dir=$(pwd)
genome_prefix="methylgrapher.genome"


# 1. Prepare genome
# It only needs to perform once for each genome.
python3 ../src/main.py PrepareGenome -gfa $input_gfa -lp ../src/lambda.fa -prefix $genome_prefix


# 2. Prepare Library
# The user (you) needs to run trim glore first.
python3 ../src/main.py PrepareLibrary -fq1 ./R1.fastq -fq2 ./R2.fastq -work_dir $work_dir


# 3. Use VG giraffe for alignment
# Use -non_directional Y if your library is non-directional
python3 ../src/main.py Align -work_dir $work_dir -index_prefix "$genome_prefix"

# 4. Methylation calling
python3 ../src/main.py MethylCall -work_dir $work_dir -index_prefix $genome_prefix

# (OPTIONAL) sorting output, but it is necessary for QC step
sort -k1n -k2n graph.methyl -o graph.methyl

# 5. QC
python3 ../src/main.py QC -work_dir $work_dir -index_prefix $genome_prefix




# Optional
# Clean-up
rm *.gaf
rm C2T.R*
rm G2A.R*




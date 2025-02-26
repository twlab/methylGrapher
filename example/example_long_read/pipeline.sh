#!/bin/bash



python3 make_example.py

# One step to run the whole pipeline
python3 ../../src/mainL.py main -basecall ./test.bam -gfa ./test.gfa -wd ./
echo ""

# Step by step
# Step 1: Prepare Fasta
python3 ../../src/mainL.py prepare_fasta -basecall ./test.bam -wd ./
echo ""

python3 ../../src/mainL.py prepare_fasta -basecall ./test.sam -wd ./
echo ""

# Step 2: Alignment
python3 ../../src/mainL.py align -gfa ./test.gfa -wd ./
echo ""

# Step 3: Methylation extraction
python3 ../../src/mainL.py extraction -gfa ./test.gfa -wd ./
echo ""




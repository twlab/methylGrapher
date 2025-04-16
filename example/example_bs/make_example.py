import os
import sys
import string
import random



test_directory = "./"


def random_h(n):
    return "".join([random.choice(string.ascii_letters) for i in range(n)])

def reverse_complement(seq):
    seq = seq.translate(str.maketrans("ATCG", "TAGC"))
    return seq[::-1]

def lambda_reader(p):
    res = ""
    f = open(p, "r")
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            continue
        res += line
    f.close()
    return res

lambda_genome = lambda_reader("lambda.fa")

genome_graph_segments = {
    "10001": "AAATTAATAGTACGCGTTGCTTGGCCCGCTGGTCATATACTTGAATTTGGTTGAAGTAATCTTGAATTGACTCGAGGTACACTAAACAACCGCCCCCGCT",
    "10002": "TCGGTTGTGACTGTTGAATGACCGAGCACATGCCCTAAGCCACGACCCCGAGCCCGTTGAGGTATCGTTTAAAGATGTCTGCTTCTTCGAATCCCTTCCA",
    "10003": "GCTTCTATATCGGTGGGGGGGACTGTCCTCATTCTCGACTACTCCTTTCTGGCCCGCATCGGAAGACTTTGACTTTGCAACTTGCTCCGTAGAACCTTCT",
    "10004": "TTTTTGCGAGCTACCCTGGACGCAGCTCAGTAGCGGGTTACCCTAAACCCAAAGGCAATCTCTACCTCAAGCCCGCGTCGCTGTCAATGCCGTGTTTGCA",
    "10005": "A",
    "10006": "T",
    "10007": "C",
    "10008": "G",


}

genome_graph_links = {
    "10001": ["10002", "10003", "10005", "10006", "10007", "10008"],
    "10002": ["10004"],
    "10003": ["10004"],
    "10005": ["10004"],
    "10006": ["10004"],
    "10007": ["10004"],
    "10008": ["10004"],
}


p = genome_graph_segments["10001"] + genome_graph_segments["10002"] + genome_graph_segments["10004"]
methylation = {}


reads1 = {}
reads2 = {}

# Simulate 100x~200x coverage, 10001 0% methylation, 10002 50% methylation, 10004 100% methylation

# 10001 - 10002
for i in range(1000):

    read1 = genome_graph_segments["10001"]
    read2 = genome_graph_segments["10002"]
    read2 = reverse_complement(read2)

    read1 = list(read1)
    read2 = list(read2)

    for j in range(len(read1)):
        # 10001 0% methylation
        if read1[j] == "C":
            read1[j] = "T"

    for j in range(len(read2)):
        # 10002 50% methylation
        if read2[j] == "G" and random.random() < 0.5:
            read2[j] = "A"


    rn = random_h(30)
    reads1[rn + "1p2r"] = "".join(read1)
    reads2[rn + "1p2r"] = "".join(read2)


# 10002 - 10004
for i in range(1000):

    read1 = genome_graph_segments["10002"]
    read2 = genome_graph_segments["10004"]
    read2 = reverse_complement(read2)

    read1 = list(read1)
    read2 = list(read2)

    for j in range(len(read1)):
        # 10002 50% methylation
        if read1[j] == "C" and random.random() < 0.5:
            read1[j] = "T"

    rn = random_h(30)
    reads1[rn + "2p4r"] = "".join(read1)
    reads2[rn + "2p4r"] = "".join(read2)



# 10001 - 10003
for i in range(1000):

    s10003 = genome_graph_segments["10003"]
    s10003m = s10003
    # Artificially add C-T in the middle of the read
    if random.random() < 0.5:
        s10003m = s10003[:4] + "T" + s10003[5:]

    read1 = reverse_complement(s10003m)
    read2 = genome_graph_segments["10001"]

    read1 = list(read1)
    read2 = list(read2)

    for j in range(len(read1)):
        # 10003 50% methylation
        if read1[j] == "C" and random.random() < 0.5:
            read1[j] = "T"

    for j in range(len(read2)):
        # 10001 0% methylation
        if read2[j] == "G":
            read2[j] = "A"


    rn = random_h(30)
    reads1[rn + "3f1p"] = "".join(read1)
    reads2[rn + "3f1p"] = "".join(read2)


# print(s10003)
# print(s10003m)



# 10001 - 10006 - 10004
for i in range(1000):

    read1 = genome_graph_segments["10001"]
    read2 = genome_graph_segments["10002"]
    read2 = reverse_complement(read2)

    read1 = list(read1)
    read2 = list(read2)

    for j in range(len(read1)):
        # 10001 0% methylation
        if read1[j] == "C":
            read1[j] = "T"

    for j in range(len(read2)):
        # 10002 0% methylation
        if read2[j] == "G" and random.random() < 0.5:
            read2[j] = "A"


    rn = random_h(30)
    reads1[rn + "1pG2r"] = "".join(read1) + "A"
    reads2[rn + "1pG2r"] = "".join(read2)

# print("AAAA")


# Simulate lambda phage genome
for i in range(2000):

    read1 = lambda_genome[100:200]
    read2 = lambda_genome[250:350]
    read2 = reverse_complement(read2)

    read1 = list(read1)
    read2 = list(read2)

    for j in range(len(read1)):
        if read1[j] == "C" and random.random() < 0.99:
            read1[j] = "T"

    for j in range(len(read2)):
        if read2[j] == "G" and random.random() < 0.99:
            read2[j] = "A"


    rn = random_h(30)
    reads1[rn] = "".join(read1)
    reads2[rn] = "".join(read2)



gg = "H	VN:Z:1.1\n"
for sid in sorted(genome_graph_segments.keys()):
    seg = f"S\t{sid}\t{genome_graph_segments[sid]}\n"
    gg += seg

for i in [2,3,5,6,7,8]:
    gg += f"W	walk{i}	0	chrx	0	201	>10001>1000{i}>10004\n"

for sid_from in sorted(genome_graph_links.keys()):
    for sid_to in genome_graph_links[sid_from]:
        link = f"L\t{sid_from}\t+\t{sid_to}\t+\t0M\n"
        gg += link


fastq1 = ""
for rn, read in reads1.items():
    rid = rn
    quality = "z" * len(read)
    f4l = f"""@{rid}
{read}
+
{quality}
"""
    fastq1 += f4l
    # print(rid, len(read))

fastq2 = ""
for rn, read in reads2.items():
    rid = rn
    quality = "z" * len(read)
    f4l = f"""@{rid}
{read}
+
{quality}
"""
    fastq2 += f4l
    # print(rid, len(read))



if __name__ == "__main__":
    pass

    reference_path = test_directory + "gg.gfa"
    index_prefix = test_directory + "gg"

    output_format = "gaf"
    fin1 = test_directory + "R1.fastq"
    fin2 = test_directory + "R2.fastq"
    fout = test_directory + f"alignment.{output_format}"
    alignment_log = test_directory + f"alignment.log"

    os.makedirs(test_directory, exist_ok=True)
    f = open(reference_path, "w")
    f.write(gg)
    f.close()

    f = open(fin1, "w")
    f.write(fastq1)
    f.close()

    f = open(fin2, "w")
    f.write(fastq2)
    f.close()


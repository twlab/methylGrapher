import os
import sys
import string
import random



test_directory = "./"


def random_h(n):
    return "".join([random.choice(string.ascii_letters) for i in range(n)])

def random_seq(n):
    return "".join([random.choice("ATCG") for i in range(n)])

def reverse_complement(seq):
    seq = seq.translate(str.maketrans("ATCG", "TAGC"))
    return seq[::-1]



genome_graph_segments = {
    "10001": random_seq(600),
    "10002": random_seq(100),
    "10003": random_seq(100),
    "10004": random_seq(600),
}

genome_graph_links = {
    "10001": ["10002", "10003"],
    "10002": ["10004"],
    "10003": ["10004"],
}




gg = "H	VN:Z:1.1\n"
for sid in sorted(genome_graph_segments.keys()):
    seg = f"S\t{sid}\t{genome_graph_segments[sid]}\n"
    gg += seg

for sid_from in sorted(genome_graph_links.keys()):
    for sid_to in genome_graph_links[sid_from]:
        s1s = "+"
        s2s = "+"

        if sid_from == "10003":
            s1s = "-"
        if sid_to == "10003":
            s2s = "-"

        link = f"L\t{sid_from}\t{s1s}\t{sid_to}\t{s2s}\t0M\n"
        gg += link


gg += f"W	walk{2}	0	chrx	0	201	>10001>10002>10004\n"
gg += f"W	walk{3}	0	chrx	0	201	>10001<10003>10004\n"



p1 = genome_graph_segments["10001"] + genome_graph_segments["10002"] + genome_graph_segments["10004"]
p2 = genome_graph_segments["10001"] + reverse_complement(genome_graph_segments["10003"]) + genome_graph_segments["10004"]

sam_file_lines = []

for i in range(100):
    for p in [p1, p2]:
        p = p.upper()

        for j in range(2):
            insert_pos = random.randrange(1, len(p)-1)
            insert_len = random.randrange(1, 10)
            insert_seq = random_seq(insert_len)
            p = p[:insert_pos] + insert_seq + p[insert_pos:]

        for j in range(2):
            delete_pos = random.randrange(10, len(p)-10)
            delete_len = random.randrange(1, 10)
            p = p[:delete_pos] + p[delete_pos+delete_len:]

        skip_cyto = 0
        mm_tag_list = []
        ml_tag_list = []

        for i in range(len(p)-1):
            dinuc = p[i:i+2]

            if dinuc == "CG":
                mm_tag_list.append(skip_cyto)
                ml_tag_list.append(random.randrange(1, 254))
                skip_cyto = 0
                continue

            if p[i] == "C":
                skip_cyto += 1
                continue

        #print(mm_tag_list)
        #print(ml_tag_list)
        assert len(mm_tag_list) == len(ml_tag_list)

        sam_line = []

        sam_line.append(random_h(10))
        sam_line.append("4")
        sam_line.append("*")
        sam_line.append("0")
        sam_line.append("255")
        sam_line.append("*")
        sam_line.append("*")
        sam_line.append("0")
        sam_line.append("0")
        sam_line.append(p)
        sam_line.append("-"*len(p))

        mm_tag = "Mm:Z:C+m," + ",".join([str(x) for x in mm_tag_list])
        ml_tag = "Ml:B:C," + ",".join([str(x) for x in ml_tag_list])

        sam_line.append(mm_tag)
        sam_line.append(ml_tag)

        sam_file_lines.append("\t".join(sam_line))









if __name__ == "__main__":

    reference_path = test_directory + "test.gfa"
    index_prefix = test_directory + "gg"

    alignment_log = test_directory + f"alignment.log"

    os.makedirs(test_directory, exist_ok=True)
    f = open(reference_path, "w")
    f.write(gg)
    f.close()

    f = open(test_directory + "test.sam", "w")
    f.write("\n".join(sam_file_lines))
    f.close()

    os.system(f"samtools view -bS {test_directory}test.sam > {test_directory}test.bam")






import os
import re
import sys
import time

import gfa
import mcall
import utility




class tmpHelp:
    def __init__(self):
        self._help_text = """
        Usage: main.py <command> [options]

        Commands:
        prepareGenome  Prepare genome for methylation calling
        callMethylation  Call methylation from nanopore reads

        Options:
        -h, --help  Show this help message and exit
        -t, --thread  Number of threads to use
        """

    def help_text(self):
        return self._help_text


help = tmpHelp()


utl = utility.Utility()






# Step1: SAM/BAM file to fasta file
# Step2: Alignment
# Step3: Methylation extraction





cytosine_methylation_regex = re.compile(r"\+(\(.*\))\-(\(.*\))")


def path_to_list(path_str):
    # Example: <29983488>29983486<29983485<29983484
    res = []
    direction = True
    segment_ID = ""
    for c in path_str:
        if c in "<>":
            if c == "<":
                direction = False
            elif c == ">":
                direction = True
        else:
            segment_ID += c
            continue

        res.append((segment_ID, direction))
        segment_ID = ""

    res.append((segment_ID, direction))
    res.pop(0)

    return res


def debug_get_fasta_by_read_name(read_name, fasta_fp):
    with open(fasta_fp) as f:
        for l in f:
            if l.startswith(">"):
                if l[1:].strip() == read_name:
                    return f.readline().strip()
    return None


def seq_reverse_complement(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]


gfa_fp = f"/scratch/wzhang/ref/graph/wgbs_bspg/bspg.wl.gfa"


# Store entire GFA in memory
class GraphicalFragmentAssemblyMemory(object):

    def __init__(self):
        self.clear()

    def clear(self):
        self._header = ""
        self._segment = {}
        self._walk = {}

        self._original_gfa_path = ""

    def parse(self, gfa_file):
        self._original_gfa_path = gfa_file
        # TODO when not keep_tag, keep_link, forbid to write GFA and access those two
        for l in open(gfa_file):
            if l[0] not in "HSL":
                continue

            l = l.strip().split("\t")

            if l[0] == "H":
                self._header = "\t".join(l[1:])

            if l[0] == "S":
                rt, segID, seq, *tags = l
                # sequence, tag, links, parent_links_count
                self._segment[segID] = seq

            else:
                continue

    def get_sequence_by_segment_ID(self, segment_ID):
        return self._segment[segment_ID]

    def get_sequence_length_by_segment_ID(self, segment_ID):
        return len(self._segment[segment_ID])

    def get_sequences_by_segment_ID(self, segment_IDs):
        res = {}
        for sID in segment_IDs:
            res[sID] = self.get_sequence_by_segment_ID(sID)
        return res

    def get_tag_by_segment_ID(self, segment_ID):
        raise NotImplementedError

    def get_parent_link_count(self, segment_ID):
        return self._segment[segment_ID][3]

    def get_parent_links(self, segment_ID):
        raise NotImplementedError

    def get_child_links(self, segment_ID):
        return self._segment[segment_ID][2]

    def all_segment_IDs(self):
        return self._segment.keys()


gfa_instance = GraphicalFragmentAssemblyMemory()
gfa_instance.parse(gfa_fp)













def bsam_to_fasta(bsam_fp, fasta_fp):
    # Explore sam file

    output_fasta = open(fasta_fp, "w")

    no_methylation_tag = 0
    with open(bsam_fp) as f:
        read_line_count = 0
        for l in f:
            # Skip header
            if l.startswith('@'):
                continue

            line = l.strip().split('\t')
            basic_info = line[:11]
            tags_list = line[11:]
            tags = {}
            for tag_str in tags_list:
                tag = [tag_str[:4], tag_str[5:]]
                if tag_str.startswith("Ml:B:C"):
                    tag = ["Ml:B:C", tag_str[7:]]
                tags[tag[0]] = tag[1]

            read_name = basic_info[0]
            read_seq = basic_info[9].upper()

            if 'Mm:Z' not in tags:
                no_methylation_tag += 1
                continue

            if "Ml:B:C" not in tags:
                no_methylation_tag += 1
                continue

            mm_tag = tags['Mm:Z']
            ml_tag = tags["Ml:B:C"]

            if mm_tag.endswith(";"):
                mm_tag = mm_tag[:-1]
            if ml_tag.endswith(";"):
                ml_tag = ml_tag[:-1]

            mm_tag = mm_tag.split(",")
            ml_tag = ml_tag.split(",")

            if mm_tag[0] != 'C+m':
                # print(mm_tag[0], "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                continue
            else:
                mm_tag = mm_tag[1:]

            if len(ml_tag) in [0, 1]:
                continue

            mm_tag = [int(x) for x in mm_tag]
            ml_tag = [int(x) for x in ml_tag]
            assert len(mm_tag) == len(ml_tag)

            cg_count = read_seq.count("C")

            cytosines_positions = []
            cpg_positions = []
            for i, nuc in enumerate(read_seq):
                if nuc == "C":
                    cytosines_positions.append(i)
                    if i < len(read_seq) - 1 and read_seq[i + 1] == "G":
                        cpg_positions.append(i)

            methylated_cytosines = []
            unmethylated_cytosines = []
            nth_cytosine = 0
            for i in range(len(mm_tag)):
                interval = mm_tag[i] + 1
                if nth_cytosine == 0:
                    nth_cytosine = interval - 1
                else:
                    nth_cytosine += interval
                # print(nth_cytosine, interval)

                methylated = ml_tag[i] > 127
                if methylated:
                    methylated_cytosines.append(cytosines_positions[nth_cytosine])
                else:
                    unmethylated_cytosines.append(cytosines_positions[nth_cytosine])

            cpg_covered_count = len(
                set(cpg_positions).intersection(set(methylated_cytosines).union(set(unmethylated_cytosines))))
            if cpg_covered_count < len(cpg_positions) - 5:
                print(read_line_count, cpg_covered_count, len(cpg_positions))
                continue
            # print(read_name)

            # print(read_line_count)
            # print(read_name)
            # print(l)
            # print(len(read_seq))
            # print(line)
            # print(basic_info[:9])
            # print(tags)
            # print(cg_count)
            # print(len(mm_tag))
            # print(len(ml_tag))

            # print(mm_tag[:40])
            # print(cytosines_positions[:80])
            # print()
            # print(methylated_cytosines[:20])
            # print(cpg_positions[:20])

            # print(len(methylated_cytosines), len(cpg_positions))
            # print(methylated_cytosines[-1], cpg_positions[-1])

            # print('\n')

            methylated_cytosines = [str(x) for x in methylated_cytosines]
            unmethylated_cytosines = [str(x) for x in unmethylated_cytosines]
            met_str = f"+({','.join(methylated_cytosines)})-({','.join(unmethylated_cytosines)})"
            fasta_entry = f">{read_name}{met_str}\n{read_seq}\n\n"

            # print(fasta_entry)
            output_fasta.write(fasta_entry)

            read_line_count += 1
            if read_line_count > 10:
                pass

    output_fasta.close()
    return None





counter_pass = 0
counter_issue_no_methylation = 0
counter_issue = 0
counter_total = 0

methylation_result = {}

time_start = time.time()
with open(output_fp) as f:
    for l in f:

        l = l.strip().split("\t")
        # print(l)

        counter_total += 1
        if counter_total % 1000 == 0:
            duration = time.time() - time_start
            print(counter_total, counter_pass, counter_issue_no_methylation, counter_issue, duration)

        read_name = l[0]

        query_len = int(l[1])
        query_start = int(l[2])
        query_end = int(l[3])

        strand = l[4]

        path_str = l[5]
        # if ">" in path_str:
        #    continue

        path_len = int(l[6])
        path_start = int(l[7])
        path_end = int(l[8])

        matched_len = int(l[9])
        alignment_block_len = int(l[10])
        mapq = int(l[11])

        tag_list = l[12:]
        tags = {}

        # print(tag_list)

        for tag_str in tag_list:
            tname, tdt, tv = tag_str.split(":")
            if tdt == "i":
                tv = int(tv)
            elif tdt == "f":
                tv = float(tv)
            tags[tname] = tv
            # print(tname, tv)

        mc = cytosine_methylation_regex.findall(read_name)
        if len(mc) != 1:
            counter_issue_no_methylation += 1
            continue

        methylated_cytosines = mc[0][0][1:-1].split(",")
        unmethylated_cytosines = mc[0][1][1:-1].split(",")

        if len(methylated_cytosines) + len(unmethylated_cytosines) == 0:
            counter_issue_no_methylation += 1
            continue

        methylated = list(map(int, list(filter(lambda x: x != "", methylated_cytosines))))
        unmethylated = list(map(int, list(filter(lambda x: x != "", unmethylated_cytosines))))

        # Example: 42M1I5M1D23M1I1M1I1M1I13M1D41M1D35M1D4M1D
        cigar_str = tags["cg"]
        cigar = []
        for c in re.findall(r"\d+[MID]", cigar_str):
            cigar.append((int(c[:-1]), c[-1]))

        alignment_block_type = []
        block_start = 0
        block_end = 0
        for length, alignment_type in cigar:
            block_end += length

            alignment_block_type.append((block_start, block_end, alignment_type))
            block_start = block_end
        alignment_block_type.append((block_start, block_end, alignment_type))

        """     
        TODO uncomment
        cigar_freq = {}
        for c in cigar:
            if c[1] not in cigar_freq:
                cigar_freq[c[1]] = 0
            cigar_freq[c[1]] += c[0]

        assert path_end - path_start == cigar_freq["M"] + cigar_freq.get("D", 0)
        assert query_end - query_start == cigar_freq["M"] + cigar_freq.get("I", 0)


        # Just to make sure the cigar interpretation is correct
        reconstructed_cigar = ""
        for c in cigar:
            reconstructed_cigar += str(c[0]) + c[1]
        assert reconstructed_cigar == cigar_str
        """

        # Just to make sure the path interpretation is correct
        path_len_constructed = 0
        path_list = path_to_list(path_str)
        for segID, direction in path_list:
            seq = gfa_instance.get_sequence_by_segment_ID(segID)
            # print(segID, len(seq))
            path_len_constructed += len(seq)
        # print(path_str)
        assert path_len_constructed == path_len

        # query_seq = debug_get_fasta_by_read_name(read_name, "./short.fasta")
        query_seq = ""
        path_seq = ""

        """
        for i, (segID, direction) in enumerate(path_list):
            seq = gfa_instance.get_sequence_by_segment_ID(segID)
            if not direction:
                seq = seq_reverse_complement(seq)
            path_seq += seq
        """

        path_segment_length = {}
        for i, (segID, direction) in enumerate(path_list):
            seql = gfa_instance.get_sequence_length_by_segment_ID(segID)
            path_segment_length[segID] = seql

        pos_in_read = -1
        pos_in_path = -1
        for alignment_block_len, alignment_block_type in cigar:

            for ai in range(alignment_block_len):
                if alignment_block_type == "M":
                    pos_in_read += 1
                    pos_in_path += 1

                    pos_in_abs_path = path_start + pos_in_path
                    pos_in_abs_read = query_start + pos_in_read

                    met_flag = None
                    if pos_in_abs_read in methylated:
                        met_flag = "M"
                    elif pos_in_abs_read in unmethylated:
                        met_flag = "U"

                    if met_flag is not None:

                        corresponding_segment_ID = None
                        corresponding_segment_pos = pos_in_abs_path
                        for segID, direction in path_list:
                            if corresponding_segment_pos < path_segment_length[segID]:
                                corresponding_segment_ID = segID
                                if not direction:
                                    # TODO need double check
                                    corresponding_segment_pos = path_segment_length[
                                                                    segID] - corresponding_segment_pos - 1
                                break

                            corresponding_segment_pos -= path_segment_length[segID]

                        segment_seq = gfa_instance.get_sequence_by_segment_ID(corresponding_segment_ID)
                        csp1 = corresponding_segment_pos
                        csp2 = corresponding_segment_pos + 2
                        if csp2 > len(segment_seq):
                            csp2 = csp1 + 1

                        if not direction:
                            csp2 = corresponding_segment_pos
                            csp1 = corresponding_segment_pos - 2
                            if csp1 < 0:
                                csp1 = 0
                        # print(path_seq[pos_in_abs_path: pos_in_abs_path+2], query_seq[pos_in_abs_read:pos_in_abs_read+2], segment_seq[csp1: csp2], met_flag)

                        if corresponding_segment_ID not in methylation_result:
                            methylation_result[corresponding_segment_ID] = {}
                        if corresponding_segment_pos not in methylation_result[corresponding_segment_ID]:
                            dir_str = "-" if direction else "+"
                            methylation_result[corresponding_segment_ID][corresponding_segment_pos] = [dir_str, 0, 0]
                        methylation_result[corresponding_segment_ID][corresponding_segment_pos][2] += 1
                        if met_flag == "M":
                            methylation_result[corresponding_segment_ID][corresponding_segment_pos][1] += 1




                elif alignment_block_type == "I":
                    pos_in_read += 1
                elif alignment_block_type == "D":
                    pos_in_path += 1
                else:
                    raise Exception("Unknown alignment block type")

        # print(path_len_constructed, path_len)
        # print(cigar)
        # print(alignment_block_type)

        # print(path_start, path_end, path_end - path_start)
        # print(query_start, query_end, query_end-query_start)
        # print(cigar_freq)

        # print(path_seq[433:533])
        # print(query_seq[:100])
        # print(path_str)
        counter_pass += 1
        # break




if __name__ == "__main__":
    args = sys.argv
    args.pop(0)

    ga_path = "GraphAligner"
    thread = 1

    # Find config file
    try:
        config = utility.ConfigParser("config.ini")
        vg_path = config.get("default", "ga_path")
        thread = config.get("default", "thread")
    except:
        pass

    try:
        thread = int(thread)
    except:
        thread = 1

    if len(args) == 0:
        print("No command specified. Use 'help' for more information.")
        print(help.help_text())
        sys.exit(0)

    command = args.pop(0)
    command = command.lower()
    if command not in ["help", "-h", "--help", "main"]:
        print(f"Unknown command: {command}")
        sys.exit(1)

    kvargs = {}
    while len(args) > 0:
        arg = args.pop(0)
        if arg.startswith("-"):
            key = arg[1:]
            value = args.pop(0)
            kvargs[key] = value
        else:
            print(f"Unknown argument: {arg}\n\n")
            print(help.help_text())
            sys.exit(1)

    if "thread" in kvargs:
        thread = int(kvargs["thread"])
        del kvargs["thread"]
    if "t" in kvargs:
        thread = int(kvargs["t"])
        del kvargs["t"]


    if command in ["help", "-h", "--help"]:
        print(help.help_text())
        sys.exit(0)



    if command == "preparegenome":
        original_gfa_file_path = kvargs["gfa"]
        prefix = kvargs["prefix"]

        fn_report = prefix + ".prepare.genome.report.txt"
        freport = open(fn_report, "w")


import os
import re
import sys
import time
import gzip
import argparse

import gfa
import utility








utl = utility.Utility()

cytosine_methylation_regex = re.compile(r"\+(\(.*\))\-(\(.*\))")
cigar_regex = re.compile(r"\d+[MID]")


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

def sam_bam_line_reader(file_path):
    lower_fp = file_path.lower()
    if lower_fp.endswith(".sam"):
        with open(file_path, 'r') as f:
            for line in f:
                yield line
    elif lower_fp.endswith(".bam"):
        samtool_cmd = utility.SystemExecute()
        samtool_out, samtool_err = samtool_cmd.execute(f"samtools view -h {file_path}", )
        for line in samtool_out:
            yield str(line, 'utf-8')
        samtool_cmd.wait()



def bam_sam_to_fasta(bsam_fp, fasta_fp):
    # Explore sam file

    output_fasta_fh = None
    if utl.isGzip(fasta_fp):
        output_fasta_fh = gzip.open(fasta_fp, 'wt')
    else:
        output_fasta_fh = open(fasta_fp, 'w')


    no_methylation_tag = 0
    read_line_count = 0


    durations = [0] * 10
    timestamps = [0] * 10
    next_report = 10


    for l in sam_bam_line_reader(bsam_fp):
        # Skip header
        if l.startswith('@'):
            continue


        if timestamps[0] == 0:
            ts = time.time()
            for i in range(10):
                timestamps[i] = ts

        durations[0] += time.time() - timestamps[0]

        if sum(durations) > next_report:
            print(durations)
            next_report += 10

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
            # print(read_line_count, cpg_covered_count, len(cpg_positions))
            # continue
            pass

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
        output_fasta_fh.write(fasta_entry)

        read_line_count += 1


        tnow = time.time()
        durations[1] += tnow - timestamps[0]
        timestamps[0] = tnow




    output_fasta_fh.close()
    return None

























































def generate_fp_within_wd(wd, file_name):
    return os.path.join(wd, file_name)

def generate_fasta_fp(wd):
    return generate_fp_within_wd(wd, "input.fasta")


def extract_methylation_single_thread(gfa_fp, wd, threads=1, debug=False):
    """
    Extract methylation information
    :param wd: working directory
    :param gfa: genome graph file path in GFA format
    :return: methylation information
    """

    counter_pass = 0
    counter_issue_no_methylation = 0
    counter_issue = 0
    counter_total = 0

    methylation_result = {}

    time_start = time.time()

    gfa_instance = gfa.GraphicalFragmentAssemblySegmentLengthMemory()
    gfa_instance.parse(gfa_fp)

    if debug:
        # Just in case the sequence is needed for debugging
        gfa_seq_instance = gfa.GraphicalFragmentAssemblyMemorySegmentOptimized()
        gfa_seq_instance.parse(gfa_fp)



    with open(generate_fp_within_wd(wd, "align.gaf")) as f:
        for l in f:

            l = l.strip().split("\t")
            # print(l)

            counter_total += 1
            if counter_total % 1000 == 0:
                duration = time.time() - time_start
                # TODO better logging
                # print(counter_total, counter_pass, counter_issue_no_methylation, counter_issue, duration)

            # Read the graph alignment file (GAF)
            read_name = l[0]

            query_len = int(l[1])
            query_start = int(l[2])
            query_end = int(l[3])

            strand = l[4]

            path_str = l[5]
            path_list = path_to_list(path_str)

            path_len = int(l[6])
            path_start = int(l[7])
            path_end = int(l[8])

            matched_len = int(l[9])
            alignment_block_len = int(l[10])
            mapq = int(l[11])

            tag_list = l[12:]
            tags = {}
            for tag_str in tag_list:
                tname, tdt, tv = tag_str.split(":")
                if tdt == "i":
                    tv = int(tv)
                elif tdt == "f":
                    tv = float(tv)
                tags[tname] = tv
                # print(tname, tv)

            # Example: 42M1I5M1D23M1I1M1I1M1I13M1D41M1D35M1D4M1D
            cigar = []
            cigar_str = tags["cg"]
            for c in cigar_regex.findall(cigar_str):
                cigar.append((int(c[:-1]), c[-1]))

            alignment_block_type = []
            block_start = 0
            block_end = 0
            for length, alignment_type in cigar:
                block_end += length

                alignment_block_type.append((block_start, block_end, alignment_type))
                block_start = block_end
            alignment_block_type.append((block_start, block_end, alignment_type))



            # Sanity checks

            # make sure there is methylation information here, and there is no conflict
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


            # Just to make sure the cigar interpretation is correct
            cigar_freq = {}
            for c in cigar:
                if c[1] not in cigar_freq:
                    cigar_freq[c[1]] = 0
                cigar_freq[c[1]] += c[0]
            assert path_end - path_start == cigar_freq["M"] + cigar_freq.get("D", 0)
            assert query_end - query_start == cigar_freq["M"] + cigar_freq.get("I", 0)

            reconstructed_cigar = ""
            for c in cigar:
                reconstructed_cigar += str(c[0]) + c[1]
            assert reconstructed_cigar == cigar_str


            # if ">" in path_str:
            #    continue


            # Just to make sure the path interpretation is correct
            path_len_constructed = 0
            path_segment_length = {}
            for segID, direction in path_list:
                seq_length = gfa_instance.get_sequence_length_by_segment_ID(segID)
                path_len_constructed += seq_length
                path_segment_length[segID] = seq_length
            assert path_len_constructed == path_len


            # Debugging
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
                                        corresponding_segment_pos = path_segment_length[segID] - corresponding_segment_pos - 1
                                    break

                                corresponding_segment_pos -= path_segment_length[segID]

                            if debug:
                                segment_seq = gfa_seq_instance.get_sequence_by_segment_ID(corresponding_segment_ID)
                                csp1 = corresponding_segment_pos
                                csp2 = corresponding_segment_pos + 2
                                if csp2 > len(segment_seq):
                                    csp2 = csp1 + 1

                                if not direction:
                                    csp2 = corresponding_segment_pos
                                    csp1 = corresponding_segment_pos - 2
                                    if csp1 < 0:
                                        csp1 = 0

                                # This is to check whether I am aligning CG within the read to CG in the segment
                                print(path_seq[pos_in_abs_path: pos_in_abs_path+2], query_seq[pos_in_abs_read:pos_in_abs_read+2], segment_seq[csp1: csp2], met_flag)

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



    for segment_ID in sorted(methylation_result):
        for pos in sorted(methylation_result[segment_ID]):
            strand, met, cov = methylation_result[segment_ID][pos]
            unmet = cov - met
            l = f"{segment_ID}\t{pos}\t{strand}\t{met}\t{unmet}\t{cov}"
            print(l)
    return None













































def prepare_fasta(basecall, wd, threads=1, debug=False):
    """
    Prepare fasta file from long read base call BAM/SAM file
    :param basecall: long read base call BAM/SAM file
    :param wd: working directory
    :return: fasta file
    """

    bam_sam_to_fasta(basecall, generate_fasta_fp(wd))
    return None


def align(gfa, wd, threads=1, debug=False):
    """
    Align long reads to genome graph
    :param wd: working directory
    :param gfa: genome graph file path in GFA format
    """

    # TODO use non-default version of GraphAligner
    cmd = f"GraphAligner -t {threads} --cigar-match-mismatch -g {gfa} -f {generate_fasta_fp(wd)} -a {generate_fp_within_wd(wd, 'align.gaf')} -x vg"
    align_exe = utility.SystemExecute()

    # TODO
    # align_exe.execute(cmd, stdout=generate_fp_within_wd(wd, 'align.log'), stderr=generate_fp_within_wd(wd, 'align.err'))
    print(cmd)

    return None

def extract_methylation(gfa_fp, wd, threads=1, debug=False):
    """
    Extract methylation information
    :param wd: working directory
    :param gfa_fp: genome graph file path in GFA format
    :return: methylation information
    """

    extract_methylation_single_thread(gfa_fp, wd, debug=debug)

    return None


def main(basecall, gfa_fp, wd, threads=1, debug=False):
    """
    One step to run all
    :param basecall: long read base call BAM/SAM file
    :param gfa_fp: genome graph file path in GFA format
    :param wd: working directory
    :return: methylation information
    """

    prepare_fasta(basecall, wd, threads=threads, debug=debug)
    align(gfa_fp, wd, threads=threads, debug=debug)
    extract_methylation(gfa_fp, wd, threads=threads, debug=debug)

    return None





if __name__ == "__main__":

    pass

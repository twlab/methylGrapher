

import os
import re
import sys
import time

import gfa
import utility



alignment_csz_tag_start_pattern = re.compile(r"^(\+|-)([ACGTN]*)")
alignment_csz_tag_end_pattern = re.compile(r"(\+|-)([ACGTN]*)$")

alignment_csz_tag_match = re.compile(r":\d*")
alignment_csz_tag_mismatch = re.compile(r"\*\w\w")

alignment_csz_tag_deletion = re.compile(r"-\w")
alignment_csz_tag_insertion = re.compile(r"\+\w")




def alignment_path_parse(path):
    res = [[], []]
    element = ""
    for i in path:
        if i in "><":

            if i == ">":
                res[1].append(1)
            else:
                res[1].append(-1)

            if element != "":
                res[0].append(element)
                element = ""

        else:
            element += i

    res[0].append(element)
    assert len(res[0]) == len(res[1])
    assert "" not in res[0]
    return res


def cs_tag_parse(alignment_tag):
    query_start_offset = 0
    query_end_offset = 0
    ref_len = 0

    tag_parsed = []
    indels = []

    element = ""
    for s in alignment_tag:
        if s in "+-:*":
            if len(element) == 0 and len(tag_parsed) == 0:
                element = s
                continue
            tag_parsed.append(element)
            element = s
            continue
        element += s
    tag_parsed.append(element)

    #print(alignment_tag)
    #print(tag_parsed)
    #print()


    if tag_parsed[0].startswith("+"):
        s = len(tag_parsed[0])-1
        query_start_offset = s
        tag_parsed.pop(0)
    if tag_parsed[-1].startswith("-"):
        e = len(tag_parsed[-1])-1
        query_end_offset = e
        tag_parsed.pop(-1)

    for i in tag_parsed:
        if i.startswith(":"):
            ref_len += int(i[1:])

        elif i.startswith("*"):
            ref_len += 1

        elif i.startswith("-"):
            l = len(i[1:])
            indel = (ref_len, l, "-")
            indels.append(indel)
            ref_len += l

        elif i.startswith("+"):
            l = len(i[1:])
            indel = (ref_len, l, "+")
            indels.append(indel)
            # ref_len -= l

    # print(alignment_csz_tag_mismatch.findall(alignment_tag))
    #ref_len += len(alignment_csz_tag_mismatch.findall(alignment_tag))
    #ref_len -= len(alignment_csz_tag_deletion.findall(alignment_tag))

    return query_start_offset, query_end_offset, ref_len, indels




def call(gfa, alignment, methylation_output, minimum_identity=50, minimum_mapq=20):

    # Alignment count, low confidence count, error count
    counter1 = 0
    counter2 = 0
    counter3 = 0

    utl = utility.Utility()

    alignment_file_handles = []
    if isinstance(alignment, list) or isinstance(alignment, tuple) or isinstance(alignment, set):
        for a in alignment:
            alignment_file_handles.append(open(a))
    else:
        alignment_file_handles.append(open(alignment))



    for alignment_file_handle in alignment_file_handles:
        lines_for_same_read = []
        for l in alignment_file_handle:
            ts = time.time()
            l = l.strip().split("\t")

            asterisk = False
            for i in [2, 3, 6, 7, 8, 9, 10, 11]:
                if l[i] == "*":
                    asterisk = True
                    break
                l[i] = int(l[i])
            if asterisk:
                continue


            if len(lines_for_same_read) == 0:
                lines_for_same_read.append(l)
                continue

            # @A00584:440:HJLVHDSX2:2:1101:3007:1031C2T1R0
            read_name = l[0][:-6]
            conversion_type = l[0][-6:-3]
            read_index = l[0][-2:]
            read_name = read_name + read_index

            l[0] = read_name
            l.append("ct:Z:" + conversion_type)

            if l[0] == lines_for_same_read[0][0]:
                lines_for_same_read.append(l)
                continue
            else:
                # Figure out the best alignment for the read pair
                best_alignments = {}

                for a in lines_for_same_read:
                    read_seq = ""
                    for tag in a[12:]:
                        if tag.startswith("bq:Z:"):
                            read_seq = tag[5:]
                            break

                    if read_seq not in best_alignments or a[9] > best_alignments[read_seq][9]:
                        best_alignments[read_seq] = a

                lines_for_same_read = [l]

            # print(best_alignments)
            for best_alignment in best_alignments.values():
                counter1 += 1

                low_confidence = False
                if best_alignment[11] < minimum_mapq:
                    low_confidence = True
                if best_alignment[9] < minimum_identity:
                    low_confidence = True
                if low_confidence:
                    counter2 += 1
                    continue


                l = best_alignment



                duration1 = time.time() - ts

                # VG giraffe always outputs + here, just in case
                assert l[4] == "+"

                query_start = l[2]
                query_end = l[3]

                # VG giraffe does not provide stable path_end
                path = alignment_path_parse(l[5])
                path_start = l[7]
                path_end = 0
                """
                path_strand = "+"
                tmp1 = path.startswith(">")
                tmp2 = path.startswith("<")
                assert not (tmp1 and tmp2)
                if tmp2:
                    path_strand = "-"
                """

                # Figure out the alignment start and end, but to do that, I need the alignment CS tag
                alignment_tag = None
                original_bs_read = None
                read_conversion_type = None
                for tag in l[12:]:
                    if tag.startswith("cs:Z:"):
                        alignment_tag = tag[5:]
                    if tag.startswith("bq:Z:"):
                        original_bs_read = tag[5:]
                    if tag.startswith("ct:Z:"):
                        rct1 = tag[5]
                        rct2 = tag[7]
                        read_conversion_type = rct1


                if alignment_tag is None or original_bs_read is None or read_conversion_type is None:
                    continue

                qso, qeo, rl, indels = cs_tag_parse(alignment_tag)
                query_start += qso
                query_end -= qeo
                path_end = path_start + rl

                duration2 = time.time() - ts


                path_sequences = []
                for pi in range(len(path[0])):
                    segmentID = path[0][pi]
                    seq = gfa.get_sequence_by_segment_ID(segmentID)
                    if path[1][pi] < 0:
                        seq = utl.seq_reverse_complement(seq)
                    path_sequences.append(seq)
                path_sequence = "".join(path_sequences)

                duration3 = time.time() - ts

                """
                TODO remove
                if path_strand == "-":
                    a = len(path_sequence) - path_end
                    b = len(path_sequence) - path_start
                    path_start = a
                    path_end = b
                """

                path_seq_portion = path_sequence[path_start:path_end]
                bs_read_portion = original_bs_read[query_start:query_end]


                # Reconstruct the alignment
                for indel in indels:
                    indel_start = indel[0]
                    indel_leng  = indel[1]
                    if indel[2] == "-":
                        bs_read_portion = bs_read_portion[:indel_start] + " "*indel_leng + bs_read_portion[indel_start:]
                    else:
                        bs_read_portion = bs_read_portion[:indel_start] + bs_read_portion[indel_start+indel_leng:]

                if len(path_seq_portion) != rl or len(path_seq_portion) != len(bs_read_portion):
                    counter3 += 1

                    continue
                    print(
                        "\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n")
                    print("Error rate:", (counter3 / counter1) * 100, "%")
                    time.sleep(0.3)
                    continue

                duration4 = time.time() - ts

                if False:
                    # TODO Just for debugging
                    identity = 0

                    #print(len(path_seq_portion), path, path_start, path_end, rl)
                    #print(path_sequences)
                    for i in range(len(path_seq_portion)):
                        if path_seq_portion[i] == bs_read_portion[i]:
                            identity += 1
                    identity = "\033[92m%0.0f%%\033[0m" % (identity / len(path_seq_portion)*100)

                    print("Identity:", identity, path)
                    #print(path_start, path_end)
                    #print(path_sequence)
                    #print(query_start, query_end)
                    #print(original_bs_read)
                    #print()
                    #path_sequence[path_start:path_end]
                    endx = min(path_end+20, len(path_sequence))
                    alignment_tag_highlight = ""
                    for i in range(len(alignment_tag)):
                        if alignment_tag[i] == "+":
                            alignment_tag_highlight += "\033[92m+\033[0m"
                        elif alignment_tag[i] == "-":
                            alignment_tag_highlight += "\033[91m-\033[0m"
                        else:
                            alignment_tag_highlight += " "
                    print(" "*(len(path_seq_portion)-5) + path_sequence[path_end-5:endx])
                    print(path_seq_portion)
                    print(bs_read_portion)

                    xxx = ""
                    for i in range(len(path_seq_portion)):
                        if path_seq_portion[i] == bs_read_portion[i]:
                            xxx += " "
                        else:
                            a = path_seq_portion[i]
                            b = bs_read_portion[i]
                            if a == "C" and b == "T":
                                xxx += "\033[92mT\033[0m"
                            elif a == "G" and b == "A":
                                xxx += "\033[92mA\033[0m"
                            else:
                                xxx += "\033[91mX\033[0m"
                    tcount = xxx.count("T")
                    acount = xxx.count("A")
                    conversion_rate = 0
                    if read_conversion_type == "C":
                        conversion_rate = tcount / path_seq_portion.count("C")
                    else:
                        conversion_rate = acount / path_seq_portion.count("G")
                    conversion_rate = conversion_rate * 100
                    print(xxx)
                    print(f"Conversion rate: {conversion_rate:.2f}")
                    #print(original_bs_read)
                    #print(alignment_tag)
                    #print(alignment_tag_highlight)
                    print()

                duration5 = time.time() - ts

                #l = len(path_seq_portion)
                pl = len(path_sequence)

                segment_index = 0
                segmentID = path[0][segment_index]
                segment_orientation = path[1][segment_index]
                segment_length = len(path_sequences[segment_index])
                path_len_so_far = segment_length
                last_segmentID = None
                #print(path)
                #print(path_sequences)
                #print(path_original_sequence)

                for i, ref_base in enumerate(path_seq_portion):
                    path_pos = i + path_start
                    # ref_base = path_seq_portion[i]
                    read_base = bs_read_portion[i]
                    category = "U"

                    if ref_base == "C"   and read_conversion_type == "C":
                        if read_base not in "CT":
                            continue
                    elif ref_base == "G" and read_conversion_type == "G":
                        if read_base not in "GA":
                            continue
                    else:
                        continue


                    while path_pos >= path_len_so_far:
                        segment_index += 1
                        segmentID = path[0][segment_index]
                        segment_length = len(path_sequences[segment_index])
                        segment_orientation = path[1][segment_index]
                        path_len_so_far += segment_length

                    if segment_orientation > 0:
                        segment_pos = segment_length - (path_len_so_far - path_pos)
                    else:
                        segment_pos = path_len_so_far - path_pos - 1

                    """
                    x = path_sequence[path_pos]
                    y = path_seq_portion[i]
                    z = path_original_sequence[segmentID][segment_pos]
                    if segment_orientation < 0:
                        z = utl.seq_reverse_complement(z)
                    assert x == y and y == z
                    """
                    base_strand = segment_orientation
                    if ref_base == "C":
                        if path_pos+1 < pl:
                            ref_base2 = path_sequence[path_pos+1]
                            if ref_base2 == "G":
                                category = "CG"
                            else:
                                if path_pos+2 < pl:
                                    ref_base3 = path_sequence[path_pos+2]
                                    if ref_base3 == "G":
                                        category = "CHG"
                                    else:
                                        category = "CHH"

                        if read_base == "T":
                            methylated = 0
                        elif read_base == "C":
                            methylated = 1
                    elif ref_base == "G":
                        base_strand = -base_strand
                        if path_pos > 0:
                            ref_base2 = path_sequence[path_pos-1]
                            if ref_base2 == "C":
                                category = "CG"
                            else:
                                if path_pos > 1:
                                    ref_base3 = path_sequence[path_pos-2]
                                    if ref_base3 == "C":
                                        category = "CHG"
                                    else:
                                        category = "CHH"

                        if read_base == "G":
                            methylated = 1
                        elif read_base == "A":
                            methylated = 0
                    else:
                        continue


                    if base_strand < 0:
                        base_strand = "-"
                    else:
                        base_strand = "+"


                    mcall_tmp = methylation_output[int(segmentID[-1])]

                    if segmentID != last_segmentID:
                        mcall_tmp.write(f"{segmentID}\t{segment_pos}\t{base_strand}\t{category}\t{methylated}\n")
                        last_segmentID = segmentID
                    else:
                        pass
                        mcall_tmp.write(f"\t{segment_pos}\t{base_strand}\t{category}\t{methylated}\n")

                    #print(i, category, ref_base, read_base, methylated)
                    #print(path_seq_portion[i], bs_read_portion[i], )


                duration = time.time() - ts
                #print(f"Total time: {duration:.7f}s")
                # print("Low confident rate: %0.2f%%, Error rate: %0.2f%%" % (100*counter2/counter1, 100*counter3/counter1))
                # print("Total time: %0.8fs" % (duration))

                #duration6 = duration - duration5
                #duration5 = duration5 - duration4
                #duration4 = duration4 - duration3
                #duration3 = duration3 - duration2
                #duration2 = duration2 - duration1
                #print(f"{int(duration1/duration*100)} {int(duration2/duration*100)} {int(duration3/duration*100)} {int(duration4/duration*100)} {int(duration5/duration*100)}  {int(duration6/duration*100)} ")
                #print()

    # Low confident count, Error count
    return counter1, counter2, counter3



def mcall_tmp_to_final(tmp_fh, out_fh):
    segmentID = ""
    methylation = {}
    for l in tmp_fh:
        l = l.strip().split("\t")
        if len(l) == 5:
            segmentID = l.pop(0)
        key = (segmentID, l[0], l[1], l[2])
        if key not in methylation:
            methylation[key] = [0, 0]
        methylation[key][int(l[3])] += 1

    for key in methylation.keys():
        unmethylated = methylation[key][0]
        methylated = methylation[key][1]
        coverage = unmethylated + methylated
        mlevel = methylated / coverage
        out_fh.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{key[3]}\t{unmethylated}\t{methylated}\t{coverage}\t{mlevel}\n")


if __name__ == "__main__":

    if len(sys.argv) > 2:
        if sys.argv[1] != "call":
            sys.exit(3)
        alignment_fp = []
        for i in range(10):
            alignment_fp.append(f"alignment.{i}.sorted.gaf")
        gfa_fp = sys.argv[2]
        out_fp = sys.argv[3]

        print("start")

        gfa_instance = gfa.GraphicalFragmentAssemblyMemory()
        gfa_instance.parse(gfa_fp)

        print("Got GFA")

        fout = open(out_fp, "w")

        call(gfa_instance, alignment_fp, fout)

        sys.exit(0)


    out1 = "./simulation/methylation.tmp.tsv"
    out2 = "./simulation/methylation.tsv"
    if True:
        # gfa_instance = gfa.GraphicalFragmentAssembly()
        gfa_instance = gfa.GraphicalFragmentAssemblyMemory()
        gfa_instance.parse("../test_data/chr20-mcpg.gfa")

        fout = open(out1, "w")

        call(gfa_instance, "./simulation/merged.gaf", fout)

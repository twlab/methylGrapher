

import os
import re
import sys
import json
import time
import utility

import multiprocessing
import multiprocessing.queues

import gfa













ref_converted_genotype_dict = {
    "C": "T",
    "G": "A",
}


# Regex Patterns
alignment_csz_tag_start_pattern = re.compile(r"^(\+|-)([ACGTN]*)")
alignment_csz_tag_end_pattern = re.compile(r"(\+|-)([ACGTN]*)$")

alignment_csz_tag_match = re.compile(r":\d*")
alignment_csz_tag_mismatch = re.compile(r"\*\w\w")

alignment_csz_tag_deletion = re.compile(r"-\w")
alignment_csz_tag_insertion = re.compile(r"\+\w")






# Basic Helper Functions
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




def get_best_alignment_from_same_read_pair_OLD(
        lines_for_same_read,
        minimum_identity=0,
        minimum_mapq=0,
        discard_multimapped=True):

    best_alignments = {}
    for a in lines_for_same_read:
        read_index = a[0][1]
        if read_index not in best_alignments:
            best_alignments[read_index] = {}
        read_seq = None
        conversion_type = None
        alignment_score = None

        for tag in a[12:]:
            if tag.startswith("bq:Z:"):
                # Getting the original
                read_seq = tag[5:]

            if tag.startswith("ct:Z:"):
                # Getting the original
                conversion_type = tag[5:]

            if tag.startswith("AS:i:"):
                # Getting the original
                alignment_score = int(tag[5:])

        assert read_seq is not None
        assert conversion_type is not None

        if alignment_score not in best_alignments[read_index]:
            best_alignments[read_index][alignment_score] = []
        best_alignments[read_index][alignment_score].append(a)



    multimapped_count = 0
    low_confidence_count = 0

    best_alignments2 = []
    for read_index in best_alignments.keys():
        highest_score = max(best_alignments[read_index].keys())
        best_alignments[read_index] = best_alignments[read_index][highest_score]

        multimapped = False
        best_alignment = best_alignments[read_index][0]
        best_alignment[5] = alignment_path_parse(best_alignment[5])

        aligned_segments = set(best_alignment[5][0])
        if len(best_alignments[read_index]) > 1:
            for a in best_alignments[read_index][1:]:
                p = set(alignment_path_parse(a[5])[0])
                aligned_segments = aligned_segments.intersection(p)
                # print(a)

            if len(aligned_segments) == 0:
                multimapped = True
                multimapped_count += 1

        if multimapped and discard_multimapped:
            continue

        for i in [2, 3, 6, 7, 8, 9, 10, 11]:
            best_alignment[i] = int(best_alignment[i])

        low_confidence = False
        if best_alignment[11] < minimum_mapq:
            low_confidence = True
        if best_alignment[9] < minimum_identity:
            low_confidence = True
        if low_confidence:
            low_confidence_count += 1
            continue


        best_alignments2.append(best_alignment)


    return best_alignments2, len(best_alignments), multimapped_count, low_confidence_count




def get_best_alignment_from_same_read_pair(
        lines_for_same_read, node_replacement_dict,
        minimum_identity=0,
        minimum_mapq=0,
        discard_multimapped=True):

    # print(minimum_identity, minimum_mapq, discard_multimapped)
    # print(node_replacement_dict)

    aligned_count = 0
    multimapped_count = 0
    low_confidence_count = 0

    best_alignments = []
    for alignment in lines_for_same_read:
        aligned_count += 1

        alignment_score = None
        multimapped = False

        for i in [1, 2, 3, 6, 7, 8, 9, 10, 11]:
            alignment[i] = int(alignment[i])
        alignment[5] = alignment_path_parse(alignment[5])
        segment_list = alignment[5][0]
        for si in range(len(segment_list)):
            segmentID = segment_list[si]
            # print(segmentID, node_replacement_dict)
            if segmentID in node_replacement_dict:
                segment_list[si] = node_replacement_dict[segmentID]
                # print("Replace: ", segmentID, node_replacement_dict[segmentID])
        

        for tag in alignment[12:]:


            if tag.startswith("AS:i:"):
                # Getting the original
                alignment_score = int(tag[5:])

            if tag.startswith("mp:i:"):
                # Getting the original
                mpi = int(tag[5:])
                if mpi == 1:
                    multimapped = True

        assert alignment_score is not None

        if multimapped:
            if discard_multimapped:
                multimapped_count += 1
                continue

        low_confidence = False
        if alignment[11] < minimum_mapq:
            low_confidence = True
        if alignment[9] < minimum_identity:
            low_confidence = True
        if low_confidence:
            low_confidence_count += 1
            continue

        best_alignments.append(alignment)

    return best_alignments, aligned_count, multimapped_count, low_confidence_count



def alignment_to_methylation(best_alignments, sequence_dict, cg_only=True, perform_gcall=False, phred_score_threshold=20):
    error_count = 0
    mcall_per_frag = set()
    gcall_per_frag = set()
    for best_alignment in best_alignments:

        l = best_alignment

        # VG giraffe always outputs + here, just in case
        assert l[4] == "+"

        query_start = l[2]
        query_end = l[3]

        # VG giraffe does not provide stable path_end
        path = l[5]
        path_start = l[7]

        # Figure out the alignment start and end, but to do that, I need the alignment CS tag
        alignment_tag = None
        original_bs_read = None
        read_conversion_type = None
        phred_score = None
        for tag in l[12:]:
            tag = tag.strip()

            if tag.startswith("cs:Z:"):
                alignment_tag = tag[5:]
            if tag.startswith("os:Z:"):
                original_bs_read = tag[5:]
            if tag.startswith("rc:Z:"):
                read_conversion_type = tag[5]
            if tag.startswith("bq:Z:"):
                phred_score = tag[5:]

        if alignment_tag is None or original_bs_read is None or read_conversion_type is None:
            error_count += 1
            continue

        qso, qeo, rl, indels = cs_tag_parse(alignment_tag)
        query_start += qso
        query_end -= qeo
        path_end = path_start + rl
        l[8] = path_end


        path_sequences = []
        for pi in range(len(path[0])):
            segmentID = path[0][pi]
            seq = sequence_dict[segmentID]
            if path[1][pi] < 0:
                seq = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            path_sequences.append(seq)
        path_sequence = "".join(path_sequences)


        path_seq_portion = path_sequence[path_start:path_end]
        bs_read_portion = original_bs_read[query_start:query_end]
        phred_score_portion = phred_score[query_start:query_end]


        # Reconstruct the alignment
        for indel in indels:
            indel_start = indel[0]
            indel_leng = indel[1]
            if indel[2] == "-":
                bs_read_portion = bs_read_portion[:indel_start] + " " * indel_leng + bs_read_portion[indel_start:]
                phred_score_portion = phred_score_portion[:indel_start] + " " * indel_leng + phred_score_portion[indel_start:]
            else:
                bs_read_portion = bs_read_portion[:indel_start] + bs_read_portion[indel_start + indel_leng:]
                phred_score_portion = phred_score_portion[:indel_start] + phred_score_portion[indel_start + indel_leng:]
        

        

        if len(path_seq_portion) != rl or len(path_seq_portion) != len(bs_read_portion) or len(bs_read_portion) != len(phred_score_portion):
            error_count += 1

            continue
            print(
                "\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n\033[91mERROR\033[0m\n")
            print("Error rate:", (counter3 / counter1) * 100, "%")
            time.sleep(0.3)
            continue


        if False:
            # Just for debugging
            identity = 0

            # print(len(path_seq_portion), path, path_start, path_end, rl)
            # print(path_sequences)
            for i in range(len(path_seq_portion)):
                if path_seq_portion[i] == bs_read_portion[i]:
                    identity += 1
            identity = "\033[92m%0.0f%%\033[0m" % (identity / len(path_seq_portion) * 100)

            print("Identity:", identity, path)
            # print(path_start, path_end)
            # print(path_sequence)
            # print(query_start, query_end)
            # print(original_bs_read)
            # print()
            # path_sequence[path_start:path_end]
            endx = min(path_end + 20, len(path_sequence))
            alignment_tag_highlight = ""
            for i in range(len(alignment_tag)):
                if alignment_tag[i] == "+":
                    alignment_tag_highlight += "\033[92m+\033[0m"
                elif alignment_tag[i] == "-":
                    alignment_tag_highlight += "\033[91m-\033[0m"
                else:
                    alignment_tag_highlight += " "
            print(" " * (len(path_seq_portion) - 5) + path_sequence[path_end - 5:endx])
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
            # print(original_bs_read)
            # print(alignment_tag)
            # print(alignment_tag_highlight)
            print()


        # l = len(path_seq_portion)
        pl = len(path_sequence)

        segment_index = 0
        segmentID = path[0][segment_index]
        segment_orientation = path[1][segment_index]
        segment_length = len(path_sequences[segment_index])
        path_len_so_far = segment_length
        # last_segmentID = None
        # print(path)
        # print(path_sequences)
        # print(path_original_sequence)

        for i, ref_base in enumerate(path_seq_portion):
            path_pos = i + path_start
            # ref_base = path_seq_portion[i]

            read_base = bs_read_portion[i]
            phred_score_letter_base = phred_score_portion[i]
            if read_base == " ":
                assert phred_score_letter_base == " "
                continue

            phred_score_int_base = utility.phred_to_int(phred_score_letter_base)


            category = "U"

            if phred_score_int_base < phred_score_threshold:
                # For both genotype and methylation calling
                continue

            if ref_base not in "CG":
                continue

            interesting = [False, False]
            # Genotyping call
            if perform_gcall:
                # The reason for genotype call is to remove the case of C->T and G->A mutation. This will cause inaccurate methylation call.
                # Why not care C->G, C->A, G->C, or G->T? Because if it happens, it doesn't interfere with the methylation call.
                if ref_base == "C" and read_conversion_type == "G":
                    interesting[0] = True
                elif ref_base == "G" and read_conversion_type == "C":
                    interesting[0] = True

            # Methylation call
            if ref_base == "C" and read_conversion_type == "C":
                if read_base in "CT":
                    interesting[1] = True
            elif ref_base == "G" and read_conversion_type == "G":
                if read_base in "GA":
                    interesting[1] = True

            if True not in interesting:
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


            if interesting[0]:
                # Genotype call
                a, b = ref_base, read_base
                if segment_orientation < 0:
                    a = utility.atcg_complement_dict.get(ref_base, "N")
                    b = utility.atcg_complement_dict.get(read_base, "N")
                gcall_per_frag.add((segmentID, segment_pos, a, b, phred_score_letter_base))
                continue



            """
            x = path_sequence[path_pos]
            y = path_seq_portion[i]
            z = path_original_sequence[segmentID][segment_pos]
            if segment_orientation < 0:
                z = utl.seq_reverse_complement(z)
            assert x == y and y == z
            """
            # Figure out context
            base_strand = segment_orientation
            if ref_base == "C":
                if path_pos + 1 < pl:
                    ref_base2 = path_sequence[path_pos + 1]
                    if ref_base2 == "G":
                        category = "CG"
                    else:
                        if path_pos + 2 < pl:
                            ref_base3 = path_sequence[path_pos + 2]
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
                    ref_base2 = path_sequence[path_pos - 1]
                    if ref_base2 == "C":
                        category = "CG"
                    else:
                        if path_pos > 1:
                            ref_base3 = path_sequence[path_pos - 2]
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

            """
            if segmentID == "10007" and methylated == 1:
                highlight = ""
                hc = "C"
                if segment_pos > 230:
                    hc = "G"

                for b in bs_read_portion:
                    if b == hc:
                        highlight += "\033[92m*\033[0m"
                    else:
                        highlight += " "
                segment_offset = segment_pos - 100
                if segment_pos > 230:
                    segment_offset = segment_pos - 250

                print(segmentID, segment_pos, base_strand, category, methylated)
                print(ref_base, read_base)
                print("Ref:  ", path_seq_portion)
                print("Read: ", bs_read_portion)
                print("      ", highlight)
                print("      ", " " * segment_offset + "*")
                print()
            """


            # print(i, category, ref_base, read_base, methylated)
            # print(path_seq_portion[i], bs_read_portion[i], )

            if cg_only and category != "CG":
                continue

            d = (segmentID, segment_pos, base_strand, category, methylated)
            # Just to verify fragment coverage
            #if d in mcall_per_frag:
            #    assert best_alignment[0][1] == "R2"
            #    print(d, best_alignments[0][6], best_alignments[1][6])
            #    print(f"{best_alignments[0][7]}-{best_alignments[0][8]} ON {best_alignments[0][5]}")
            #    print(f"{best_alignments[1][7]}-{best_alignments[1][8]} ON {best_alignments[1][5]}")
            mcall_per_frag.add(d)

    mcall_per_frag = list(sorted(mcall_per_frag))
    return mcall_per_frag, gcall_per_frag



def alignment_parse(work_dir, node_replacement_dict, minimum_identity=50, minimum_mapq=20, discard_multimapped=True, pid="0_0"):
    # Alignment count, low confidence count, error count, multi-mapped alignment count
    counter1 = 0
    counter2 = 0
    counter3 = 0
    counter4 = 0

    ts = time.time()


    alignment_fp = work_dir + "/alignment.gaf"


    lines_for_same_read = []
    with open(alignment_fp) as alignment_fh:
        for l in alignment_fh:
            l = l.strip().split("\t")

            query_name = l[0]
            if len(lines_for_same_read) == 0:
                lines_for_same_read.append(l)
                continue

            if query_name == lines_for_same_read[0][0]:
                lines_for_same_read.append(l)
                continue

            # Figure out the best alignment for the read pair.
            best_alignments, aligned_count, mpc, lcc = get_best_alignment_from_same_read_pair(
                lines_for_same_read, node_replacement_dict,
                minimum_identity=minimum_identity,
                minimum_mapq=minimum_mapq,
                discard_multimapped=discard_multimapped
            )

            # Reset the buffer
            lines_for_same_read = [l]

            counter1 += aligned_count
            counter2 += lcc
            counter4 += mpc

            yield best_alignments


    best_alignments, aligned_count, mpc, lcc = get_best_alignment_from_same_read_pair(
        lines_for_same_read, node_replacement_dict,
        minimum_identity=minimum_identity,
        minimum_mapq=minimum_mapq,
        discard_multimapped=discard_multimapped
    )
    counter1 += aligned_count
    counter2 += lcc
    counter4 += mpc
    yield best_alignments

    # print(counter1, counter2, counter3, counter4)
    total_alignment_count, low_confidence_count, error_count, multimapped_count = counter1, counter2, counter3, counter4

    logging = f"""
Total aligned reads: {total_alignment_count}
Total multi-mapped reads: {multimapped_count}
Low confidence count (low quality alignment in terms of block length and MapQ): {low_confidence_count}
"""
    # logging = f"{total_alignment_count}\t{multimapped_count}\t{low_confidence_count}"
    #  Inconsistent alignment entries: {error_count}

    # print(logging)

    log_fh = open(f"{work_dir}/report.txt", "a")
    log_fh.write(logging)
    log_fh.close()








# Deprecated
def get_methylation_output_tmp_fps_from_work_dir(work_dir):
    raise RuntimeError
    res = []
    for i in range(0, 10):
        fn = f"{work_dir}/mcall.{i}.tmp"
        res.append(fn)
    return res


# Single thread version
def call_single(
        gfa_fp, work_dir,
        cg_only=True,
        minimum_identity=50, minimum_mapq=20, discard_multimapped=True
):

    gfa_instance = gfa.GraphicalFragmentAssemblyMemory()
    gfa_instance.parse(gfa_fp)

    methylation_output = []
    for i in range(10):
        mo = f"{work_dir}/mcall.0.{i}.tmp"
        methylation_output.append(open(mo, "w"))


    for best_alignments in alignment_parse(
            work_dir,
            minimum_identity=minimum_identity, minimum_mapq=minimum_mapq, discard_multimapped=discard_multimapped
    ):

            segments = set()
            for best_alignment in best_alignments:
                segments = segments.union(best_alignment[5][0])
            sequence_dict = gfa_instance.get_sequences_by_segment_ID(segments)


            ts = time.time()
            mcall_per_frag, gcall_per_frag = alignment_to_methylation(best_alignments, sequence_dict, cg_only=cg_only)
            #print(f"Call: {time.time() - ts:.8f}")
            #print()


            # Merge the methylation calls for the same read pair/fragment.
            # So the coverage is by fragment rather than by read.
            mcall_per_frag = list(sorted(mcall_per_frag))
            last_segmentID = None
            for (segmentID, segment_pos, base_strand, category, methylated) in mcall_per_frag:
                mcall_tmp = methylation_output[int(segmentID[-1])]

                if segmentID != last_segmentID:
                    mcall_tmp.write(f"{segmentID}\t{segment_pos}\t{base_strand}\t{category}\t{methylated}\n")
                    last_segmentID = segmentID
                else:
                    mcall_tmp.write(f"\t{segment_pos}\t{base_strand}\t{category}\t{methylated}\n")



    return None








class WorkerConterAndTimer(object):

    def __init__(self, process_name="Worker", interval=10000):
        self._counter = 0
        self._interval = interval
        self._process_name = process_name

    def start(self):
        self._ts = time.time()
        self._ts_last = self._ts

    def count(self):
        self._counter += 1
        if self._counter % self._interval == 0:
            speed_long = self._counter / (time.time() - self._ts)
            speed_short = self._interval / (time.time() - self._ts_last)
            self._ts_last = time.time()


            print(f"{self._process_name}: {speed_short:.1f}/sec {speed_long:.1f}/sec (ave)", file=sys.stderr)




# Worker Deamon Processes Functions
def alignment_parse_worker(pid, work_dir, node_replacement_dict_fp, output_queue, minimum_identity=50, minimum_mapq=20, discard_multimapped=True, genotyping_cytosine=False, batch_size=4096):
    wct = WorkerConterAndTimer(process_name=f"Alignment Parsing {pid}")
    wct.start()


    node_replacement_dict = {}
    if genotyping_cytosine:
        with open(node_replacement_dict_fp) as node_replacement_dict_fh:
            node_replacement_dict_raw = json.load(node_replacement_dict_fh)

        for v in node_replacement_dict_raw.values():
            node_replacement_dict.update(v)

    batch = []
    for r in alignment_parse(
        work_dir, node_replacement_dict,
        minimum_identity=minimum_identity, minimum_mapq=minimum_mapq, discard_multimapped=discard_multimapped,
        pid=f"{pid}"
    ):
        wct.count()

        if len(r) == 0:
            # No need to pass empty list
            continue

        batch.append(r)
        if len(batch) >= batch_size:
            output_queue.put(batch)
            batch = []

    if len(batch) > 0:
        output_queue.put(batch)
        batch = []


    print(f"Alignment Parser finished", file=sys.stderr)
    return


def gfa_worker(pid, gfa_fp, input_queue, output_queue, batch_size=4096):
    wct = WorkerConterAndTimer(process_name=f"Segment Sequence Fetching {pid}")



    print(f"GFA Worker started", file=sys.stderr)
    gfa_instance = gfa.GraphicalFragmentAssemblyMemory()
    gfa_instance.parse(gfa_fp)
    print(f"GFA Worker finished parsing GFA", file=sys.stderr)

    wct.start()

    batch = []
    while True:
        best_alignment_batch = input_queue.get()
        if best_alignment_batch is None:
            break

        for best_alignments in best_alignment_batch:
            segments = set(best_alignments[0][5][0])
            for best_alignment in best_alignments[1:]:
                segments = segments.union(set(best_alignment[5][0]))
            sequence_dict = gfa_instance.get_sequences_by_segment_ID(segments)

            wct.count()

            batch.append((best_alignments, sequence_dict))
            if len(batch) >= batch_size:
                output_queue.put(batch)
                batch = []


    if len(batch) > 0:
        output_queue.put(batch)

    print(f"GFA Worker finished", file=sys.stderr)
    return




def alignment_to_methylation_worker(pid, work_dir, cg_only, input_queue, genotyping_cytosine=False):

    wct = WorkerConterAndTimer(process_name=f"MCall_{pid}", interval=1000)

    methylation_output = []
    for i in range(100):
        mo = f"{work_dir}/mcall.{pid}.{i}.tmp"
        methylation_output.append(open(mo, "w"))

    genotype_output = []
    if genotyping_cytosine:
        for i in range(100):
            go = f"{work_dir}/gcall.{pid}.{i}.tmp"
            genotype_output.append(open(go, "w"))

    result = []
    result_gcall = []

    init = False
    while True:
        # print(f"MCall Worker {pid} waiting", file=sys.stderr)
        pair_batch = input_queue.get()
        # print(f"MCall Worker {pid} got batch, {pair_batch == None}", file=sys.stderr)
        if pair_batch is None:
            break
        if not init:
            wct.start()
            init = True

        for best_alignments, sequence_dict in pair_batch:
            wct.count()

            mcall_per_frag, gcall_per_frag = alignment_to_methylation(best_alignments, sequence_dict, cg_only=cg_only, perform_gcall=genotyping_cytosine)
            result += list(mcall_per_frag)
            result_gcall += list(gcall_per_frag)

        last_segmentID = None
        for (segmentID, segment_pos, base_strand, category, methylated) in sorted(result):
            # print((segmentID, segment_pos, base_strand, category, methylated))
            mcall_tmp = methylation_output[int(segmentID[-2:])]

            if segmentID != last_segmentID:
                mcall_tmp.write(f"{segmentID}\t{segment_pos}\t{base_strand}\t{category}\t{methylated}\n")
                last_segmentID = segmentID
            else:
                mcall_tmp.write(f"\t{segment_pos}\t{base_strand}\t{category}\t{methylated}\n")


        if not genotyping_cytosine:
            continue
        last_segmentID = None
        for gcall0 in result_gcall:
            segmentID = gcall0[0]
            gcall_fh = genotype_output[int(segmentID[-2:])]

            if segmentID != last_segmentID:
                gcall_fh.write("\t".join(list(map(str, gcall0[0:]))) + "\n")
                last_segmentID = segmentID
            else:
                gcall_fh.write("\t"+"\t".join(list(map(str, gcall0[1:]))) + "\n")

        result = []
        result_gcall = []

    print(f"MCall Worker {pid} finished", file=sys.stderr)
    return


# Deprecated. The inter-process communication is too expensive
def writer_worker(input_queue, methylation_output_fp):

    wct = WorkerConterAndTimer(process_name="Writer", interval=100)
    wct.start()

    methylation_output = []
    for mo in methylation_output_fp:
        methylation_output.append(open(mo, "w"))

    while True:

        mcall_per_frags = input_queue.get()
        if mcall_per_frags is None:
            break
        wct.count()
        # print(mcall_per_frag)
        # print(f"Writer got batch, {len(mcall_per_frags)}")

        last_segmentID = None
        for (segmentID, segment_pos, base_strand, category, methylated) in mcall_per_frags:
            # print((segmentID, segment_pos, base_strand, category, methylated))
            mcall_tmp = methylation_output[int(segmentID[-1])]

            if segmentID != last_segmentID:
                mcall_tmp.write(f"{segmentID}\t{segment_pos}\t{base_strand}\t{category}\t{methylated}\n")
                last_segmentID = segmentID
            else:
                mcall_tmp.write(f"\t{segment_pos}\t{base_strand}\t{category}\t{methylated}\n")

    print(f"Writer finished", file=sys.stderr)





# Debugging Process, just to clear the queue
def debug_worker(input_queue, output_queue):

    counter = 0
    start_ts = time.time()

    while True:

        mcall_per_frag = input_queue.get()
        if mcall_per_frag is None:
            break



# Multi process version
def call_parallel(
        # Essential
        gfa_fp, node_replacement_dict_fp, work_dir,
        cg_only=True, genotyping_cytosine=False,

        # Alignment Filtering Parameters
        minimum_identity=50, minimum_mapq=20, discard_multimapped=True,

        # Performance Parameters
        process_count=3, batch_size=4096, gfa_worker_num=1, alignment_parse_worker_num=1
):

    # maxsize is needed to prevent too much memory usage
    maxsize = 100
    # q0 = multiprocessing.Queue(maxsize=maxsize) // Not needed any more
    q1 = multiprocessing.Queue(maxsize=maxsize)
    q2 = multiprocessing.Queue(maxsize=maxsize)


    mcall_worker_num = int(process_count)
    if mcall_worker_num < 1:
        mcall_worker_num = 1

    for i in [process_count, gfa_worker_num, mcall_worker_num]:
        assert i > 0
        assert isinstance(i, int)


    alignment_parse_process = multiprocessing.Process(
        name=f"MGAR0",
        target=alignment_parse_worker,
        args=(i, work_dir, node_replacement_dict_fp, q1, minimum_identity, minimum_mapq, discard_multimapped, genotyping_cytosine),
        kwargs={"batch_size": batch_size}
    )
    alignment_parse_process.start()


    gfa_process_pool = []
    for i in range(gfa_worker_num):
        gfa_process = multiprocessing.Process(
            name=f"MGGFA{i}",
            target=gfa_worker,
            args=(i, gfa_fp, q1, q2),
            kwargs={"batch_size": batch_size}
        )
        gfa_process.start()
        gfa_process_pool.append(gfa_process)

    mcall_process_pool = []
    for i in range(mcall_worker_num):
        mp = multiprocessing.Process(
            name=f"MGMCall{i}",
            target=alignment_to_methylation_worker,
            args=(i, work_dir, cg_only, q2, genotyping_cytosine)
        )
        mp.start()
        mcall_process_pool.append(mp)



    # writter_process = multiprocessing.Process(name=f"MGWriter", target=writer_worker, args=(q3, methylation_output_fp))
    # writter_process.start()

    #debug_process = multiprocessing.Process(name=f"DEBUG", target=debug_worker, args=(q3, q5))
    #debug_process.start()

    #processes = [alignment_parse_process, writter_process] + gfa_process_pool + mcall_process_pool
    #for p in processes:
    #    p.join()

    # I am not whether it is a good idea to manually manage the processes... But anyway.
    finished = False

    send_kill_signal_to_gfa_workers = False
    send_kill_signal_to_mcall_workers = False

    while not finished:
        time.sleep(1)

        """
        processes = [alignment_parse_process, writter_process] + gfa_process_pool + mcall_process_pool
        for p in processes:
            alive = ""
            if p.is_alive():
                alive = "alive"
            print(p.name, alive, file=sys.stderr)

        print(file=sys.stderr)

        for qi, q in enumerate([q1, q2, q3]):
            print(f"Queue {qi} Empty: {q.empty()}", file=sys.stderr)

            while not q.empty():
                print(q.get(), file=sys.stderr)

        print(file=sys.stderr)
        """



        alignment_parse_process.join()

        if not send_kill_signal_to_gfa_workers:
            for p in gfa_process_pool:
                q1.put(None)
            send_kill_signal_to_gfa_workers = True
        for p in gfa_process_pool:
            p.join()

        if not send_kill_signal_to_mcall_workers:
            for p in mcall_process_pool:
                q2.put(None)
            send_kill_signal_to_mcall_workers = True
        for p in mcall_process_pool:
            p.join()

        # if not send_kill_signal_to_writer_workers:
        #    q3.put(None)
        #    send_kill_signal_to_writer_workers = True
        # writter_process.join()

        finished = True


    return None







# Second step
def mcall_tmp_to_final(work_dir, final_out, max_pid=1000):
    # Deprecated: This is the old version

    segmentID = ""
    for i in range(100):

        methylation = {}
        for pid in range(max_pid):
            fn = f"{work_dir}/mcall.{pid}.{i}.tmp"
            if not os.path.exists(fn):
                continue
            # print(fn)

            fh = open(fn)
            for l in fh:
                l = l.strip().split("\t")
                if len(l) == 5:
                    segmentID = l.pop(0)
                key = (segmentID, l[0], l[1], l[2])
                if key not in methylation:
                    methylation[key] = [0, 0]
                methylation[key][int(l[3])] += 1

            fh.close()

            os.remove(fn)

        for key in methylation.keys():
            unmethylated = methylation[key][0]
            methylated = methylation[key][1]
            coverage = unmethylated + methylated
            mlevel = methylated / coverage
            final_out.write(
                f"{key[0]}\t{key[1]}\t{key[2]}\t{key[3]}\t{unmethylated}\t{methylated}\t{coverage}\t{mlevel}\n"
            )

    return None






def heterozygous_prob(n, k):
    if k > n - k:  # Take advantage of symmetry
        k = n - k
    result = 1
    for i in range(1, k + 1):
        result = result * (n - i + 1) // i

    for d in range(n):
        if result % 2 == 0 or result > 1e40:
            result = result // 2
        else:
            result = result / 2

    return result


def merge_genotype_func(work_dir, split_i, max_pid=100, coverage_threshold=5):
    homozygous_ref_prior_prob = 1 / 3
    homozygous_alt_prior_prob = 1 / 3
    heterozygous_ref_alt_prior_prob = 1 / 3

    # This step piles up the data
    pileup = {}
    res = {}
    for pid_i in range(max_pid):
        f2 = f"{work_dir}/gcall.{pid_i}.{split_i}.tmp"

        if not os.path.exists(f2):
            continue
        last_segment = None
        with open(f2) as f:
            for l in f:
                l = l.strip().split()
                if len(l) == 5:
                    last_segment = l.pop(0)

                if last_segment not in pileup:
                    pileup[last_segment] = {}

                pos = int(l[0])
                ref = l[1]
                alt = l[2]
                qual = l[3]

                k = (pos, ref)
                if k not in pileup[last_segment]:
                    pileup[last_segment][k] = {}

                if alt not in pileup[last_segment][k]:
                    pileup[last_segment][k][alt] = []

                pileup[last_segment][k][alt].append(qual)

        os.remove(f2)

    # Genotyping call
    i = 0
    genotype_interesting_case_pileup = []
    for segment in pileup.keys():
        for (pos, ref), alt_phred_scores in pileup[segment].items():
            i += 1

            genotypes_sorted = sorted(alt_phred_scores.keys(), key=lambda x: len(alt_phred_scores[x]), reverse=True)
            top2_genotypes = genotypes_sorted[:2]
            if not ref_converted_genotype_dict[ref] in top2_genotypes:
                continue

            combined_coverage = 0
            for genotype in top2_genotypes:
                combined_coverage += len(alt_phred_scores[genotype])

            if combined_coverage < coverage_threshold:
                continue

            ref_count = len(alt_phred_scores.get(ref, []))
            ref_coverted_count = len(alt_phred_scores.get(ref_converted_genotype_dict[ref], []))

            # check whether the ref_converted genotype is most likely homozygous?
            # This is a short cut, if all reads support ref_converted, then no need to run bayes model
            homozygous_ref_converted_shortcut = False
            if ref_coverted_count == combined_coverage:
                # Homozygous converted genotype, false positive
                homozygous_ref_converted_shortcut = True

                # ref, alt, cov, FP, likelihood of homozygous ref, likelihood of heterozygous, likelihood of homozygous alt
                res[(segment, pos)] = [ref, top2_genotypes[0], combined_coverage, True, "*", "*", "*"]
                continue

            if ref not in top2_genotypes:
                res[(segment, pos)] = [ref, "|".join(top2_genotypes), combined_coverage, True, "*", "*", "*"]
                continue

            assert ref in top2_genotypes
            assert ref_converted_genotype_dict[ref] in top2_genotypes

            likelihoods = [homozygous_ref_prior_prob, homozygous_alt_prior_prob, heterozygous_ref_alt_prior_prob]
            for qual in alt_phred_scores[ref]:
                likelihoods[1] *= utility.phred_to_prob(qual)
                likelihoods[0] *= (1 - utility.phred_to_prob(qual))

            for qual in alt_phred_scores[ref_converted_genotype_dict[ref]]:
                likelihoods[1] *= (1 - utility.phred_to_prob(qual))
                likelihoods[0] *= utility.phred_to_prob(qual)

            ref_and_converted_cov = ref_count + ref_coverted_count
            p_het = heterozygous_prob(ref_and_converted_cov, ref_count)
            likelihoods[2] *= p_het

            ls = sum(likelihoods)
            if ls == 0:
                ls = 1
            likelihoods = [l / ls for l in likelihoods]

            if likelihoods[1] > likelihoods[0] and likelihoods[1] > likelihoods[2]:
                res[(segment, pos)] = [ref, "|".join(top2_genotypes), combined_coverage, True] + likelihoods
                continue

            # print(segment, pos, ref, alt_phred_scores)
            # sys.exit(3)

            # l = [segment, pos, "+", "GT", ref_count, ref_coverted_count, ref_count+ref_coverted_count, ref_coverted_count/(ref_count+ref_coverted_count)]
            # genotype_interesting_case_pileup.append("\t".join(map(str, l)) + "\n")
            # out_fh.write("\t".join(map(str, l)) + "\n")
            # print(l)

    return res


def merge_methylation_func(work_dir, split_i, max_pid=100):
    methylation = {}
    segmentID = ""
    for pid in range(max_pid):
        # working with methylation now
        fn = f"{work_dir}/mcall.{pid}.{split_i}.tmp"
        if not os.path.exists(fn):
            continue

        with open(fn) as fh:
            for l in fh:
                l = l.strip().split("\t")
                if len(l) == 5:
                    segmentID = l.pop(0)
                key = (segmentID, l[0], l[1], l[2])
                if key not in methylation:
                    methylation[key] = [0, 0]
                methylation[key][int(l[3])] += 1

        os.remove(fn)

    return methylation


def extraction_merge_and_cleanup_by_split(work_dir, split_i, max_pid=100, coverage_threshold=5):
    genotypes = merge_genotype_func(work_dir, split_i, max_pid=max_pid, coverage_threshold=coverage_threshold)
    methylation = merge_methylation_func(work_dir, split_i, max_pid=max_pid)

    mres = []
    gres = []

    # print(len(genotypes), len(methylation))
    # print(list(genotypes.keys())[0])
    # print(list(methylation.keys())[0])

    for key in methylation.keys():
        k2 = (key[0], int(key[1]))
        if k2 in genotypes:
            # print("Genotype and methylation overlap", key, k2)
            continue

        unmethylated = methylation[key][0]
        methylated = methylation[key][1]
        coverage = unmethylated + methylated
        mlevel = methylated / coverage
        mres.append(f"{key[0]}\t{key[1]}\t{key[2]}\t{key[3]}\t{unmethylated}\t{methylated}\t{coverage}\t{mlevel}\n")

    for key in genotypes.keys():
        gv = genotypes[key]
        line = list(key) + gv
        gres.append("\t".join(map(str, line)) + "\n")

    return mres, gres


def extraction_merge_and_cleanup_by_split_single_arg(arg):
    # print(arg)
    assert len(arg) == 4
    work_dir, split_i, max_pid, coverage_threshold = arg
    return extraction_merge_and_cleanup_by_split(work_dir, split_i, max_pid=max_pid,
                                                 coverage_threshold=coverage_threshold)


# Single threaded version
def extraction_merge_and_cleanup(work_dir, max_pid=100, coverage_threshold=5):
    methylation_out_fp = os.path.join(work_dir, "graph.methyl")
    genotype_out_fp = os.path.join(work_dir, "genotype.info.txt")

    genotyped = False

    mfh = open(methylation_out_fp, "w")
    gfh = open(genotype_out_fp, "w")

    for split_i in range(100):
        # print(split_i)
        mres, gres = extraction_merge_and_cleanup_by_split(work_dir, split_i, max_pid=max_pid,
                                                           coverage_threshold=coverage_threshold)

        for mline in mres:
            mfh.write(mline)

        for gline in gres:
            genotyped = True
            gfh.write(gline)

    mfh.close()
    gfh.close()

    if not genotyped:
        os.remove(genotype_out_fp)

    return


def cytosine_CG_validation(index_prefix, wd):
    cpg_table_fp = f"{index_prefix}.cpg.tsv"

    genotype_info_fp = f"{wd}/genotype.info.txt"
    input_cytosine_table_fp = f"{wd}/graph.methyl"

    assert os.path.exists(cpg_table_fp)
    assert os.path.exists(input_cytosine_table_fp)
    if not os.path.exists(genotype_info_fp):
        return

    invalid_cytosines = set()
    with open(genotype_info_fp) as fh:
        for l in fh:
            segID, pos, *others = l.strip().split("\t")
            pos = int(pos)

            invalid_cytosines.add((segID, pos))

    with open(cpg_table_fp) as fh:
        for l in fh:
            cpgid, segid1, pos1, segid2, pos2, *others = l.strip().split("\t")
            c1 = (segid1, int(pos1))
            c2 = (segid2, int(pos2))

            invalid = c1 in invalid_cytosines or c2 in invalid_cytosines
            if invalid:
                invalid_cytosines.add(c1)
                invalid_cytosines.add(c2)

    f1p = f"{wd}/graph.methyl.val"
    f1h = open(f1p, "w")
    with open(input_cytosine_table_fp) as f2h:
        for l in f2h:
            segID, pos, *others = l.strip().split("\t")
            pos = int(pos)

            if (segID, pos) not in invalid_cytosines:
                f1h.write(l)
    f1h.close()

    # mv f1p f2h
    os.remove(input_cytosine_table_fp)
    os.rename(f1p, input_cytosine_table_fp)

    return


# multithreaded version
def extraction_merge_and_cleanup_mp(work_dir, cpu_count, max_pid=100, coverage_threshold=5):
    pool = multiprocessing.Pool(processes=cpu_count)
    args = [(work_dir, split_i, max_pid, coverage_threshold) for split_i in range(100)]
    results = pool.imap_unordered(
        extraction_merge_and_cleanup_by_split_single_arg,
        args
    )
    pool.close()

    methylation_out_fp = os.path.join(work_dir, "graph.methyl")
    genotype_out_fp = os.path.join(work_dir, "genotype.info.txt")
    mfh = open(methylation_out_fp, "w")
    gfh = open(genotype_out_fp, "w")
    genotyped = False

    for mres, gres in results:
        for mline in mres:
            mfh.write(mline)

        for gline in gres:
            genotyped = True
            gfh.write(gline)

    if not genotyped:
        os.remove(genotype_out_fp)

    pool.join()

    return





def mcall_main(
        work_dir, index_prefix,
        cg_only=True, genotyping_cytosine=False,
        minimum_identity=50, minimum_mapq=20, discard_multimapped=True,
        process_count=1, alignment_parse_worker_num=1, gfa_worker_num=1, batch_size=4096

):
    #print("DEV version")
    #sys.stderr.write("DEV version")

    gfa_fp = f"{index_prefix}.wl.gfa"
    node_replacement_dict_fp = f"{index_prefix}.wl.node.replacement.json"

    start_ts = time.time()
    if process_count <= -1:
        # Not maintained
        raise RuntimeError("Single process mode is not supported any more. Please use multi-process mode.")
        call_single(
            gfa_fp, work_dir,
            cg_only=cg_only,
            minimum_identity=minimum_identity,
            minimum_mapq=minimum_mapq,
            discard_multimapped=discard_multimapped
        )

    else:
        call_parallel(
            # Essential
            gfa_fp, node_replacement_dict_fp, work_dir,
            cg_only=cg_only, genotyping_cytosine=genotyping_cytosine,

            # Alignment Filtering Parameters
            minimum_identity=minimum_identity,
            minimum_mapq=minimum_mapq,
            discard_multimapped=discard_multimapped,

            # Performance Parameters
            process_count=process_count, alignment_parse_worker_num=1, gfa_worker_num=gfa_worker_num,
            batch_size=batch_size
        )

    step1_time = time.time() - start_ts

    # final_out = open(f"{work_dir}/graph.methyl", "w")
    # mcall_tmp_to_final(work_dir, final_out)
    # final_out.close()

    cpu_count = process_count
    if process_count > 10:
        # Each process takes roughly 5gb of memory, using too many processes will cause high memory usage
        cpu_count = 10

    extraction_merge_and_cleanup_mp(work_dir, cpu_count, max_pid=100, coverage_threshold=5)

    if cg_only and genotyping_cytosine:
        cytosine_CG_validation(index_prefix, work_dir)

    step2_time = time.time() - start_ts - step1_time

    return step1_time, step2_time





if __name__ == "__main__":

    mcall_main













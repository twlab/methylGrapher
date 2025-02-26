import os
import sys
import gzip
import resource
import multiprocessing
import time

import mcall
import utility




tmp_alignment_file_count = 1000


# The OS may have limit of how many file can you open at the same time.
soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
assert tmp_alignment_file_count+500 < hard_limit
resource.setrlimit(resource.RLIMIT_NOFILE, (tmp_alignment_file_count+500, hard_limit))












conversion_types = ["C2T", "G2A"]


def alignment_clenup(work_dir):
    # Converted FASTQ files removal
    for ct in conversion_types:
        for ri in ["R1", "R2"]:
            for fmt in ["fastq", "fastq.gz"]:
                fn = f"{work_dir}/{ct}.{ri}.{fmt}"
                if os.path.exists(fn):
                    os.remove(fn)

    # TMP alignment files removal
    for i in range(tmp_alignment_file_count):
        alignment_out = f"{work_dir}/alignment.{i}.gaf"
        if os.path.exists(alignment_out):
            os.remove(alignment_out)

    return


def alignment(work_dir="./", index_prefix="", output_format="gaf", thread=1, directional=True, compress=True, vg_path="vg"):
    index_prefix_ct = index_prefix + ".wl.C2T"
    index_prefix_ga = index_prefix + ".wl.G2A"

    read1_alignment_type = ["C2T"]
    read2_alignment_type = ["G2A"]

    if not directional:
        read1_alignment_type = conversion_types
        read2_alignment_type = conversion_types

    alignment_file_path = []
    alignment_outs = {}
    for i in range(tmp_alignment_file_count):
        alignment_out = f"{work_dir}/alignment.{i}.{output_format}"
        alignment_file_path.append(alignment_out)

        fh = open(alignment_out, "w")
        alignment_outs[i] = fh

    for ref_type in conversion_types:
        for read_type1 in read1_alignment_type:
            for read_type2 in read2_alignment_type:
                if read_type1 == read_type2:
                    continue

                # TODO change log for single end mode
                print(f"Aligning R1({read_type1}) & R2({read_type2}) on reference({ref_type})")

                fq1 = f"{work_dir}/{read_type1}.R1.fastq"
                fq2 = f"{work_dir}/{read_type2}.R2.fastq"

                if compress:
                    fq1 += ".gz"
                    fq2 += ".gz"

                if ref_type == "G2A":
                    index_prefix = index_prefix_ga
                else:
                    index_prefix = index_prefix_ct

                giraffe_input = f"-f {fq1}"
                if os.path.exists(fq2):
                    giraffe_input += f" -f {fq2}"

                alignment_log = f"{work_dir}/alignment.Ref_{ref_type}.R1_{read_type1}.R2_{read_type2}.log"


                dist_fp = f"{index_prefix}.dist"
                gbz_fp = f"{index_prefix}.giraffe.gbz"

                min1_fp = f"{index_prefix}.min"
                min2_fp = f"{index_prefix}.shortread.withzip.min"
                zipcode_fp = f"{index_prefix}.shortread.zipcodes"


                assert os.path.exists(dist_fp)
                assert os.path.exists(gbz_fp)

                assert os.path.exists(min1_fp) or os.path.exists(min2_fp)
                if os.path.exists(min2_fp):
                    assert os.path.exists(zipcode_fp)


                index_params = f"-Z {gbz_fp} -d {dist_fp}"
                if os.path.exists(min1_fp):
                    # v1.62.0 and lower
                    index_params += f" -m {min1_fp}"
                else:
                    # v1.63.0 and higher
                    index_params += f" -m {min2_fp} -z {zipcode_fp}"


                cmd = f"{vg_path} giraffe -p -t {thread} -o {output_format} -M 2 --named-coordinates {index_params} {giraffe_input}"
                # print(cmd)
                with open(alignment_log, "w") as alignment_log_fh:
                    alignment_log_fh.write("Command used: \n")
                    alignment_log_fh.write(cmd + "\n\n")

                se = utility.SystemExecute()
                fout, flog = se.execute(cmd, stdout=None, stderr=alignment_log)
                for line in fout:
                    line = line.decode("utf-8")

                    # Skip unaligned
                    l = line.strip().split("\t")
                    asterisk = False
                    for i in [2, 3, 6, 7, 8, 9, 10, 11]:
                        if l[i] == "*":
                            asterisk = True
                            break
                    if asterisk:
                        continue

                    # Just to split alignment to multiple files, reduce memory usage for sorting
                    ind = int(l[0].split("_")[2])
                    alignment_out_fh = alignment_outs[ind]
                    alignment_out_fh.write(line)

                se.wait()

    for fh in alignment_outs.values():
        fh.close()

    # Not needed any more
    # for fn in alignment_file_path:
    # cmd = f"sort -o {fn} {fn}"

    # se = utility.SystemExecute()
    # fout, flog = se.execute(cmd, stdout=None, stderr=None)
    # se.wait()

    return


def tmp_gaf_processing(tmp_gaf_fp):
    reads = {}
    result_gaf_str = ''
    with open(tmp_gaf_fp) as f:
        for i, l in enumerate(f):
            l = l.strip().split('\t')
            newl = l[:]

            query_name_complex = l[0].split('_')
            query_name = query_name_complex[0]
            read_conversion = query_name_complex[1][0] + query_name_complex[1][2]
            original_seq = query_name_complex[3]

            newl[0] = query_name

            pop_i = -1
            r1 = False
            r2 = False
            for j, e in enumerate(newl):
                if e.startswith("fn:Z:"):
                    r1 = True
                    pop_i = j
                if e.startswith("fp:Z:"):
                    r2 = True
                    pop_i = j

            assert not (r1 and r2)
            newl.pop(pop_i)

            ri = 1 if r1 else 2
            ri_tag = f'ri:i:{ri}'
            newl.append(ri_tag)

            newl.append(f'os:Z:{original_seq}')

            newl.append(f'rc:Z:{read_conversion}')

            if query_name not in reads:
                reads[query_name] = [[], []]

            reads[query_name][ri - 1].append(newl)

    counter = [0, 0, 0, 0, 0]
    for query_name, read_pair_alignments in reads.items():
        # print(f'Processing {query_name}')

        for read_alignments in read_pair_alignments:
            counter[0] += 1

            if len(read_alignments) == 0:
                counter[1] += 1
                continue

            elif len(read_alignments) == 1:
                counter[2] += 1
                read_alignment = read_alignments[0]
                result_gaf_str += '\t'.join(read_alignment) + '\n'
                continue

            else:
                # Determine best alignment
                best_alignments = []
                best_score = -1

                for read_alignment in read_alignments:
                    mapq = int(read_alignment[11])

                    if mapq > best_score:
                        best_score = mapq

                for read_alignment in read_alignments:
                    mapq = int(read_alignment[11])
                    if mapq == best_score:
                        best_alignments.append(read_alignment)

                if len(best_alignments) == 1:
                    counter[3] += 1
                    result_gaf_str += '\t'.join(best_alignments[0]) + '\n'
                    continue

                # Determine multimapping from now on
                gaf_line_set = set()
                for read_alignment in best_alignments:
                    gaf_line = '\t'.join(read_alignment) + '\n'
                    gaf_line_set.add(gaf_line)

                # Are all the alignments the same? If so, just output one, and it is not multimapping
                if len(gaf_line_set) == 1:
                    counter[3] += 1
                    result_gaf_str += '\t'.join(best_alignments[0]) + '\n'
                    continue

                # distinguish by alignment score.
                best_score = -1
                best_alignments_by_as = []
                for read_alignment in best_alignments:
                    for e in read_alignment:
                        if e.startswith('AS:i:'):
                            score = int(e.split(':')[-1])
                            if score > best_score:
                                best_score = score

                for read_alignment in best_alignments:
                    for e in read_alignment:
                        if e.startswith('AS:i:'):
                            score = int(e.split(':')[-1])
                            if score == best_score:
                                best_alignments_by_as.append(read_alignment)

                if len(best_alignments_by_as) == 1:
                    counter[3] += 1
                    result_gaf_str += '\t'.join(best_alignments_by_as[0]) + '\n'
                    continue

                # Do alignment share any common segments? If they do, it is multi-path alignment.
                # And I do not consider them as multimapping
                shared_segments = set()
                init = True
                for read_alignment in best_alignments_by_as:
                    path_str = read_alignment[5]
                    path = mcall.alignment_path_parse(path_str)
                    # print(path)

                    if init:
                        shared_segments = set(path[0])
                        init = False

                    shared_segments = shared_segments.intersection(set(path[0]))

                if len(shared_segments) > 0:
                    result_gaf_str += '\t'.join(best_alignments_by_as[0]) + '\n'
                    counter[3] += 1
                    continue

                # Hmm, still multimapping. Let's just output one of them.
                result_gaf_str += '\t'.join(best_alignments_by_as[0]) + "\tmp:i:1" + '\n'
                counter[4] += 1

    # counter2 = counter[:]
    # for ic in range(len(counter2)):
    #    counter2[ic] = counter2[ic] / counter[0] * 100

    return result_gaf_str


def tmp_gaf_processing_worker(pid, input_queue, result_queue):
    while True:
        tmp_gaf_fp = input_queue.get()
        if tmp_gaf_fp is None:
            result_queue.put(None)
            break

        result_gaf_str = tmp_gaf_processing(tmp_gaf_fp)
        result_queue.put(result_gaf_str)
    return

def merger_gaf_writer_worker(pid, worker_num, result_queue, output_fp):
    counter = 0
    with open(output_fp, 'w') as fh:
        while True:
            result_gaf_str = result_queue.get()
            if result_gaf_str is None:
                counter += 1
                if counter == worker_num:
                    break
                continue
            fh.write(result_gaf_str)
    return


def alignment_merge_main(working_dir, worker_num=20):
    output_gaf_fp = os.path.join(working_dir, 'alignment.gaf')

    input_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()

    pool = []
    for i in range(worker_num):
        p = multiprocessing.Process(target=tmp_gaf_processing_worker, args=(i, input_queue, result_queue))
        p.start()
        pool.append(p)

    writer_worker = multiprocessing.Process(target=merger_gaf_writer_worker, args=(0, worker_num, result_queue, output_gaf_fp))
    writer_worker.start()

    for i in range(tmp_alignment_file_count):
        input_gaf_fp = os.path.join(working_dir, f'alignment.{i}.gaf')
        input_queue.put(input_gaf_fp)

    for i in range(worker_num):
        input_queue.put(None)

    for p in pool:
        p.join()

    writer_worker.join()

    return


def alignment_main(fq1, fq2, work_dir, index_prefix, compress=True, thread=1, directional=True, vg_path="vg"):
    utility.fastq_converter(fq1, fq2, work_dir,
                            compress=compress,
                            thread=thread,
                            directional=directional,
                            split_num=tmp_alignment_file_count)

    alignment(work_dir=work_dir,
              index_prefix=index_prefix,
              output_format="gaf",
              thread=thread,
              directional=directional,
              compress=compress,
              vg_path=vg_path)

    alignment_merge_main(work_dir, worker_num=thread)

    # Just to wait a bit for alignment_merge_main to finish and garbage collection
    time.sleep(5)

    alignment_clenup(work_dir)

    return


if __name__ == '__main__':
    working_dir = sys.argv[1]
    alignment_merge_main(working_dir, worker_num=20)
    sys.exit(0)



































































































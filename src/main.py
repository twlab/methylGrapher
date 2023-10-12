

__author__ = "Wenjin Zhang"
__copyright__ = "Copyright 2023, Ting Wang Lab"
__credits__ = ["Juan Macias"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Wenjin Zhang"
__email__ = "wenjin@wustl.edu"


import os
import sys

import gfa
import mcall
import utility


conversion_types = ["C2T", "G2A"]


# TODO
# 1. Add fastuniq
# Remove OLD unused code


help_txt = """
Usage: python main.py <command> <arguments>
Commands:
    help
    PrepareGenome
    PrepareLibrary
    Align
    MethylCall

Help:
    python main.py help

PrepareGenome:
    It adds lambda phage genome to your genome graph, converts a GFA file into fully G->A and C->T converted GFA file, and indexes it for vg giraffe alignment.
    python main.py PrepareGenome 
    # Input options
    -gfa <gfa_file_path> 
    -lp <lambda_phage_genome_path> 
    # Output options
    -prefix <output_prefix> 
    -compress <Y/N>
    # Computing options
    -t <number_of_thread(s)> 

PrepareLibrary:
    Attention: The user should run Trim Glore first. 
    It first deduplicates your BS library (FASTQ file(s)), and then convert them into fully G->A and C->T converted FASTQ file.
    For single-end reads, just provide FASTQ file path to -fq1 argument.
    python main.py PrepareLibrary 
    # Input options
    -fq1 <fastq_file_path> 
    -fq2 <fastq_file_path> 
    # Output options
    -work_dir <work_directory> 
    -compress <Y/N> (default: Y)
    # Computing options
    -t <number_of_thread(s)> 
    -directional <Y/N> (default: Y) 

Align:
    VG Giraffe alignment, please provide work directory and index prefix.
    python main.py Align 
    -index_prefix <prefix> 
    -work_dir <work_directory> 
    -directional <Y/N> (default: Y)

MethylCall:
    Methylation call from vg giraffe alignment result.
    python main.py MethylCall 
    -work_dir <work_directory>

""".strip()


if __name__ == "__main__":
    args = sys.argv
    args.pop(0)

    # Find config file
    config = utility.ConfigParser("config.ini")
    vg_path = config.get("default", "vg_path")
    thread = config.get("default", "thread")
    try:
        thread = int(thread)
    except:
        thread = 1

    if vg_path in ["", None]:
        vg_path = "vg"


    if len(args) == 0:
        print("No command specified. Use 'help' for more information.")
        sys.exit(0)

    command = args.pop(0)
    command = command.lower()
    kvargs = {}
    while len(args) > 0:
        arg = args.pop(0)
        if arg.startswith("-"):
            key = arg[1:]
            value = args.pop(0)
            kvargs[key] = value
        else:
            print(f"Unknown argument: {arg}")
            sys.exit(1)

    if "thread" in kvargs:
        thread = int(kvargs["thread"])
        del kvargs["thread"]
    if "t" in kvargs:
        thread = int(kvargs["t"])
        del kvargs["t"]


    if command in ["help", "-h", "--help"]:
        print(help_txt)
        sys.exit(0)



    if command == "preparegenome":
        original_gfa_file_path = kvargs["gfa"]
        lambda_ref = kvargs["lp"]
        prefix = kvargs["prefix"]

        fn_report = prefix + ".prepare.genome.report.txt"
        freport = open(fn_report, "w")

        original_gfa_with_lambda = prefix + ".wl.gfa"
        gfa_c2t = prefix + ".wl.C2T.gfa"
        gfa_g2a = prefix + ".wl.G2A.gfa"

        lambda_segment_id = gfa.add_lambda_genome_to_gfa(original_gfa_file_path, original_gfa_with_lambda, lambda_ref)
        freport.write(f"Insert lambda phage genome into genome graph as segment: {lambda_segment_id}\n")


        g = gfa.GraphicalFragmentAssemblyMemory()
        g.parse(original_gfa_with_lambda, keep_link=True)
        g.write_converted(gfa_c2t, "C", "T")

        del g

        g = gfa.GraphicalFragmentAssemblyMemory()
        g.parse(original_gfa_with_lambda, keep_link=True)
        g.write_converted(gfa_g2a, "G", "A")

        del g

        for gfa_file_path in [gfa_c2t, gfa_g2a]:
            index_prefix = gfa_file_path[:-4]
            cmd = f"{vg_path} autoindex -g {gfa_file_path} -p {index_prefix} -w giraffe -t {thread} >> {fn_report} 2>&1"
            # os.system(cmd)

            se = utility.SystemExecute()
            fout, flog = se.execute(cmd, stdout=None, stderr=None)
            for line in fout:
                line = line.decode("utf-8")
                freport.write(line)
            for line in flog:
                line = line.decode("utf-8")
                freport.write(line)
            se.wait()

        freport.close()
        sys.exit(0)



    # Convert genome graph (in gfa format) to full C->T and G->A converted gfa file.
    if command == "preparelibrary":
        fq1 = kvargs["fq1"]
        fq2 = kvargs.get("fq2", None)
        work_dir = kvargs["work_dir"]

        compress = kvargs.get("compress", "Y")
        if compress.lower() == "y":
            compress = True
        else:
            compress = False

        directional = kvargs.get("directional", "Y")
        if directional.lower() == "y":
            directional = True
        else:
            directional = False
        utility.fastq_converter(fq1, fq2, work_dir, compress=compress, thread=thread, directional=directional)


        sys.exit(0)


    # Align converted FASTQ files to converted GFA files.
    if command == "align":
        work_dir = kvargs["work_dir"]
        index_prefix = kvargs["index_prefix"]

        index_prefix_ct = index_prefix + ".wl.C2T"
        index_prefix_ga = index_prefix + ".wl.G2A"

        read1_alignment_type = ["C2T"]
        read2_alignment_type = ["G2A"]

        if kvargs.get("directional", "y") in "nN":
            read1_alignment_type = conversion_types
            read2_alignment_type = conversion_types


        output_format = "gaf"

        alignment_file_path = []
        alignment_outs = {}
        for i in range(10):
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

                    fq1 = f"{work_dir}/{read_type1}.R1.fastq.gz"
                    fq2 = f"{work_dir}/{read_type2}.R2.fastq.gz"

                    if ref_type == "G2A":
                        index_prefix = index_prefix_ga
                    else:
                        index_prefix = index_prefix_ct

                    giraffe_input = f"-f {fq1}"
                    if os.path.exists(fq2):
                        giraffe_input += f" -f {fq2}"


                    alignment_log = f"{work_dir}/alignment.Ref_{ref_type}.R1_{read_type1}.R2_{read_type2}.log"

                    cmd = f"vg giraffe -p -t {thread} -o {output_format} -M 2 --named-coordinates -Z {index_prefix}.giraffe.gbz -m {index_prefix}.min -d {index_prefix}.dist {giraffe_input}"
                    # print(cmd)
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
                        ind = int(l[0][-3])
                        alignment_out_fh = alignment_outs[ind]
                        alignment_out_fh.write(line)
                    se.wait()

        for fh in alignment_outs.values():
            fh.close()


        for fn in alignment_file_path:
            cmd = f"sort -o {fn} {fn}"

            se = utility.SystemExecute()
            fout, flog = se.execute(cmd, stdout=None, stderr=None)
            se.wait()

        sys.exit(0)



    # Call methylation
    if command == "methylcall":
        work_dir = kvargs["work_dir"]
        gfa_file = kvargs["index_prefix"] + ".wl.gfa"

        # TODO make them parameters
        minimum_identity = 20
        minimum_mapq = 0
        discard_multimapped = True

        alignment_fp = []
        for i in range(10):
            alignment_fp.append(f"{work_dir}/alignment.{i}.gaf")

        gfa_instance = gfa.GraphicalFragmentAssemblyMemory()
        gfa_instance.parse(gfa_file)

        methylation_tmp_file_names = []
        for i in range(0, 10):
            fn = f"{work_dir}/mcall.{i}.tmp"
            methylation_tmp_file_names.append(fn)

        methylation_tmp_output = []
        for fn in methylation_tmp_file_names:
            methylation_tmp_output.append(open(fn, "w"))

        total_alignment_count, low_confidence_count, error_count, multimapped_count = mcall.call(
            gfa_instance,
            alignment_fp,
            methylation_tmp_output,
            minimum_identity=minimum_identity,
            minimum_mapq=minimum_mapq,
            discard_multimapped=discard_multimapped
        )
        freport = open(f"{work_dir}/report.txt", "a")
        freport.write(f"Total aligned reads: {total_alignment_count}\n")
        freport.write(f"Total multi-mapped reads: {multimapped_count}\n")
        freport.write(f"Low confidence count (low quality alignment in terms of block length and MapQ): {low_confidence_count}\n")
        freport.write(f"Inconsistent alignment entries: {error_count}\n")

        del gfa_instance

        for fh in methylation_tmp_output:
            fh.close()

        methylation_tmp_output = []
        for fn in methylation_tmp_file_names:
            methylation_tmp_output.append(open(fn))
        final_out = open("graph2.methyl", "w")

        for fn in methylation_tmp_output:
            mcall.mcall_tmp_to_final(fn, final_out)


        for fn in methylation_tmp_file_names:
            os.remove(fn)

        sys.exit(0)

    if command == "qc":
        import postprocessingQC


        # TODO run QC
        work_dir = kvargs["work_dir"]
        lp_id = kvargs["lp_id"]

        mcalled = f"{work_dir}/graph.methyl"

        r = postprocessingQC.estimate_conversion_rate(mcalled, lp_id)
        print(f"Lambda Phage Segment ID: {lp_id}")
        print(f"Conversion rate: {r*100:.2f}%")

        sys.exit(0)











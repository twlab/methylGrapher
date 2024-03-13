

__author__ = "Wenjin Zhang"
__copyright__ = "Copyright 2023, Ting Wang Lab"
__credits__ = ["Juan Macias"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = "Wenjin Zhang"
__email__ = "wenjin@wustl.edu"


import os
import sys

import gfa
import mcall
import utility


# Global variables
conversion_types = ["C2T", "G2A"]



help = utility.HelpDocument()





def alignment(work_dir="./", index_prefix="", output_format="gaf", thread=1, directional=True):

    index_prefix_ct = index_prefix + ".wl.C2T"
    index_prefix_ga = index_prefix + ".wl.G2A"

    read1_alignment_type = ["C2T"]
    read2_alignment_type = ["G2A"]

    if not directional:
        read1_alignment_type = conversion_types
        read2_alignment_type = conversion_types

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

    return







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
        print(help.help_text())
        sys.exit(0)

    command = args.pop(0)
    command = command.lower()
    if command not in ["preparegenome", "preparelibrary", "align", "methylcall", "qc", "help", "-h", "--help", "simple"]:
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

        original_gfa_with_lambda = prefix + ".wl.gfa"
        gfa_c2t = prefix + ".wl.C2T.gfa"
        gfa_g2a = prefix + ".wl.G2A.gfa"


        lambda_ref = kvargs.get("lp", None)
        lambda_segment_id = gfa.add_lambda_genome_to_gfa(original_gfa_file_path, original_gfa_with_lambda, lambda_ref)

        if lambda_ref is None:
            freport.write("NO lambda phage genome provided. \nNOT able to estimate conversion rate.")
        else:
            freport.write(f"Insert lambda phage genome into genome graph as segment: {lambda_segment_id}\n")


        graph_trim_flag = kvargs.get("trim", "Y").lower() in "yestrue"
        g = gfa.GraphicalFragmentAssemblyMemory()
        g.parse(original_gfa_with_lambda, keep_link=True)
        g.write_converted(gfa_c2t, "C", "T", SNV_trim=graph_trim_flag)

        del g

        g = gfa.GraphicalFragmentAssemblyMemory()
        g.parse(original_gfa_with_lambda, keep_link=True)
        g.write_converted(gfa_g2a, "G", "A", SNV_trim=graph_trim_flag)

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
        work_dir = kvargs.get("work_dir", "./")

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
        work_dir = kvargs.get("work_dir", "./")
        index_prefix = kvargs["index_prefix"]
        output_format = "gaf"

        directional = kvargs.get("directional", "Y")
        if directional.lower() == "y":
            directional = True
        else:
            directional = False

        alignment(work_dir=work_dir, index_prefix=index_prefix, output_format=output_format, thread=thread, directional=directional)

        sys.exit(0)



    # Call methylation
    if command == "methylcall":
        work_dir = kvargs.get("work_dir", "./")
        gfa_file = kvargs["index_prefix"] + ".wl.gfa"

        minimum_identity = 20
        minimum_mapq = 0
        discard_multimapped = True

        if "minimum_identity" in kvargs:
            minimum_identity = int(kvargs["minimum_identity"])
        if "minimum_mapq" in kvargs:
            minimum_mapq = int(kvargs["minimum_mapq"])
        if "discard_multimapped" in kvargs:
            discard_multimapped = kvargs["discard_multimapped"].lower() in "yestrue"

        batch_size = 4096
        if "batch_size" in kvargs:
            batch_size = int(kvargs["batch_size"])

        assert minimum_identity >= 0
        assert minimum_mapq >= 0

        gfa_worker_num = 1
        if thread > 20:
            gfa_worker_num = 2

        mcall.mcall_main(
            work_dir, gfa_file,
            minimum_identity=minimum_identity, minimum_mapq=minimum_mapq, discard_multimapped=True,
            process_count=thread, alignment_parse_worker_num=1, gfa_worker_num=gfa_worker_num, batch_size=batch_size
        )

        sys.exit(0)


    # One command to run all
    if command == "main":

        work_dir = kvargs.get("work_dir", "./")

        fq1 = kvargs["fq1"]
        fq2 = kvargs.get("fq2", None)

        index_prefix = kvargs["index_prefix"]

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

        alignment_output_format = "gaf"

        minimum_identity = 20
        minimum_mapq = 0
        discard_multimapped = True

        if "minimum_identity" in kvargs:
            minimum_identity = int(kvargs["minimum_identity"])
        if "minimum_mapq" in kvargs:
            minimum_mapq = int(kvargs["minimum_mapq"])
        if "discard_multimapped" in kvargs:
            discard_multimapped = kvargs["discard_multimapped"].lower() in "yestrue"
        assert minimum_identity >= 0
        assert minimum_mapq >= 0



        batch_size = 4096
        if "batch_size" in kvargs:
            batch_size = int(kvargs["batch_size"])

        gfa_worker_num = 1
        if thread > 20:
            gfa_worker_num = 2

        utility.fastq_converter(fq1, fq2, work_dir, compress=compress, thread=thread, directional=directional)

        alignment(work_dir=work_dir, index_prefix=index_prefix, output_format=alignment_output_format, thread=thread, directional=directional)

        mcall.mcall_main(
            work_dir, index_prefix + ".wl.gfa",
            minimum_identity=minimum_identity, minimum_mapq=minimum_mapq, discard_multimapped=True,
            process_count=thread, alignment_parse_worker_num=1, gfa_worker_num=gfa_worker_num, batch_size=batch_size
        )

        for ct in conversion_types:
            for ri in ["R1", "R2"]:
                fn = f"{work_dir}/{ct}.{ri}.fastq.gz"
                if os.path.exists(fn):
                    os.remove(fn)

        sys.exit(0)












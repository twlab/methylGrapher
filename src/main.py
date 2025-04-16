#!/usr/bin/env python

__author__ = "Wenjin Zhang"
__copyright__ = "Copyright 2023-2025, Ting Wang Lab"
__credits__ = ["Juan Macias"]
__license__ = "MIT"
__version__ = "0.2.0"
__maintainer__ = "Wenjin Zhang"
__email__ = "wenjin@wustl.edu"



import os
import sys

import gfa
import json
import mcall
import utility
import alignments

# Global variables
conversion_types = ["C2T", "G2A"]

help = utility.HelpDocument()



if __name__ == "__main__":
    args = sys.argv
    args.pop(0)

    vg_path = "vg"
    thread = 1

    # Find config file
    try:
        config = utility.ConfigParser("config.ini")
        vg_path = config.get("default", "vg_path")
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
    if command not in ["preparegenome", "align", "methylcall", "conversionrate", "mergecpg", "help", "-h", "--help", "main", "mergegaf", "vg_check"]:
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
    if "vg_path" in kvargs:
        vg_path = kvargs["vg_path"]

    if command in ["help", "-h", "--help"]:
        print(help.help_text())
        sys.exit(0)

    # Convert genome graph (in gfa format) to full C->T and G->A converted gfa file.
    if command == "preparegenome":
        original_gfa_file_path = kvargs["gfa"]
        prefix = kvargs["prefix"]

        fn_report = prefix + ".prepare.genome.report.txt"
        freport = open(fn_report, "w")

        found, vg_log = utility.vg_binary_check(vg_path=vg_path)
        freport.write(vg_log+"\n\n")

        original_gfa_with_lambda = prefix + ".wl.gfa"
        gfa_c2t = prefix + ".wl.C2T.gfa"
        gfa_g2a = prefix + ".wl.G2A.gfa"
        graph_all_cpg_fp = prefix + ".cpg.tsv"

        lambda_ref = kvargs.get("lp", None)
        lambda_segment_id = gfa.add_lambda_genome_to_gfa(original_gfa_file_path, original_gfa_with_lambda, lambda_ref)

        if lambda_ref is None:
            freport.write("NO lambda phage genome provided. \nNOT able to estimate conversion rate.")
        else:
            freport.write(f"Insert lambda phage genome into genome graph as segment: {lambda_segment_id}\n")

        graph_trim_flag = kvargs.get("trim", "N").lower() in "yestrue"


        if graph_trim_flag:
            g = gfa.GraphicalFragmentAssemblyMemory()
            g.parse(original_gfa_with_lambda, keep_link=True)
            g.write_converted(gfa_c2t, "C", "T", SNV_trim=graph_trim_flag)

            del g

            g = gfa.GraphicalFragmentAssemblyMemory()
            g.parse(original_gfa_with_lambda, keep_link=True)
            g.write_converted(gfa_g2a, "G", "A", SNV_trim=graph_trim_flag)

            del g

        else:
            utility.gfa_converter(original_gfa_with_lambda, prefix+".wl", compress=False)


        g = gfa.GraphicalFragmentAssemblyMemory()
        g.parse(original_gfa_with_lambda, keep_link=True)

        replacement_ndoes = {}
        replacement_ndoes["CT"] = g.get_replacement_SNV("C", "T")
        replacement_ndoes["GA"] = g.get_replacement_SNV("G", "A")

        with open(prefix + ".wl.node.replacement.json", "w") as nr_json_fh:
            json.dump(replacement_ndoes, nr_json_fh, indent=4)



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

        utility.get_all_cpg_from_graph(original_gfa_file_path, graph_all_cpg_fp)

        freport.close()
        sys.exit(0)

    # Alignment and related stuff...
    if command == "align":

        fq1 = kvargs["fq1"]
        fq2 = kvargs.get("fq2", None)

        work_dir = kvargs.get("work_dir", "./")
        index_prefix = kvargs["index_prefix"]
        output_format = "gaf"

        compress = kvargs.get("compress", "N")
        if compress.lower() == "y":
            compress = True
        else:
            compress = False

        directional = kvargs.get("directional", "Y")
        if directional.lower() == "y":
            directional = True
        else:
            directional = False

        alignments.alignment_main(fq1, fq2, work_dir, index_prefix,
                                  compress=compress,
                                  thread=thread,
                                  directional=directional,
                                  vg_path=vg_path
                                  )
        sys.exit(0)

    # Call methylation
    if command == "methylcall":
        work_dir = kvargs.get("work_dir", "./")
        # gfa_file = kvargs["index_prefix"] + ".wl.gfa"

        minimum_identity = 20
        minimum_mapq = 0
        discard_multimapped = True
        cg_only = True
        genotyping_cytosine = False


        if "minimum_identity" in kvargs:
            minimum_identity = int(kvargs["minimum_identity"])
        if "minimum_mapq" in kvargs:
            minimum_mapq = int(kvargs["minimum_mapq"])
        if "discard_multimapped" in kvargs:
            discard_multimapped = kvargs["discard_multimapped"].lower() in "yestrue"
        if "cg_only" in kvargs:
            cg_only = kvargs["cg_only"].lower() in "yestrue"
        if "genotyping_cytosine" in kvargs:
            genotyping_cytosine = kvargs["genotyping_cytosine"].lower() in "yestrue"

        batch_size = 4096
        if "batch_size" in kvargs:
            batch_size = int(kvargs["batch_size"])

        assert minimum_identity >= 0
        assert minimum_mapq >= 0

        gfa_worker_num = 1
        if thread > 20:
            gfa_worker_num = 2

        mcall.mcall_main(
            work_dir, kvargs["index_prefix"],
            cg_only=cg_only, genotyping_cytosine=genotyping_cytosine,
            minimum_identity=minimum_identity, minimum_mapq=minimum_mapq, discard_multimapped=discard_multimapped,
            process_count=thread, alignment_parse_worker_num=1, gfa_worker_num=gfa_worker_num, batch_size=batch_size
        )

        sys.exit(0)

    # Conversion rate estimation
    if command == "conversionrate":
        work_dir = kvargs.get("work_dir", "./")
        index_prefix = kvargs["index_prefix"]

        s = utility.estimate_conversion_rate_print(index_prefix, work_dir)
        print(s)
        sys.exit(0)

    # merge CpG
    if command == "mergecpg":
        work_dir = kvargs.get("work_dir", "./")
        index_prefix = kvargs["index_prefix"]

        graph_all_cpg_fp = index_prefix + ".cpg.tsv"
        cytosine_fp = work_dir + "/graph.methyl"
        genotype_cytosine_fp = work_dir + "/genotype.info.txt"
        merged_cytosine_in_cpg_context_fp = work_dir + "/graph.cpg.tsv"
        validated_cytosine_fp = work_dir + "/graph.validated.CG.methyl"

        utility.merge_graph_cytosines(
            graph_all_cpg_fp,
            cytosine_fp,
            genotype_cytosine_fp,
            merged_cytosine_in_cpg_context_fp,
            full_position=False
        )

    # One command to run all
    if command == "main":

        work_dir = kvargs.get("work_dir", "./")

        fq1 = kvargs["fq1"]
        fq2 = kvargs.get("fq2", None)

        index_prefix = kvargs["index_prefix"]

        compress = kvargs.get("compress", "N")
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
        cg_only = True
        genotyping_cytosine = False

        if "minimum_identity" in kvargs:
            minimum_identity = int(kvargs["minimum_identity"])
        if "minimum_mapq" in kvargs:
            minimum_mapq = int(kvargs["minimum_mapq"])
        if "discard_multimapped" in kvargs:
            discard_multimapped = kvargs["discard_multimapped"].lower() in "yestrue"
        if "cg_only" in kvargs:
            cg_only = kvargs["cg_only"].lower() in "yestrue"
        if "genotyping_cytosine" in kvargs:
            genotyping_cytosine = kvargs["genotyping_cytosine"].lower() in "yestrue"

        assert minimum_identity >= 0
        assert minimum_mapq >= 0

        batch_size = 4096
        if "batch_size" in kvargs:
            batch_size = int(kvargs["batch_size"])

        gfa_worker_num = 1
        if thread > 20:
            gfa_worker_num = 2

        alignments.alignment_main(fq1, fq2, work_dir, index_prefix, compress=compress, thread=thread, directional=directional)

        mcall.mcall_main(
            work_dir, index_prefix,
            cg_only=cg_only, genotyping_cytosine=genotyping_cytosine,
            minimum_identity=minimum_identity, minimum_mapq=minimum_mapq, discard_multimapped=discard_multimapped,
            process_count=thread, alignment_parse_worker_num=1, gfa_worker_num=gfa_worker_num, batch_size=batch_size
        )

        sys.exit(0)


    # Check if vg is installed and in the path
    if command == "vg_check":
        found_vg_binary, log = utility.vg_binary_check(vg_path=vg_path)
        print(log)

        sys.exit(0)









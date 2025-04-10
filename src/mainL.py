

import os
import re
import sys
import time
import argparse


import utility
import longread








if __name__ == "__main__":

    utl = utility.Utility()


    # Arugment parser
    parser = argparse.ArgumentParser(description='methylGrapher for long reads')

    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')


    # main
    parser_main = subparsers.add_parser('main', help='One step to run all')

    # 1 prepare fasta
    parser_prepare_fasta = subparsers.add_parser('prepare_fasta', help='Prepare fasta file')

    # 2 align
    parser_align = subparsers.add_parser('align', help='Align long reads to genome graph')

    # 3 methylation extraction
    parser_extraction = subparsers.add_parser('extraction', help='Extract methylation information')


    all_parser = [parser_main, parser_prepare_fasta, parser_align, parser_extraction]
    for p in all_parser:
        p.add_argument('-t', help='number of threads', default=1, type=int)
        p.add_argument('-debug', help='debug mode')
        p.add_argument('-wd', help='working directory', default="./")
        p.add_argument('-verbose', help='verbosity')



    for p in [parser_main, parser_prepare_fasta]:
        p.add_argument('-basecall', help='Long read base call BAM/SAM file', required=True)

    for p in [parser_main, parser_align]:
        p.add_argument('-binary_path', help='Binary path for aligner', default="vg")
        p.add_argument('-aligner', help='Aligner, either vg or GraphAligner')
        p.add_argument('-additional_alignment_params', help='Additional alignment parameters')

    for p in [parser_main, parser_extraction]:
        p.add_argument('-discard_cg_mismatch', help='discard methylation information if read aligns to non CG sites')
        p.add_argument('-mapq_threshold', help='discard read alignment if mapq < threshold', default=10, type=int)

    for p in [parser_main, parser_align, parser_extraction]:
        p.add_argument('-gfa', help='Genome graph file path in GFA format', required=True)



    args = parser.parse_args()


    # Default parameters
    wd = "./"
    thread = 1
    debug = False
    discard_mismatched_cg = True
    verbose = False

    basecall_fp = None
    gfa_fp = None

    mapq_threshold = 10


    aligner_binary_path = "vg"
    aligner = "vg"
    additional_alignment_params = ""




    # Parse arguments
    if "wd" in args and args.wd:
        wd = args.wd
    if "t" in args and args.t:
        thread = int(args.t)
    if "debug" in args and args.debug:
        debug = True

    if "basecall" in args and args.basecall:
        basecall_fp = args.basecall

    if "gfa" in args and args.gfa:
        gfa_fp = args.gfa
        if not os.path.exists(gfa_fp):
            print("Genome graph file does not exist", file=sys.stderr)
            sys.exit(1)

    if "verbose" in args and args.verbose:
        verbose = utl.argument_boolean(args.verbose)

    if "aligner" in args and args.aligner:
        aligner = args.aligner
        if aligner not in ["vg", "GraphAligner"]:
            print("Aligner must be either vg or GraphAligner", file=sys.stderr)
            sys.exit(1)

    if "binary_path" in args and args.binary_path:
        aligner_binary_path = args.binary_path

    if "additional_alignment_params" in args and args.additional_alignment_params:
        additional_alignment_params = args.additional_alignment_params

    if "discard_cg_mismatch" in args and args.discard_cg_mismatch:
        discard_mismatched_cg = utl.argument_boolean(args.discard_cg_mismatch)

    if "mapq_threshold" in args and args.mapq_threshold:
        mapq_threshold = int(args.mapq_threshold)
        if mapq_threshold < 0:
            print("Mapq threshold must be greater than 0", file=sys.stderr)
            sys.exit(1)



    if debug:
        verbose = True


    if args.command == 'main':
        print("Running methylGrapherL with main function", file=sys.stderr)


        # Log parameters
        print("Long read base call file path: ", basecall_fp, file=sys.stderr)
        print("Genome graph file path: ", gfa_fp, file=sys.stderr)
        print("Working directory: ", wd, file=sys.stderr)
        print("Number of threads: ", thread, file=sys.stderr)
        print("Debug mode: ", debug, file=sys.stderr)
        print("Verbose: ", verbose, file=sys.stderr)
        print("Discard mismatched CG: ", discard_mismatched_cg, file=sys.stderr)

        # Run main function
        longread.main(basecall_fp, gfa_fp, wd, threads=thread, debug=debug, mapq_threshold=mapq_threshold, discard_cg_mismatch=discard_mismatched_cg, verbose=verbose)


    elif args.command == 'prepare_fasta':
        print("Prepare fasta", file=sys.stderr)

        # Log parameters
        print("Long read base call file path: ", basecall_fp, file=sys.stderr)
        print("Working directory: ", wd, file=sys.stderr)
        print("Number of threads: ", thread, file=sys.stderr)
        print("Debug mode: ", debug, file=sys.stderr)

        # Run prepare fasta function
        longread.prepare_fasta(basecall_fp, wd)



    elif args.command == 'align':
        print("Align", file=sys.stderr)

        # Log parameters
        print("Genome graph file path: ", gfa_fp, file=sys.stderr)
        print("Working directory: ", wd, file=sys.stderr)
        print("Number of threads: ", thread, file=sys.stderr)
        print("Debug mode: ", debug, file=sys.stderr)

        # Run align function
        # binary_path=None, aligner="GraphAligner", additional_alignment_params=""
        longread.align(gfa_fp, wd, threads=thread, debug=debug, binary_path=aligner_binary_path, aligner=aligner, additional_alignment_params=additional_alignment_params)



    elif args.command == 'extraction':

        # Log parameters
        print("Genome graph file path: ", gfa_fp, file=sys.stderr)
        print("Working directory: ", wd, file=sys.stderr)
        print("Number of threads: ", thread, file=sys.stderr)
        print("Debug mode: ", debug, file=sys.stderr)
        print("Verbose: ", verbose, file=sys.stderr)
        print("Discard mismatched CG: ", discard_mismatched_cg, file=sys.stderr)
        print("MAPQ threshold: ", mapq_threshold, file=sys.stderr)


        longread.extract_methylation(gfa_fp, wd, threads=1, debug=False, discard_cg_mismatch=discard_mismatched_cg, verbose=verbose, mapq_threshold=mapq_threshold)

    else:
        print('No command provided')




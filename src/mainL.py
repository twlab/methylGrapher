


import os
import re
import sys
import argparse


import longread








if __name__ == "__main__":

    # Arugment parser
    parser = argparse.ArgumentParser(description='methylGrapher for long reads')

    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')


    # main
    parser_main = subparsers.add_parser('main', help='One step to run all')
    parser_main.add_argument('-basecall', help='Long read base call BAM/SAM file', required=True)
    parser_main.add_argument('-gfa', help='Genome graph file path in GFA format', required=True)
    parser_main.add_argument('-t', help='number of threads', default=1, type=int)
    parser_main.add_argument('-debug', help='debug mode')
    parser_main.add_argument('-wd', help='working directory', default="./")


    # 1 prepare fasta
    parser_prepare_fasta = subparsers.add_parser('prepare_fasta', help='Prepare fasta file')
    parser_prepare_fasta.add_argument('-basecall', help='Long read base call BAM/SAM file', required=True)
    parser_prepare_fasta.add_argument('-wd', help='working directory', default="./")
    parser_prepare_fasta.add_argument('-t', help='number of threads', default=1, type=int)
    parser_prepare_fasta.add_argument('-debug', help='debug mode')


    # 2 align
    parser_align = subparsers.add_parser('align', help='Align long reads to genome graph')
    parser_align.add_argument('-gfa', help='Genome graph file path in GFA format', required=True)
    parser_align.add_argument('-wd', help='Long read fasta file', default="./")
    parser_align.add_argument('-t', help='number of threads', default=1, type=int)
    parser_align.add_argument('-debug', help='debug mode')

    # 3 methylation extraction
    parser_extraction = subparsers.add_parser('extraction', help='Extract methylation information')
    parser_extraction.add_argument('-gfa', help='Genome graph file path in GFA format', required=True)
    parser_extraction.add_argument('-wd', help='working directory', default="./")
    parser_extraction.add_argument('-t', help='number of threads', default=1, type=int)
    parser_extraction.add_argument('-debug', help='debug mode')


    # TODO add function execution for subcommands



    args = parser.parse_args()

    if args.command == 'main':
        print("Running methylGrapherL with main function", file=sys.stderr)


        wd = "./"
        thread = 1
        debug = False
        basecall_fp = args.basecall
        gfa_fp = args.gfa

        if args.wd:
            wd = args.wd
        if args.t:
            thread = int(args.t)
        if args.debug:
            debug = True



        # Log parameters
        print("Long read base call file path: ", basecall_fp, file=sys.stderr)
        print("Genome graph file path: ", gfa_fp, file=sys.stderr)
        print("Working directory: ", wd, file=sys.stderr)
        print("Number of threads: ", thread, file=sys.stderr)
        print("Debug mode: ", debug, file=sys.stderr)

        # Run main function
        longread.main(basecall_fp, gfa_fp, wd)


    elif args.command == 'prepare_fasta':
        print("Prepare fasta", file=sys.stderr)

        wd = "./"
        thread = 1
        debug = False
        basecall_fp = args.basecall

        if args.wd:
            wd = args.wd
        if args.t:
            thread = int(args.t)
        if args.debug:
            debug = True

        # Log parameters
        print("Long read base call file path: ", basecall_fp, file=sys.stderr)
        print("Working directory: ", wd, file=sys.stderr)
        print("Number of threads: ", thread, file=sys.stderr)
        print("Debug mode: ", debug, file=sys.stderr)

        # Run prepare fasta function
        longread.prepare_fasta(basecall_fp, wd)



    elif args.command == 'align':
        print("Align")

        wd = "./"
        thread = 1
        debug = False
        gfa_fp = args.gfa

        if args.wd:
            wd = args.wd
        if args.t:
            thread = int(args.t)
        if args.debug:
            debug = True


        # Log parameters
        print("Genome graph file path: ", gfa_fp, file=sys.stderr)
        print("Working directory: ", wd, file=sys.stderr)
        print("Number of threads: ", thread, file=sys.stderr)
        print("Debug mode: ", debug, file=sys.stderr)



    elif args.command == 'extraction':
        wd = "./"
        thread = 1
        debug = False
        gfa_fp = args.gfa

        if args.wd:
            wd = args.wd
        if args.t:
            thread = int(args.t)
        if args.debug:
            debug = True

        # Log parameters
        print("Genome graph file path: ", gfa_fp, file=sys.stderr)
        print("Working directory: ", wd, file=sys.stderr)
        print("Number of threads: ", thread, file=sys.stderr)
        print("Debug mode: ", debug, file=sys.stderr)

    else:
        print('No command provided')



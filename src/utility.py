
import os
import re
import sys
import gzip
import time
# import ctypes
import random
import string
import hashlib
import subprocess
import multiprocessing


atcg_complement_dict = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N"
}


phred_lookup_int = {chr(i): i - 33 for i in range(33, 255)}
def phred_to_int(qual_char):
    assert qual_char in phred_lookup_int
    return phred_lookup_int[qual_char]


phred_lookup_prob = {char: 10 ** (-(phred_score / 10.0)) for char, phred_score in phred_lookup_int.items()}
def phred_to_prob(qual_char):
    assert qual_char in phred_lookup_prob
    return phred_lookup_prob[qual_char]


def sequence_reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


# Read methylGrapher result
def graph_methyl_cytosine_reader(fp):
    res = {}

    with open(fp) as fh:
        for l in fh:
            l = l.strip().split("\t")
            segID, pos, strand, context, unmet, met, cov, ml = l

            met = int(met)
            cov = int(cov)

            if context != "CG":
                continue

            if segID not in res:
                res[segID] = {}

            res[segID][pos] = (met, cov)

    return res


class Utility(object):
    lambda_phage_ref = None

    def __init__(self):
        self._nuclotides_translate = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        try:

            cfp = ConfigParser("config.ini")
            cfp.find()
            cfp.parse()
            self.lambda_phage_ref = cfp.get("default", "lambda_phage_ref")
        except:
            pass

        return

    @staticmethod
    def random_sequence(length):
        # Just for testing purposes
        return ''.join([random.choice("atcg") for _ in range(length)])

    def seq_reverse_complement(self, seq):
        return seq.translate(self._nuclotides_translate)[::-1]

    @staticmethod
    def isGzip(file_path):
        res = False
        for ext in [".gz", ".gzip", ".GZ", ".GZIP"]:
            if file_path.endswith(ext):
                res = True

        return res

    @staticmethod
    def get_random_sha1():
        rstr = ''.join(random.choices(string.ascii_lowercase + string.digits, k=100))
        res = hashlib.sha1(rstr.encode('utf-8')).hexdigest()
        return res

    @staticmethod
    def get_random_md5():
        rstr = ''.join(random.choices(string.ascii_lowercase + string.digits, k=100))
        res = hashlib.md5(rstr.encode('utf-8')).hexdigest()
        return res

    @staticmethod
    def bs_convert(read1, error_rate, methylation_level, convert_base="C"):
        assert error_rate >= 0 and error_rate <= 1
        assert methylation_level >= 0 and methylation_level <= 1
        assert convert_base in ["C", "G"]
        res = ""
        for nt in read1:
            new_nt = nt
            if random.random() < error_rate:
                new_nt = random.choice(['A', 'T', 'C', 'G', ""])

            if new_nt == convert_base:
                if random.random() > methylation_level:
                    new_nt = "T" if new_nt == "C" else "A"

            res += new_nt
        return res

    @staticmethod
    def argument_boolean(arg):
        if arg.lower() in ['true', 't', 'yes', 'y']:
            return True
        elif arg.lower() in ['false', 'f', 'no', 'n']:
            return False
        else:
            raise Exception(f"Invalid boolean argument: {arg}")


def gfa_converter(input_gfa, output_file_prefix, compress=True, thread=1):
    utils = Utility()

    for conversion_str in ['C2T', 'G2A']:
        conversion = conversion_str.split('2')

        output_gfa = f"{output_file_prefix}.{conversion_str}.gfa"

        output_file_handle = None
        if compress:
            output_gfa += '.gz'
            output_file_handle = gzip.open(output_gfa, 'wt')
        else:
            output_file_handle = open(output_gfa, 'w')

        input_compress = False

        if utils.isGzip(input_gfa):
            input_compress = True

        # print(thread)
        if input_compress:
            input_file_handle = gzip.open(input_gfa, 'rt')
        else:
            input_file_handle = open(input_gfa, 'r')

        for l in input_file_handle:
            l = l.strip().split('\t')

            if l[0] == 'S':
                sc = l[2].upper().replace(conversion[0], conversion[1])
                l[2] = sc

            newl = '\t'.join(l) + '\n'

            output_file_handle.write(newl)

    return 0


class ConfigParser(object):

    def __init__(self, filename, search=True):
        self._path = filename
        self._data = {}
        if search:
            self.find()
        self.parse()

    def parse(self):
        with open(self._path, "r") as f:
            block_name = ""
            key = ""
            for line in f:
                line = line.strip()
                if len(line) == 0 or line.startswith("#"):
                    continue

                if line.startswith("["):
                    block_name = line[1:-1]
                    self._data[block_name] = {}
                    continue
                else:
                    if "=" in line:
                        key, value = line.split("=")
                        key = key.strip()
                        value = value.strip()
                        self._data[block_name][key] = value
                    else:
                        self._data[block_name][key] = None
        return

    def get(self, block_name, key):
        res = None
        try:
            res = self._data[block_name][key]
        except KeyError:
            pass
        return res

    def find(self):
        for p in ["./", "~", os.path.dirname(__file__)]:
            path = os.path.join(p, self._path)
            if os.path.exists(path):
                # print("Find config file: %s" % path)
                self._path = path


class SystemExecute(object):

    def __init__(self):
        self._pool = []

    def execute(self, cmd, stdout=None, stderr=None):
        if stdout is None:
            stdout = subprocess.PIPE
        else:
            if not os.path.exists(os.path.dirname(stdout)):
                os.makedirs(os.path.dirname(stdout))
            stdout = open(stdout, "a")

        if stderr is None:
            stderr = subprocess.PIPE
        else:
            if not os.path.exists(os.path.dirname(stderr)):
                os.makedirs(os.path.dirname(stderr))
            stderr = open(stderr, "a")

        sp = subprocess.Popen(cmd, shell=True, stdout=stdout, stderr=stderr)
        self._pool.append(sp)

        stdout, stderr = sp.stdout, sp.stderr
        return stdout, stderr

    def wait(self):
        for sp in self._pool:
            sp.wait()
        self._pool = []


class GFASimple:

    def __init__(self):
        self.clear()

    def clear(self):
        self._segment = {}
        self._hg38_segments = set()

        self._link_forward = {}
        self._link_backward = {}

        self._segment_offset = {}

    def parse(self, gfa_file):

        gfa_file_handle = None
        if gfa_file.endswith(".gz"):
            gfa_file_handle = gzip.open(gfa_file, "rt")
        else:
            gfa_file_handle = open(gfa_file)

        for l in gfa_file_handle:
            if l[0] not in "SL":
                continue

            l = l.strip().split("\t")

            if l[0] == "S":
                rt, segID, seq, *tags = l
                ishg38 = False
                chrom, offset = None, None
                for t in tags:
                    if t.startswith("SN:Z:GRCh38."):
                        ishg38 = True
                        chrom = t[12:]
                    if t.startswith("SO:i:"):
                        offset = int(t[5:])

                if len(tags) > 2:
                    self._segment_offset[segID] = (chrom, offset)
                if ishg38:
                    self._hg38_segments.add(segID)

                self._segment[segID] = seq.upper()


            elif l[0] == "L":
                rt, segID1, dir1, segID2, dir2, *tags = l
                if segID1 not in self._link_forward:
                    self._link_forward[segID1] = {}
                self._link_forward[segID1][segID2] = (dir1, dir2)

        gfa_file_handle.close()

        return


def get_all_cpg_from_graph(gfa_fp, out_fp):
    gfa_simple_object = GFASimple()
    gfa_simple_object.parse(gfa_fp)

    cpg_list_file_path = out_fp

    cpg_list_file_handle = open(cpg_list_file_path, "w")

    within_seg_cpg_index = 0
    within_seg_cytosines = set()
    for segID in gfa_simple_object._segment:
        seq = gfa_simple_object._segment[segID]

        for pos in range(len(seq)):
            context = seq[pos:pos + 2]
            if context != "CG":
                continue

            cpg_type = "Other"
            isref = segID in gfa_simple_object._hg38_segments
            if isref:
                cpg_type = "hg38"

                chrom, offset = gfa_simple_object._segment_offset[segID]
                hg38_pos = f"({chrom}:{offset + pos})"
                cpg_type = cpg_type + hg38_pos

            line = [within_seg_cpg_index, segID, pos, segID, pos + 1, cpg_type]
            line_str = "C" + "\t".join(map(str, line)) + "\n"
            cpg_list_file_handle.write(line_str)

            p1 = (segID, pos)
            p2 = (segID, pos + 1)

            within_seg_cpg_index += 1
            within_seg_cytosines.add(p1)
            within_seg_cytosines.add(p2)

    agct_reverse_complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A"
    }

    aaa, bbb = 0, 0
    edge_cpg_index = 0
    edge_cytosines = set()
    for segID1 in gfa_simple_object._link_forward:
        for segID2 in gfa_simple_object._link_forward[segID1]:
            dir1, dir2 = gfa_simple_object._link_forward[segID1][segID2]

            seq1 = gfa_simple_object._segment[segID1]
            seq2 = gfa_simple_object._segment[segID2]

            context1 = seq1[-1]
            context2 = seq2[0]

            pos1 = len(seq1) - 1
            pos2 = 0

            if dir1 == "-":
                if seq1[0] not in "ACGT":
                    continue
                context1 = agct_reverse_complement[seq1[0]]
                pos1 = 0

            if dir2 == "-":
                if seq2[-1] not in "ACGT":
                    continue
                context2 = agct_reverse_complement[seq2[-1]]
                pos2 = len(seq2) - 1

            p1 = (segID1, pos1)
            p2 = (segID2, pos2)

            if context1 + context2 != "CG":
                continue

            # Debugging
            if False:
                print(segID1, segID2, dir1, dir2, len(seq1), len(seq2))
                print("\t", seq1, seq2)
                print(p1, p2)
                print()

            isref1 = segID1 in gfa_simple_object._hg38_segments
            isref2 = segID2 in gfa_simple_object._hg38_segments

            cpg_type = "SV"
            if isref1 and isref2:
                cpg_type = "hg38"
            elif isref1 or isref2:
                cpg_type = "SNV"

            if cpg_type == "hg38":

                aaa += 1
                chrom1, offset1 = gfa_simple_object._segment_offset[segID1]
                chrom2, offset2 = gfa_simple_object._segment_offset[segID2]
                assert chrom1 == chrom2
                chrom = chrom1
                hg38_pos1 = offset1 + pos1
                hg38_pos2 = offset2 + pos2
                hg38_pos = f"({chrom}:{min([hg38_pos1, hg38_pos2])})"
                diff = hg38_pos1 - hg38_pos2

                if abs(diff) > 1:
                    bbb += 1
                    # print(abs(diff), chrom, offset1, pos1, dir1, offset2, pos2, dir2, bbb / (aaa+bbb))
                # hg38_pos =
                cpg_type = cpg_type + hg38_pos

            line = [edge_cpg_index, segID1, pos1, segID2, pos2, cpg_type]
            line_str = "E" + "\t".join(map(str, line)) + "\n"
            cpg_list_file_handle.write(line_str)

            edge_cpg_index += 1
            edge_cytosines.add(p1)
            edge_cytosines.add(p2)

    cpg_list_file_handle.close()


def merge_graph_cytosines(cpg_fp, cytosine_fp, genotype_cytosine_fp, out_fp, full_position=False):
    cytosine_data = graph_methyl_cytosine_reader(cytosine_fp)

    false_positive_cytosines = set()
    if os.path.exists(genotype_cytosine_fp):
        with open(genotype_cytosine_fp) as fh:
            for l in fh:
                l = l.strip().split("\t")
                segID, pos, *others = l
                false_positive_cytosines.add((segID, pos))


    with open(out_fp, "w") as fh_out:
        for l in open(cpg_fp):
            cpg_ind, segID1, pos1, segID2, pos2, *tags = l.strip().split("\t")
            met = 0
            cov = 0

            cytosines = [(segID1, pos1), (segID2, pos2)]

            skip = False
            for cytosine in cytosines:
                if cytosine in false_positive_cytosines:
                    print("SKIP")
                    skip = True
            if skip:
                continue

            for segID, pos in cytosines:
                if segID not in cytosine_data:
                    continue
                if pos not in cytosine_data[segID]:
                    continue

                m, c = cytosine_data[segID][pos]
                met += m
                cov += c

            if cov == 0:
                continue

            line = [cpg_ind, met, cov]
            if full_position:
                line = [cpg_ind, segID1, pos1, segID2, pos2, met, cov]
            line_str = "\t".join(map(str, line)) + "\n"
            fh_out.write(line_str)

    return


def fastq_converter_worker_function(input_fastq_fp, output_fastq_fp, conversion_str, read_counts, split_num):
    utils = Utility()

    conversion = conversion_str.split('2')

    input_fastq_fh = None
    input_compress = utils.isGzip(input_fastq_fp)
    if input_compress:
        input_fastq_fh = gzip.open(input_fastq_fp, 'rt')
    else:
        input_fastq_fh = open(input_fastq_fp, 'r')

    output_fastq_fh = None
    output_compress = utils.isGzip(output_fastq_fp)
    if output_compress:
        output_fastq_fh = gzip.open(output_fastq_fp, 'wt')
    else:
        output_fastq_fh = open(output_fastq_fp, 'w')

    seq = ""
    phred = ""
    original_query_name = ""

    for i, l in enumerate(input_fastq_fh):

        if i % 4 == 0:
            if i == 0:
                assert l.startswith('@')
                original_query_name = l.strip()
                continue

            if " " in original_query_name:
                original_qn1 = original_query_name.split(' ', 1)[0]
            else:
                original_qn1 = original_query_name

            reminder = int(i / 4) % split_num
            converted_seq = seq.replace(conversion[0], conversion[1])

            newl = f"{original_qn1}_{conversion_str}_{reminder}_{seq}\n{converted_seq}\n+\n{phred}\n"
            output_fastq_fh.write(newl)

            seq = ""
            phred = ""
            original_query_name = ""

            assert l.startswith('@')
            original_query_name = l.strip()
        elif i % 4 == 1:
            seq = l.strip()
        elif i % 4 == 2:
            assert l.strip() == "+"
        elif i % 4 == 3:
            phred = l.strip().upper()

    # Edge case. The last read
    if " " in original_query_name:
        original_qn1 = original_query_name.split(' ', 1)[0]
    else:
        original_qn1 = original_query_name

    reminder = 1
    converted_seq = seq.replace(conversion[0], conversion[1])

    newl = f"{original_qn1}_{conversion_str}_{reminder}_{seq}\n{converted_seq}\n+\n{phred}\n"
    output_fastq_fh.write(newl)

    read_counts.append(int(i / 4) + 1)

    return 0


def fastq_converter(fq1, fq2, workdir, compress=True, thread=1, directional=True, split_num=1000):
    report_file_handle = open(f"{workdir}/report.txt", 'a')

    pool = []
    manager = multiprocessing.Manager()
    read_counts = manager.list()

    read_count = 0
    for fq in (fq1, fq2):
        read_count += 1

        if fq == None:
            continue

        for conversion_str in ['C2T', 'G2A']:

            if directional:
                if read_count == 1 and conversion_str == 'G2A':
                    continue
                if read_count == 2 and conversion_str == 'C2T':
                    continue

            input_fastq = fq
            output_fastq = f"{workdir}/{conversion_str}.R{read_count}.fastq"

            if compress:
                output_fastq += '.gz'

            p = multiprocessing.Process(
                target=fastq_converter_worker_function,
                args=(input_fastq, output_fastq, conversion_str, read_counts, split_num)
            )

            p.start()
            pool.append(p)

    for p in pool:
        p.join()

    # print(read_counts)
    assert len(set(read_counts)) == 1
    report_file_handle.write(f"Total reads for R1: {read_counts[0]}\n")
    if fq2 != None:
        report_file_handle.write(f"Total reads for R2: {read_counts[0]}\n\n")

    report_file_handle.close()

    return 0



def estimate_conversion_rate(index_prefix, work_dir):
    index_report_fp = f"{index_prefix}.prepare.genome.report.txt"
    # print(index_report_fp)
    lambda_segment_id = None
    with open(index_report_fp) as index_report_fh:
        for l in index_report_fh:
            if l.startswith("Insert lambda phage genome into genome graph as segment: "):
                lambda_segment_id = l.strip().split()[-1]
                break

    if lambda_segment_id == None:
        raise Exception(
            "Cannot find lambda phage segment id in the index report file. Did you run PrepareGenome with spike-in genome?")

    conversion_rate_by_context = {}
    graph_methyl_fp = f"{work_dir}/graph.methyl"
    # print(lambda_segment_id)
    with open(graph_methyl_fp) as graph_methyl_fh:
        for l in graph_methyl_fh:
            l = l.strip().split()
            if l[0] != lambda_segment_id:
                continue

            context = l[3]
            unmet = int(l[4])
            cov = int(l[6])

            if context not in conversion_rate_by_context:
                conversion_rate_by_context[context] = [0, 0]

            conversion_rate_by_context[context][0] += unmet
            conversion_rate_by_context[context][1] += cov

    res = {}
    unmet_total = 0
    cov_total = 0
    for context in conversion_rate_by_context:
        unmet, cov = conversion_rate_by_context[context]
        conversion_rate_by_context[context] = unmet / cov
        unmet_total += unmet
        cov_total += cov

        res[context] = conversion_rate_by_context[context]

    res["overall"] = unmet_total / cov_total

    return res


def estimate_conversion_rate_print(index_prefix, work_dir):
    conversion_rate = estimate_conversion_rate(index_prefix, work_dir)

    overall = conversion_rate["overall"] * 100
    cg_cr  = "Not available"
    chg_cr = "Not available"
    chh_cr = "Not available"

    if "CG" in conversion_rate:
        cg_cr = f"{conversion_rate['CG'] * 100:.2f}%"
    if "CHG" in conversion_rate:
        chg_cr = f"{conversion_rate['CHG'] * 100:.2f}%"
    if "CHH" in conversion_rate:
        chh_cr = f"{conversion_rate['CHH'] * 100:.2f}%"

    res = (f"Overall conversion rate: {overall}\n"
           f"CG context conversion rate: {cg_cr}\n"
           f"CHG context conversion rate: {chg_cr}\n"
           f"CHH context conversion rate: {chh_cr}\n")

    return res


def vg_binary_check(vg_path=None):
    fine = False
    if vg_path == None:
        vg_path = "vg"

        try:
            config = ConfigParser("config.ini")
            vg_path = config.get("default", "vg_path")
        except:
            pass

    cmd = f"{vg_path} version"

    vgp = SystemExecute()
    stdout, stderr = vgp.execute(cmd)
    vgp.wait()

    res = f"Checking vg binary @ {vg_path}\n"

    for l in stdout:
        res += l.decode('utf-8')

    for l in stderr:
        res += l.decode('utf-8')

    x = len(re.compile(r"vg version v\d*.\d*.\d*").findall(res))
    if x == 0:
        res += "Error: Cannot find vg version\n"
    else:
        res += "Success: vg binary seems to be fine\n"
        fine = True

    return fine, res


class HelpDocument(object):

    def __init__(self):
        pass

    def header(self):
        return """███╗   ███╗███████╗████████╗██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗  █████╗ ██████╗ ██╗  ██╗███████╗██████╗ 
████╗ ████║██╔════╝╚══██╔══╝██║  ██║╚██╗ ██╔╝██║     ██╔════╝ ██╔══██╗██╔══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗
██╔████╔██║█████╗     ██║   ███████║ ╚████╔╝ ██║     ██║  ███╗██████╔╝███████║██████╔╝███████║█████╗  ██████╔╝
██║╚██╔╝██║██╔══╝     ██║   ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔══██╗██╔══██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗
██║ ╚═╝ ██║███████╗   ██║   ██║  ██║   ██║   ███████╗╚██████╔╝██║  ██║██║  ██║██║     ██║  ██║███████╗██║  ██║
╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝
""".strip()

    def version(self):
        return "0.2.0"

    def help_text_raw(self):
        return """
Usage: python main.py <command> <arguments>
Commands:
    Help
    PrepareGenome
    Main
    Align
    MethylCall
    ConversionRate
    MergeCpG

Help:
    methylGrapher help
    Or
    Read the detailed documentation: https://twlab.github.io/methylGrapher/build/html/

PrepareGenome:
    It adds lambda phage genome to your genome graph, converts a GFA file into fully G->A and C->T converted GFA file, and indexes it for vg giraffe alignment.
    methylGrapher PrepareGenome 
    # Input options
    -gfa <gfa_file_path> 
    -lp <lambda_phage_genome_path>
    
    # Output options
    -prefix <output_prefix> 
    -compress <Y/N>
    
    # Computing options
    -t <number_of_thread(s)> 

Main:
    It automatically executes Align and MethylCall in a sequence.
    methylGrapher Main
    # Align options: refer to Align section
    # MethylCall options: refer to MethylCall section

Align:
    VG Giraffe alignment, please provide work directory and index prefix.
    methylGrapher Align 
    # Required parameters
    -index_prefix <prefix> 
    -fq1 <fastq_file_path> 
    -fq2 <fastq_file_path> 
    -work_dir <work_directory> 
    
    # Computing options
    -t <number_of_thread(s)> (default: 1)
    -directional <Y/N> (default: Y)
    -compress <Y/N> (default: N)

MethylCall:
    Methylation call from vg giraffe alignment result.
    methylGrapher MethylCall 
    -work_dir <work_directory>
    
    # Alignment filtering options
    -discard_multimapped <Y/N> (default: Y)
    -minimum_identity <minimum_identity> (default: 20)
    -minimum_mapq <minimum_mapq> (default: 0)
    
    # Methylation calling options
    -cg_only <Y/N> (default: Y), only output methylation call in CG context
    
    # Genotype calling options
    -genotyping_cytosine <Y/N> (default: N)
    
    # Computing options
    -t <number_of_thread(s)> (default: 1)
    -batch_size <batch_size> (default: 4096)

ConversionRate:
    Estimate conversion rate from methylation extraction. MethylCall must be executed before this step.
    methylGrapher ConversionRate 
    -index_prefix <prefix> 
    -work_dir <work_directory>

MergeCpG:
    Merge cytosine methylation call (graph.methyl) into CpG methylation call. During graph indexing, all CpG locations are identified and stored in a separate TSV file. The graph CpG locations are stored under {index_prefix}cpg.tsv, with CpG id and both cytosine location on graph coordinate. MergeCpG function will merge the cytosine methylation call (graph.methyl) into CpG methylation call using graph CpG id.
    methylGrapher MergeCpG 
    -index_prefix <prefix> 
    -work_dir <work_directory>


""".strip()

    def help_text(self):
        return f"\n\n{self.header()}\n\nDeveloped by Wenjin @ Wanglab (WUSTL)\nVersion: {self.version()}\n\n{self.help_text_raw()}\n\n\n\n\n\n\n\n"

    def abc(self):
        return """"""


if __name__ == "__main__":
    pass

    help = HelpDocument()
    print(help.help_text())



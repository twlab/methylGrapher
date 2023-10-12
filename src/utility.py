import os
import sys
import gzip
import time
# import ctypes
import random
import string
import hashlib
import subprocess
import multiprocessing

import pgzip




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





def gfa_converter(input_gfa, output_file_prefix, compress=True, thread=1):
    utils = Utility()

    for conversion_str in ['C2T', 'G2A']:
        conversion = conversion_str.split('2')

        output_gfa = f"{output_file_prefix}.{conversion_str}.gfa"

        output_file_handle = None
        if compress:
            output_gfa += '.gz'
            output_file_handle = pgzip.open(output_gfa, 'wt')
        else:
            output_file_handle = open(output_gfa, 'w')

        input_compress = False

        if utils.isGzip(input_gfa):
            input_compress = True

        # print(thread)
        if input_compress:
            input_file_handle = pgzip.open(input_gfa, 'rt', thread=thread)
        else:
            input_file_handle = open(input_gfa, 'r')

        for l in input_file_handle:
            l = l.strip().split('\t')

            if l[0] == 'S':
                l[2] = l[2].replace(conversion[0], conversion[1])
                l[2] = l[2].replace(conversion[0].lower(), conversion[1].lower())

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
            stdout = open(stdout, "w")

        if stderr is None:
            stderr = subprocess.PIPE
        else:
            if not os.path.exists(os.path.dirname(stderr)):
                os.makedirs(os.path.dirname(stderr))
            stderr = open(stderr, "w")

        sp = subprocess.Popen(cmd, shell=True, stdout=stdout, stderr=stderr)
        self._pool.append(sp)

        stdout, stderr = sp.stdout, sp.stderr
        return stdout, stderr

    def wait(self):
        for sp in self._pool:
            sp.wait()
        self._pool = []


def fastq_converter(fq1, fq2, workdir, compress=True, thread=1, directional=True):
    utils = Utility()

    report_file_handle = open(f"{workdir}/report.txt", 'a')
    # read_record_file_handle = pgzip.open(f"{workdir}/read_record.gz", 'wt', thread=thread)

    read_count = 0
    for fq in (fq1, fq2):
        read_count += 1

        for conversion_str in ['C2T', 'G2A']:

            if directional:
                if read_count == 1 and conversion_str == 'G2A':
                    continue
                if read_count == 2 and conversion_str == 'C2T':
                    continue

            conversion = conversion_str.split('2')


            input_fastq = fq
            if fq == None:
                continue
            output_fastq = f"{workdir}/{conversion_str}.R{read_count}.fastq"

            input_fastq_fh = None
            input_compress = utils.isGzip(input_fastq)
            if input_compress:
                input_fastq_fh = pgzip.open(input_fastq, 'rt', thread=thread)
            else:
                input_fastq_fh = open(input_fastq, 'r')


            output_file_handle = None
            if compress:
                output_fastq += '.gz'
                output_file_handle = pgzip.open(output_fastq, 'wt', thread=thread)
            else:
                output_file_handle = open(output_fastq, 'w')

            seq = ""
            for i, l in enumerate(input_fastq_fh):
                newl = l
                if i % 4 == 1:
                    # Fully converted
                    seq = l[:]
                    newl = newl.replace(conversion[0], conversion[1])
                    newl = newl.replace(conversion[0].lower(), conversion[1].lower())
                elif i % 4 == 3:
                    # Encode original sequence as PHRED score
                    newl = seq
                elif i % 4 == 0:
                    assert l.startswith('@')
                    original_header = l.strip()[1:]
                    original_header = original_header.split(' ')

                    x = str(int(i/4))[-1]
                    original_header[0] = f"@{original_header[0]}{conversion_str}{x}R{read_count}"
                    newl = ' '.join(original_header) + '\n'


                output_file_handle.write(newl)

                #if i % 4 in [0, 1] and conversion_str == 'C2T':
                #    read_record_file_handle.write(newl)
            i+=1

        report_file_handle.write(f"Total reads for R{read_count}: {int(i/4)}\n")

    report_file_handle.close()

    return 0


if __name__ == "__main__":
    pass



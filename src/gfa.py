
import os
import sys
import time
import sqlite3




class GFAError(RuntimeError):

    def __init__(self, msg=""):
        self._msg = msg
        return

    def __str__(self):
        return self._msg

class UnsupportedGFAError(GFAError):

    def __init__(self, msg=""):
        super().__init__(msg=msg)


class NotBluntEndLink(UnsupportedGFAError):

    def __init__(self, cigar):
        super().__init__(msg=f"{cigar} is not supported, only blunt ends are supported")


class GraphicalFragmentAssemblyAbstract(object):

    # It defines the interface for the GFA class.
    # Abstract class, it cannot be instantiated.

    def __init__(self):
        self.clear()

    def clear(self):
        # clear/initiate all internal data structures
        raise NotImplementedError

    def parse(self, gfa_file):
        raise NotImplementedError

    def write(self, gfa_file):
        raise NotImplementedError

    def get_sequence_by_segment_ID(self, segment_ID):
        raise NotImplementedError

    def get_tag_by_segment_ID(self, segment_ID):
        raise NotImplementedError

    def get_parent_links(self, segment_ID):
        raise NotImplementedError

    def get_child_links(self, segment_ID):
        raise NotImplementedError

    def get_walk(self, something):
        raise NotImplementedError

    def get_path(self, something):
        raise NotImplementedError

    def get_sequence_by_defined_path(self, something):
        raise NotImplementedError


    def all_segment_IDs(self):
        raise NotImplementedError

    def get_replacement_SNV(self, base1, base2):
        # Remove base2 SNV when base1 is present

        # replacement key: base1(removed) SNV
        replacement = {}
        replacement_snv = base1 + base2
        # print(base1, base2)


        for segID in self.all_segment_IDs():
            potential_children_link = []

            tmp = self.get_child_links(segID)
            if len(tmp) < 2:
                continue

            for cl in self.get_child_links(segID):
                cID, cl_strand, cl_overlap = cl

                if self.get_parent_link_count(cID) != 1:
                    continue

                if len(self.get_child_links(cID)) != 1:
                    continue

                c_seq = self.get_sequence_by_segment_ID(cID)
                # TODO
                # if cl_strand in ["+-", "-+"]:
                # c_seq = utils.reverse_complement(c_seq)
                if c_seq not in replacement_snv:
                    continue
                cl = [cID, cl_strand, cl_overlap, c_seq]
                potential_children_link.append(cl)

            if len(potential_children_link) != 2:
                continue

            candidate = {}
            #print(segID, potential_children_link)
            for pcl in potential_children_link:
                cID, cl_strand, cl_overlap, snv = pcl
                grandchildrenlink = self.get_child_links(cID)[0]
                if grandchildrenlink[1] != "++":
                    continue
                #print(grandchildrenlink)
                if snv not in candidate:
                    candidate[snv] = []
                candidate[snv].append((cID, grandchildrenlink[0]))

            for snv1 in candidate.get(base1, []):
                for snv2 in candidate.get(base2, []):
                    if snv1[1] == snv2[1]:
                        replacement[snv2[0]] = snv1[0]
                        #print(snv2[0], snv1[0])

            # print(candidate)

        return replacement








# Store entire GFA in memory
class GraphicalFragmentAssemblyMemory(GraphicalFragmentAssemblyAbstract):

    def __init__(self):
        super().__init__()

    def clear(self):
        self._header = ""
        self._segment = {}
        self._walk = {}

        self._original_gfa_path = ""

    def parse(self, gfa_file, keep_tag=False, keep_link=False):
        self._original_gfa_path = gfa_file
        # TODO when not keep_tag, keep_link, forbid to write GFA and access those two
        for l in open(gfa_file):
            if l[0] not in "HSL":
                continue

            l = l.strip().split("\t")

            if l[0] == "H":
                self._header = "\t".join(l[1:])

            if l[0] == "S":
                rt, segID, seq, *tags = l
                # sequence, tag, links, parent_links_count
                if not keep_tag:
                    tags = []
                self._segment[segID] = [seq, tags, [], 0]

            elif l[0] == "L":
                if not keep_link:
                    continue
                rt, pID, s1, cID, s2, cigar, *tags = l
                tags = "\t".join(tags)
                if cigar not in ["*", "0M"]:
                    raise NotBluntEndLink(cigar)
                self._segment[pID][2].append([cID, s1 + s2, cigar])
                self._segment[cID][3] += 1
            elif l[0] == "P":
                pass
            elif l[0] == "W":
                pass
            else:
                continue


    def get_sequence_by_segment_ID(self, segment_ID):
        return self._segment[segment_ID][0]

    def get_sequences_by_segment_ID(self, segment_IDs):
        res = {}
        for sID in segment_IDs:
            res[sID] = self.get_sequence_by_segment_ID(sID)
        return res

    def get_tag_by_segment_ID(self, segment_ID):
        raise NotImplementedError

    def get_parent_link_count(self, segment_ID):
        return self._segment[segment_ID][3]

    def get_parent_links(self, segment_ID):
        raise NotImplementedError

    def get_child_links(self, segment_ID):
        return self._segment[segment_ID][2]

    def all_segment_IDs(self):
        return self._segment.keys()

    def write_converted(self, output_file, replace_base1, replace_base2, SNV_trim=True):

        # Remove base2 SNV when base1 is present
        replacement = {}
        if SNV_trim:
            replacement = self.get_replacement_SNV(replace_base1, replace_base2)

        replaced_gfa_file = open(output_file, "w")
        with replaced_gfa_file as f:
            f.write("H\t" + self._header + "\n")
            for segID in self._segment:
                if segID in replacement:
                    continue
                seq, tags, links, parent_links_count = self._segment[segID]
                seq = seq.replace(replace_base1, replace_base2)
                tags = "\t".join(tags)
                f.write(f"S\t{segID}\t{seq}\t{tags}\n")
                for l in links:
                    cID, s1s2, cigar = l
                    if cID in replacement:
                        continue
                    s1 = s1s2[0]
                    s2 = s1s2[1]
                    f.write(f"L\t{segID}\t{s1}\t{cID}\t{s2}\t{cigar}\n")

            with open(self._original_gfa_path) as f:
                for l in f:
                    # TODO: support for P
                    if l[0] not in "W":
                        continue

                    if not SNV_trim:
                        replaced_gfa_file.write(l)
                        continue

                    # Graph Trimming
                    l = l.strip().split("\t")
                    # l[6] = l[6][:100]
                    walk = l[6]
                    walk_splited = []
                    node = ["", ""]
                    for c in walk:
                        if c in "><":
                            if node[0] != "" and node[1] != "":
                                if node[1] in replacement:
                                    node[1] = replacement[node[1]]
                                walk_splited.append(node[0]+node[1])
                            node = ["", ""]
                            node[0] = c
                            continue
                        node[1] += c
                    if node[1] in replacement:
                        node[1] = replacement[node[1]]
                    walk_splited.append(node[0]+node[1])

                    l[6] = "".join(walk_splited)

                    # replaced_gfa_file.write(l)
                    replaced_gfa_file.write("\t".join(l) + "\n")


class GraphicalFragmentAssemblyMemorySegmentOptimized(GraphicalFragmentAssemblyAbstract):

    def __init__(self):
        super().__init__()

    def clear(self):
        self._numeric_segment_ID = False
        self._segment = {}
        self._SN_segments = {
            "A": set(),
            "T": set(),
            "C": set(),
            "G": set(),
        }

        self._original_gfa_path = ""

    def parse(self, gfa_file, keep_tag=False, keep_link=False):
        if keep_link:
            raise NotImplementedError
        if keep_tag:
            raise NotImplementedError
        self._original_gfa_path = gfa_file

        """
        Well, it does take less memory, but it is slower in look up
        numeric_segment_ID = True
        for l in open(gfa_file):
            if not l.startswith("S"):
                continue
            l = l.strip().split("\t")
            rt, segID, seq, *tags = l

            try:
                int(segID)
            except ValueError:
                numeric_segment_ID = False
                break

        self._numeric_segment_ID = numeric_segment_ID
        """

        for l in open(gfa_file):
            if not l.startswith("S"):
                continue

            l = l.strip().split("\t")

            rt, segID, seq, *tags = l

            #if self._numeric_segment_ID:
            #    segID = int(segID)
            seq = seq.upper()

            if seq in self._SN_segments:
                self._SN_segments[seq].add(segID)
            else:
                self._segment[segID] = seq



    def get_sequence_by_segment_ID(self, segment_ID):
        for sn in "ATCG":
            if segment_ID in self._SN_segments[sn]:
                return sn
        return self._segment[segment_ID]

    def get_sequences_by_segment_ID(self, segment_IDs):
        res = {}
        for sID in segment_IDs:
            res[sID] = self.get_sequence_by_segment_ID(sID)
        return res

    def get_sequence_length_by_segment_ID(self, segment_ID):
        return len(self.get_sequence_by_segment_ID(segment_ID))


# Store segment length in memory
class GraphicalFragmentAssemblySegmentLengthMemory(object):

    def __init__(self):
        self.clear()

    def clear(self):
        self._segment_length = {}

    def parse(self, gfa_file):
        with open(gfa_file) as gfa_fh:
            for l in gfa_fh:
                if l[0] not in "S":
                    continue

                l = l.strip().split("\t")

                if l[0] == "S":
                    rt, segID, seq, *tags = l
                    # sequence, tag, links, parent_links_count
                    self._segment_length[segID] = len(seq)


    def get_sequence_length_by_segment_ID(self, segment_ID):
        return self._segment_length[segment_ID]




# Store GFA in sqlite3 database
class GraphicalFragmentAssemblySQL(GraphicalFragmentAssemblyAbstract):

    def __init__(self):
        super().__init__()

    def clear(self):
        self.connection = None
        self.cursor_obj = None

    # Connect to database
    def parse(self, gfa_db_path, read_only=False):

        if read_only:
            gfa_db_path = "file:" + gfa_db_path + "?mode=ro"
        self.connection = sqlite3.connect(gfa_db_path, uri=True)
        self.cursor = self.connection.cursor()
        return None

    def execute(self, query, *args):
        return self.cursor.execute(query, *args)

    def create_new_table(self):
        query = []

        query.append(
            'CREATE TABLE IF NOT EXISTS metadata ('
            'key TEXT, '
            'value TEXT '
            ')'
        )

        query.append('CREATE TABLE IF NOT EXISTS segment ('
            'name TEXT PRIMARY KEY, '
            'sequence TEXT, '
            'SN TEXT, '
            'SR INT, '
            'SO INT, '
            'tag TEXT'
            ')'
        )

        query.append(
            'CREATE TABLE IF NOT EXISTS link ('
            'parentID TEXT, '
            'childID TEXT, '
            'strand TEXT, '
            'overlap TEXT, '
            'tag TEXT'
            ')'
        )

        query.append(
            'CREATE TABLE IF NOT EXISTS path ('
            'name TEXT PRIMARY KEY, '
            'path TEXT, '
            'overlap TEXT'
            ')'
        )

        query.append(
            'CREATE TABLE IF NOT EXISTS walk ('
            'sample TEXT, '
            'haplotype INT, '
            'sequenceID TEXT, '
            'start INT, '
            'end INT, '
            'walk TEXT'
            ')'
        )

        query.append(
            "CREATE INDEX link_index ON link (parentID, childID)"
        )
        query.append(
            "CREATE INDEX parentID_index ON link (parentID)"
        )
        query.append(
            "CREATE INDEX childID_index ON link (childID)"
        )

        for q in query:
            self.cursor.execute(q)


    def gfa_to_db(self, gfa_file):
        i = 0
        for l in open(gfa_file):
            i += 1
            l = l.strip().split("\t")

            if l[0] == "H":
                pass
                #print(l)

            if l[0] == "S":
                rt, segID, seq, *tags = l
                sn = None
                sr = None
                so = None

                other_tags = []
                for t in tags:
                    if t.startswith("SN:Z:"):
                        sn = t[5:]
                    elif t.startswith("SR:i:"):
                        sr = int(t[5:])
                    elif t.startswith("SO:i:"):
                        so = int(t[5:])
                    else:
                        other_tags.append(t)

                tags = "\t".join(other_tags)
                if tags == "":
                    tags = None
                # print(segID, len(seq), tags)
                # print(l)
                # print()
                self.execute(
                    'INSERT INTO segment (name, sequence, tag, sn, sr, so) VALUES (?, ?, ?, ?, ?, ?)',
                    (segID, seq, tags, sn, sr, so)
                )
            elif l[0] == "L":
                rt, pID, s1, cID, s2, cigar, *tags = l
                tags = "\t".join(tags)
                if cigar in ["*", "0M"]:
                    cigar = None
                # print(pID, cID, s1, s2, cigar, tags)
                self.execute(
                    'INSERT INTO link (parentID, childID, strand, overlap, tag) VALUES (?, ?, ?, ?, ?)',
                    (pID, cID, s1 + s2, cigar, tags)
                )

            elif l[0] == "P":
                rt, name, path, cigar = l
                if cigar in ["*", "0M"]:
                    cigar = None

                # pos = l[2].count("+")
                # neg = l[2].count("-")
                # print("P", l[1], l[2][:10], l[3])

                self.execute(
                    'INSERT INTO path (name, path, overlap) VALUES (?, ?, ?)',
                    (name, path, cigar)
                )
            elif l[0] == "W":
                # pos = walk.count(">")
                # neg = walk.count("<")

                # print("W", len(l), l[1:6])
                # print("POS vs NEG", pos, neg)
                for j in [2, 4, 5]:
                    l[j] = int(l[j])

                self.execute(
                    'INSERT INTO walk (sample, haplotype, sequenceID, start, end, walk) VALUES (?, ?, ?, ?, ?, ?)',
                    (l[1:])
                )
            else:
                continue

        self.connection.commit()


    def get_sequence_by_segment_ID(self, segment_ID):
        return self.get_sequences_by_segment_ID([segment_ID])[segment_ID]


    def get_sequences_by_segment_ID_query_constructor(self, segment_IDs):
        segment_IDs = tuple(segment_IDs).__str__().replace(",)", ")")
        q = f"""
            SELECT name, sequence FROM segment
            WHERE name IN {segment_IDs}
        """
        return q
    def get_sequences_by_segment_ID(self, segment_IDs):
        q = self.get_sequences_by_segment_ID_query_constructor(segment_IDs)

        res = {}
        for i in self.execute(q):
            name, seq = i
            res[name] = seq
        return res

    def get_tag_by_segment_ID(self, segment_ID):
        raise NotImplementedError

    def get_parent_links(self, segment_ID):
        q = f"""
            SELECT parentID, strand, overlap FROM link
            WHERE childID == '{segment_ID}'
        """

        res = []
        for i in self.execute(q):
            res.append(tuple(i))
        return res

    def get_child_links(self, segment_ID):
        q = f"""
            SELECT childID, strand, overlap FROM link
            WHERE parentID == '{segment_ID}'
        """

        res = []
        for i in self.execute(q):
            res.append(tuple(i))
        return res

    def find_parallel_segment_on_ref(self, segment_ID):

        fail = False

        # Find children that are on reference
        q1 = f"""
            SELECT link.childID, link.strand, segment.SN, segment.SO, segment.sequence FROM link
            LEFT JOIN segment 
            ON segment.name == link.childID
            WHERE link.parentID == '{segment_ID}' and segment.SR == 0
        """

        child_res = None
        for i in self.execute(q1):
            if child_res is not None:
                fail = True
                break
            child_res = list(i)

        if fail or child_res is None:
            return None

        q2 = f"""
            SELECT link.parentID, link.strand, segment.SN, segment.SO, segment.sequence FROM link
            LEFT JOIN segment 
            ON segment.name == link.parentID
            WHERE link.childID == '{segment_ID}' and segment.SR == 0
        """

        parent_res = None
        for i in self.execute(q2):
            if parent_res is not None:
                fail = True
                break
            parent_res = list(i)

        if fail or parent_res is None:
            return None

        q3 = f"""
                SELECT link.childID, link.strand, segment.SN, segment.SO, segment.sequence FROM link
                LEFT JOIN segment 
                ON segment.name == link.childID
                WHERE (link.parentID == '{parent_res[0]}' AND segment.SR == 0)
            """

        carried_over_segment_candidate = []
        for i in self.execute(q3):
            carried_over_segment_candidate.append(list(i))

        q4 = f"""
                SELECT link.parentID, link.strand, segment.SN, segment.SO, segment.sequence FROM link
                LEFT JOIN segment 
                ON segment.name == link.parentID
                WHERE (link.childID == '{child_res[0]}' AND segment.SR == 0)
            """

        carried_over_segment = []
        for i in self.execute(q4):
            sid = list(i)[0]
            for cosc in carried_over_segment_candidate:
                if cosc[0] == sid:
                    carried_over_segment.append(cosc)


        print("LIFT OVER")
        print(parent_res, "\n", child_res)
        print(carried_over_segment)
        res = []
        return res



def add_lambda_genome_to_gfa(igfa, ogfa, lambda_ref):
    figfa = open(igfa)
    fogfa = open(ogfa, "w")

    if lambda_ref is None:
        for l in figfa:
            fogfa.write(l)
        return

    lseq = ""
    for l in open(lambda_ref):
        if l[0] == ">":
            continue
        lseq += l.strip()

    max_segid = 0
    for l in figfa:
        if l[0] == "S":
            rt, segID, seq, *tags = l.strip().split("\t")
            max_segid = max(max_segid, int(segID))
        fogfa.write(l)

    max_segid += 1
    fogfa.write(f"S\t{max_segid}\t{lseq}\n")
    fogfa.write(f"W	LambdaPhage	0	lambdaphagegenome	0	{len(lseq)}	>{max_segid}\n")

    return max_segid




if __name__ == "__main__":

    if len(sys.argv) == 4:
        add_lambda_genome_to_gfa(sys.argv[1], sys.argv[2], sys.argv[3])
        sys.exit(0)

    if len(sys.argv) == 5:
        assert sys.argv[1] == "convert"
        base1, base2 = sys.argv[2]
        gfa_path = sys.argv[3]
        out_gfa_path = sys.argv[4]

        gfa = GraphicalFragmentAssemblyMemory()
        gfa.parse(gfa_path, keep_link=True)
        gfa.write_converted(out_gfa_path, base1, base2)

    sys.exit(0)

    gfa_path1 = "../test_data/chr20-mcpg.gfa"
    gfa_db_path1 = "../test_data/chr20-mcpg.db"

    gfa_path2 = "../test_data/hprc.gfa"
    gfa_db_path2 = "../test_data/hprc.db"

    gfa_path = gfa_path1
    gfa_db_path = gfa_db_path1

    if True:
        gfa_min = GraphicalFragmentAssemblyMemory()
        gfa_min.parse(gfa_path, keep_link=True)

        a = gfa_min.get_sequence_by_segment_ID("554")
        b = gfa_min.get_sequences_by_segment_ID(["555", "556"])

        #c = gfa_min.get_parent_links("554")
        d = gfa_min.get_child_links("554")


        print(a)
        print(b)
        #print(c)
        print(d)

        # gfa_min.replace_SNV("C", "T")
        gfa_min.write_converted("x.gfa", "C", "T")

        print("Finish")
        time.sleep(100)
    else:
        gfa_db = GraphicalFragmentAssemblySQL()
        gfa_db.parse(gfa_db_path)
        # gfa_db.create_new_table()
        # gfa_db.gfa_to_db(gfa_path)



        gfa_db = GraphicalFragmentAssemblySQL()
        gfa_db.parse(gfa_db_path, read_only=True)


        a = gfa_db.get_sequence_by_segment_ID("554")
        b = gfa_db.get_sequences_by_segment_ID(["555", "556"])

        c = gfa_db.get_parent_links("554")
        d = gfa_db.get_child_links("554")

        p = gfa_db.get_sequence_by_defined_path(["554", "555", "556"])

        print(a)
        print(b)
        print(c)
        print(d)



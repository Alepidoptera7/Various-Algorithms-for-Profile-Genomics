import sys
import os
import subprocess
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO


class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
    print (head,seq)
    '''

    def __init__(self, fname):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            #header = ''
            #sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
                header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class tree_maker:
    """

    """

    def __init__(self):
        self.codon_dict = {
            'TCA': 'S',  # Serina
            'TCC': 'S',  # Serina
            'TCG': 'S',  # Serina
            'TCT': 'S',  # Serina
            'TTC': 'F',  # Fenilalanina
            'TTT': 'F',  # Fenilalanina
            'TTA': 'L',  # Leucina
            'TTG': 'L',  # Leucina
            'TAC': 'Y',  # Tirosina
            'TAT': 'Y',  # Tirosina
            'TAA': '*',  # Stop
            'TAG': '*',  # Stop
            'TGC': 'C',  # Cisteina
            'TGT': 'C',  # Cisteina
            'TGA': '*',  # Stop
            'TGG': 'W',  # Triptofano
            'CTA': 'L',  # Leucina
            'CTC': 'L',  # Leucina
            'CTG': 'L',  # Leucina
            'CTT': 'L',  # Leucina
            'CCA': 'P',  # Prolina
            'CCC': 'P',  # Prolina
            'CCG': 'P',  # Prolina
            'CCT': 'P',  # Prolina
            'CAC': 'H',  # Histidina
            'CAT': 'H',  # Histidina
            'CAA': 'Q',  # Glutamina
            'CAG': 'Q',  # Glutamina
            'CGA': 'R',  # Arginina
            'CGC': 'R',  # Arginina
            'CGG': 'R',  # Arginina
            'CGT': 'R',  # Arginina
            'ATA': 'I',  # Isoleucina
            'ATC': 'I',  # Isoleucina
            'ATT': 'I',  # Isoleucina
            'ATG': 'M',  # Methionina
            'ACA': 'T',  # Treonina
            'ACC': 'T',  # Treonina
            'ACG': 'T',  # Treonina
            'ACT': 'T',  # Treonina
            'AAC': 'N',  # Asparagina
            'AAT': 'N',  # Asparagina
            'AAA': 'K',  # Lisina
            'AAG': 'K',  # Lisina
            'AGC': 'S',  # Serina
            'AGT': 'S',  # Serina
            'AGA': 'R',  # Arginina
            'AGG': 'R',  # Arginina
            'GTA': 'V',  # Valina
            'GTC': 'V',  # Valina
            'GTG': 'V',  # Valina
            'GTT': 'V',  # Valina
            'GCA': 'A',  # Alanina
            'GCC': 'A',  # Alanina
            'GCG': 'A',  # Alanina
            'GCT': 'A',  # Alanina
            'GAC': 'D',  # Acido Aspartico
            'GAT': 'D',  # Acido Aspartico
            'GAA': 'E',  # Acido Glutamico
            'GAG': 'E',  # Acido Glutamico
            'GGA': 'G',  # Glicina
            'GGC': 'G',  # Glicina
            'GGG': 'G',  # Glicina
            'GGT': 'G'  # Glicina
        }

        #adapter sequences
        self.p5 = 'AATGATACGGCGACCACCGAGATCT'
        self.p7 = 'CAAGCAGAAGACGGCATACGAGAT'

        #for assembly
        self.reverse_complement_dict = {}

        #for baits
        self.baits_by_head_dict = {}

        self.seq_dict_by_head = {}
        self.head_dict_by_seq = {}

        self.duplicate_reads_by_file_set = set()
        self.all_reads_by_file_set = set()

        self.duplicate_sequence_dict = {}
        self.read_bait_correspondence_dict = {}

        self.unique_sequence_dict = {}
        self.duplicate_header_dict = {}

        self.fourmer_frequency_dict = {}
        self.fourmer_frequency_by_file_dict = {}
        self.fourmer_sum_dict = {}

        self.forward_correspondence_dict = {}

        self.sink_dict = {}
        self.start_list = []
        self.end_list = []

        self.path_dict = {}
        self.direct_path_dict = {}

        self.assembly_seq_dict = {}

    def file_converter(self):
        """Designed to open a folder on desktop containing .fq.gz data from sample reads.

         These .fq.gz file names used to develop .fa files via bash subprocesses.

         Input: files in path
         Output: converted .fa files for use
        """

        path = "C://Users//Quin The Conquoror!//Desktop//Bobs Files fq.gz/Bobs files .fa"
        dir_list = os.listdir(path)

        for infile_name in dir_list:
            print(infile_name)
            outfile_name = infile_name.split('.')[0] + '.fa'
            print(outfile_name)

            #run the seqtk bash command which converts .fq.gz files to .fa for use
            #the .fa files are placed in a folder on the desktop for later use.
            subprocess_command = "seqtk seq -a " + infile_name + " > " + outfile_name
            print(subprocess_command)
            subprocess.run([subprocess_command], cwd=path, shell=True)

    def protein_translator(self, seq):
        """"""
        protein = ''
        codon_list = [seq[i:i + 3] for i in range(0, len(seq), 3)]

        for codon in codon_list:
            if codon in self.codon_dict.keys():
                protein += self.codon_dict[codon]

        print(protein)

    def gc_content(self, seq):
        """Develops GC content by read and file."""

        return seq.count("C") + seq.count("G") / len(seq)

    def reverse_complement(self, seq):
        """Develops the reverse complement of a parameter sequence."""

        return seq.lower().replace('a', 'T').replace('t', 'A').replace('g', 'C').replace('c', 'G').replace(' ', '')[::-1]

    def paired_end_assembly(self):
        """Assembles Illumina Hiseq paired end reads.

        3 cases:
        1) R1 ends in R2
        2) R1 begins in R2
        3) No overlap
        """

        for read_header in self.seq_dict_by_head.keys():
            if "2:N:0" in read_header:
                #provides a dictionary for the R2 reads
                self.reverse_complement_dict[read_header] = self.reverse_complement(self.seq_dict_by_head[read_header])

        assembly_case_dict = {i: [] for i in range(0, 3)}

        #searching for overlap between paired sequences
        for read_header_r1, read_header_r2 in zip(self.seq_dict_by_head.keys(), self.reverse_complement_dict.keys()):
            if read_header_r1[:-15] == read_header_r2[:-15]:
                for i in range(0, len(self.seq_dict_by_head[read_header_r1]), 2):

                    if self.seq_dict_by_head[read_header_r1][:-i] == self.reverse_complement_dict[read_header_r2][i:]:
                        print("case 1 agreement: ", self.seq_dict_by_head[read_header_r1][:-i])
                        print("Read headers: ", read_header_r1, read_header_r2)
                        print("Original seqs: ", self.seq_dict_by_head[read_header_r1], self.reverse_complement_dict[read_header_r2])
                        print("Assuy7ujhjyufuihfkjhhoiyhjgdcghfsreqadembled Seq:", self.reverse_complement_dict[read_header_r2][:i] + self.seq_dict_by_head[read_header_r1])

                        case1_seq = self.reverse_complement_dict[read_header_r2][:i] + self.seq_dict_by_head[read_header_r1]

                        print("paired seq len: ", len(case1_seq))
                        #assembly_case_dict[1].append()
                        print(" ")

                    if self.seq_dict_by_head[read_header_r1][i:] == self.reverse_complement_dict[read_header_r2][:-i]:
                        print("case 2 agreement:", self.seq_dict_by_head[read_header_r1][i:])
                        print("Read headers: ", read_header_r1, read_header_r2)
                        print("Original seqs: ", self.seq_dict_by_head[read_header_r1], self.reverse_complement_dict[read_header_r2])
                        print("Assembled seq: ", self.seq_dict_by_head[read_header_r1][:i] + self.reverse_complement_dict[read_header_r2])

                        case2_seq = self.seq_dict_by_head[read_header_r1][:i] + self.reverse_complement_dict[read_header_r2]

                        print("paired seq len: ", len(case2_seq))
                        #assembly_case_dict[2].append()
                        print(" ")

                    if self.seq_dict_by_head[read_header_r1][:-i] != self.reverse_complement_dict[read_header_r2][i:] and self.seq_dict_by_head[read_header_r1][i:] != self.reverse_complement_dict[read_header_r2][:-i]:
                        #print(read_header_r1, read_header_r2)
                        #print
                        #print("case 3")
                        #assembly_case_dict[3].append(((read_header_r1, read_header_r2), (self.seq_dict_by_head[read_header_r1], self.reverse_complement_dict[read_header_r2])))
                        pass



    def unique_and_duplicate_sequence_detection(self):
        """Detects read sequences shared by more than one header,
         and uses this to find all unique read sequences."""

        #find sequences associated with multiple headers -- duplicate sequences
        for head in self.seq_dict_by_head.keys():
            duplicate_sequence_list = []
            for seq in self.head_dict_by_seq.keys():
                if self.seq_dict_by_head[head] == seq and self.head_dict_by_seq[seq] != head:
                    if head not in self.duplicate_sequence_dict.keys():
                        self.duplicate_sequence_dict[head] = [(head, self.head_dict_by_seq[seq], seq)]

                    else:
                        duplicate_sequence_list.append((head, self.head_dict_by_seq[seq], seq))

            self.duplicate_sequence_dict[head].append(duplicate_sequence_list)

        #find unique sequences
        for head in self.seq_dict_by_head.keys():
            if self.seq_dict_by_head[head] not in self.duplicate_sequence_dict.keys():
                self.unique_sequence_dict[head] = self.seq_dict_by_head[head]

        print(self.unique_sequence_dict)

    def read_bait_association(self):
        """ The purpose of this def is to associate RNA baits with experimental read sequences.
        This is accomplished by searching for subsets of bait sequences within the read sequences.
        Baits must be broken into kmers, and the sequences must be searched for these kmers.
        Associated sequences will be paired. Sequences hybridized with baits may be exons.

        Specify the maximum number of hits, and the minimum length of acceptable hybrid sequences.
        """

        minimum_seq_length = 10
        read_number = 0

        for read_head in self.seq_dict_by_head.keys():
            read_bait_association_list = []
            for bait_head in self.baits_by_head_dict.keys():

                for j in range(len(self.seq_dict_by_head[read_head])):

                    if self.seq_dict_by_head[read_head][j:] in self.baits_by_head_dict[bait_head]:
                        if len(self.seq_dict_by_head[read_head][j:]) > minimum_seq_length:
                            read_bait_association_list.append((bait_head, self.seq_dict_by_head[read_head][j:], self.baits_by_head_dict[bait_head]))

                    if self.seq_dict_by_head[read_head][:-j] in self.baits_by_head_dict[bait_head]:
                        if len(self.seq_dict_by_head[read_head][:-j]) > minimum_seq_length:
                            read_bait_association_list.append((bait_head, self.seq_dict_by_head[read_head][:-j], self.baits_by_head_dict[bait_head]))

                    if j % 2 == 0:
                        if self.seq_dict_by_head[read_head][j:-j] in self.baits_by_head_dict[bait_head]:
                            if len(self.seq_dict_by_head[read_head][j:-j]) > minimum_seq_length:
                                read_bait_association_list.append((bait_head, self.seq_dict_by_head[read_head][j:-j], self.baits_by_head_dict[bait_head]))

                    if j % 3 == 0:
                        if self.seq_dict_by_head[read_head][j:-(j - 1)] in self.baits_by_head_dict[bait_head]:
                            if len(self.seq_dict_by_head[read_head][j:-(j - 1)]) > minimum_seq_length:
                                read_bait_association_list.append((bait_head, self.seq_dict_by_head[read_head][j:-(j - 1)], self.baits_by_head_dict[bait_head]))

            self.read_bait_correspondence_dict[read_head] = read_bait_association_list

            read_number += 1
            print("READ ", read_number, " : ", read_head, ", correspondences: ", len(self.read_bait_correspondence_dict[read_head]))

    def fourmer_frequency(self, file_name, head, seq, seq_number):
        """Splits all read sequences into strings of length 4.
        Matching 4-mers are counted and divided by 256 to develop a frequency value.
        This value is used as a specificity threshold."""

        fourmer_list = [seq[i:i + 4] for i in range(0, len(seq), 4)]

        if head == '':
            head = 'first read'

        self.fourmer_frequency_dict[head] = {fourmer: fourmer_list.count(fourmer)/256 for fourmer in fourmer_list}

        fourmer_sum_list = [fourmer_list.count(fourmer)/256 for fourmer in fourmer_list]
        self.fourmer_sum_dict[head] = sum(fourmer_sum_list)

        #print("READ ", seq_number,
        #     " : ", head, ", kmer Pr Values: ", self.fourmer_frequency_dict[head])

        for fourmer in fourmer_list:

            if fourmer not in self.fourmer_frequency_by_file_dict.keys():
                self.fourmer_frequency_by_file_dict[file_name] = {fourmer: fourmer_list.count(fourmer)/256 for fourmer in fourmer_list}

            else:
                self.fourmer_frequency_by_file_dict[file_name] = {fourmer: (self.fourmer_frequency_by_file_dict[file_name][fourmer] + (fourmer_list.count(fourmer) / 256))/seq_number
                                                                  for fourmer in fourmer_list}

    def sequence_correspondence_per_taxa(self):
        """Develops a forwards correspondence graph used in acyclic graphing alignment per taxa."""
        seq_num = 0
        sequence = ''

        for read_head in self.seq_dict_by_head.keys():
            forward_correspondence_list = []
            backwards_correspondence_list = []
            seq_num += 1
            self.forward_correspondence_dict[read_head] = []
            self.sink_dict[read_head] = []
            for other_read_head in self.seq_dict_by_head.keys():
                if read_head != other_read_head:
                    for j in range(8, len(self.seq_dict_by_head[read_head]), 8):
                        if self.seq_dict_by_head[read_head][j:] == self.seq_dict_by_head[other_read_head][:len(self.seq_dict_by_head[read_head][j:])]:

                            forward_correspondence_list.append(other_read_head)

                        if self.seq_dict_by_head[read_head][:j] == self.seq_dict_by_head[other_read_head][len(self.seq_dict_by_head[read_head][j:]):]:
                            backwards_correspondence_list.append(other_read_head)

            #print("in paths:", len(backwards_correspondence_list), ", out paths:", len(forward_correspondence_list))

            self.forward_correspondence_dict[read_head] = forward_correspondence_list

            self.sink_dict[read_head].append((backwards_correspondence_list, forward_correspondence_list))
        print(sequence)
        #print("READ", read_head, " alignments: ", self.alignment_correspondence_graph[read_head])

    def sink_finder(self):
        """Finds alignment origins by comparing number of inward edges to number of outward edges."""

        # finding sinks
        for read_head in self.sink_dict.keys():

            if len(self.sink_dict[read_head][0][1]) > len(self.sink_dict[read_head][0][0]):
                self.start_list.append((read_head, self.forward_correspondence_dict[read_head]))

            if len(self.sink_dict[read_head][0][1]) < len(self.sink_dict[read_head][0][0]):
                self.end_list.append((read_head, self.forward_correspondence_dict[read_head]))

        for possible_start in self.start_list:
            # hypothetical path made for each possible start node
            self.path_dict[possible_start[0]] = [possible_start[1]]

    def breadth_first_search(self):
        """Performs a breadth first search to find all possible alignments."""

        self.direct_path_dict = {start: [] for start in self.path_dict.keys()}

        #extracts the next node from the current
        for start in self.path_dict.keys():
            path_list = [(start, self.fourmer_sum_dict[start], start)]
            direct_path_list = []

            while self.forward_correspondence_dict[path_list[-1][0]]:

                for current_node in self.forward_correspondence_dict[path_list[-1][0]]:
                    if self.forward_correspondence_dict[current_node]:
                        direct_path_list.append((path_list[-1][0], current_node))

                    else:
                        self.direct_path_dict[start].append((len(direct_path_list), direct_path_list))
                        direct_path_list = []

                    path_list.append((current_node, self.fourmer_sum_dict[current_node], path_list[-1][0]))

            self.path_dict[start] = path_list
            #print("start: ", start)
        #print(self.direct_path_dict)

    def direct_path_assembly(self):
        """finds the optimal path from the end to the beginning by accessing fourmer frequency.
        The longest path with the highest specificity is chosen for assembly.
        """

        for start in self.direct_path_dict.keys():
            self.assembly_seq_dict[start] = []
            for path in self.direct_path_dict[start]:
                assembly_seq_list = []
                if path[0] > 0:
                    for node in path[1]:
                        assembly_seq_list.append((node[0], self.seq_dict_by_head[node[0]], self.fourmer_sum_dict[node[0]]))

                else:
                    break

                if assembly_seq_list:
                    self.assembly_seq_dict[start].append(assembly_seq_list)

                #print(start)
                #print(self.assembly_seq_dict[start])
                #print(" ")

    def alignment_per_taxa(self):
        """Returns optimally aligned sequence.
        Check if the current sequence in the current path fits to the next sequence in the path."""

        for start in self.assembly_seq_dict.keys():
            path_sum = 0

            for path in self.assembly_seq_dict[start]:
                sequence = ''

                for i in range(len(path)-1):
                    for j in range(len(path[i][1])):
                        if path[i][1][j:] == path[i+1][1][:len(path[i+1][1])-j]:
                            sequence += path[i][1][:j]
                            path_sum += path[i][2]
                            break

                seq_list = [tup[1] for tup in path]
                dup_count_list = [seq_list.count(seq) for seq in seq_list]

                sequence += path[-1][1]
                path_sum += path[-1][2]

                print("full path", path)
                print("sequence: ", sequence)
                print("path length:", len(path), "path sum:", path_sum)
                print("duplicate sequences in path", dup_count_list)
                #print("NCBI BLAST results:", self.blast(sequence))
                print(" ")

    def bait_parse(self):
        """Handles bait file."""

        bait_file_name = 'pagonaviticepts.fa'
        bait_seqs = FastAreader(bait_file_name)

        for head, seq in bait_seqs.readFasta():
            self.baits_by_head_dict[head] = seq

    def file_parse(self):
        """In parsing the fasta files, the bait sequences must first be assessed.

        Afterwards, the read files are unpacked sequentially.
        """

        path = "C://Users//Quin The Conquoror!//Desktop//Bobs Files fq.gz/Bobs files .fa"
        dir_list = os.listdir(path)

        for read_file_name in dir_list:
            read_file_data = FastAreader(read_file_name)
            print("file name: ", read_file_name)
            seq_number = 0
            self.fourmer_frequency_by_file_dict[read_file_name] = {}

            for head, seq in read_file_data.readFasta():
                #self.fourmer_frequency(read_file_name, head, seq, seq_number)
                self.seq_dict_by_head[head] = seq
                self.head_dict_by_seq[seq] = head
                seq_number += 1
                #if seq_number == 2500:
                #    break

    def blast(self, cluster_sequence):
        """Develops a BLAST request for each clustered sequence via BioPython.

        Design taken from Biopython Manual.
        """
        E_VALUE_THRESH = 0.0000000000000000001

        result_handle = NCBIWWW.qblast("blastn", "nt", sequence=cluster_sequence, expect=1)

        blast_record = NCBIXML.read(result_handle)

        count = 0

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    count += 1
                    print(f"sequence: {alignment.title}, length: {alignment.length}nt, e-value: {hsp.expect}")
                    #print(hsp.query[0:75] + "...")
                    #print(hsp.match[0:75] + "...")
                    #print(hsp.sbjct[0:75] + "...")

        print(f"There are {count} similar sequences in Blast output\n")

    def print_out(self):
        """Prints"""

        print("Number of reads: ", len(self.seq_dict_by_head))
        print("Number of baits: ", len(self.baits_by_head_dict))
        #print("normalized fourmer frequency Pr values: ", self.fourmer_frequency_by_file_dict)

    def driver(self):

        # self.file_converter()

        self.file_parse()
        self.paired_end_assembly()

        #self.print_out()

        #self.unique_and_duplicate_sequence_detection(head, seq)
        #self.protein_translator(seq)

        #when and how are baits used? Where are the exons?
        #self.read_bait_association()

        #self.sequence_correspondence_per_taxa()
        #self.sink_finder()
        #self.breadth_first_search()
        #self.direct_path_assembly()
        #self.alignment_per_taxa()


def main():
    """Access class and call driver member function."""

    class_access = tree_maker()
    class_access.driver()


if __name__ == '__main__':
    main()
#!/usr/bin/env python
"""
bam2fastq
"""

import shlex
import argparse
import subprocess
import time
import errno


def bitfield(num):
    """
    # returns bit array as little endian.  I add 4096 so that it will always
    # have all of the bit in sam flag field.
    """
    return [digit == '1' for digit in bin(int(num)+4096)][::-1]


def reverse_complement(seq):
    """
    reverse complements DNA sequence
    """
    flip = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return ''.join([flip[x] for x in seq[::-1]])


class SamLine(object):
    """
    SamLine class: for storing individual sam records
    """
    def __init__(self, read, read_group=None):
        self.read = read
        fields = read.split()
        self.name = fields[0]
        if read_group is None:
            for attrib in fields[11:]:
                if attrib[0:2] == "RG":
                    self.read_group = attrib.split(':')[-1]
                    break
            # self.read_group = read[0:5] + '.' + read[16]
        else:
            self.read_group = read_group
        self.flag = bin(int(fields[1]))
        self.unmapped = self.flag[-3] == "1"
        self.next_unmapped = self.flag[-4] == "1"
        self.reverse = self.flag[-5] == "1"
        self.secondary = self.flag[-9] == "1"
        self.seq = fields[9]
        self.qual = fields[10]

    def __str__(self):
        return self.read

    @property
    def sam(self):
        """
        returns the sam line in a FASTQ format
        """

        flag = 128 * (self.flag[-8] == "1") + 64 * (self.flag[-7] == "1") + 13
        if self.reverse:
            seq = reverse_complement(self.seq)
            qual = self.qual[::-1]
        else:
            seq = self.seq
            qual = self.qual
        attrib = "RG:Z:{rg}".format(rg=self.read_group)
        out = [self.name, str(flag), "*", "0",
               "0", "*", "*", "0", "0",
               seq, qual, attrib]
        return "\t".join(out)

    @property
    def fastq(self):
        """
        returns the sam line in a FASTQ format
        """
        fields = self.read.split("\t")
        seq = fields[9]
        qual = fields[10]
        if self.reverse:
            return "@%s\n%s\n+\n%s\n" % (self.name, reverse_complement(seq),
                                         qual[::-1])
        else:
            return "@%s\n%s\n+\n%s\n" % (self.name, seq, qual)

    def __lt__(self, other):
        return self.name < other.name

    def __eq__(self, other):
        return self.name == other.name


class Bam(object):
    """
    Class for the reading of the bam file
    """
    def __init__(self, bam):
        self.bam = bam
        self.header = None
        self.contigs = self._get_contigs()
        self.read_groups = self._get_read_groups()
        self.rgs = self._read_group_trim()
        self.files = None

    def _read_group_trim(self):
        """
        create the trimmed version of the readgroups
        """
        rgs = []
        for read_group in self.read_groups:
            fields = read_group.split(" ")
            for field in fields:
                if field[0:2] == "ID":
                    rgs.append(field[4:])
                    break
        return rgs

    def _get_header(self):
        """
        Get the header lines from the bam file
        """
        # a list of lists. First is by row and second is by column
        self.header = []
        command_line = shlex.split("samtools view -H")
        command_line.append(self.bam)
        samtools = subprocess.Popen(command_line,
                                    stdout=subprocess.PIPE)
        for line in samtools.stdout:
            header_fields = line.split("\t")
            self.header.append(header_fields)

    def _get_contigs(self):
        """
        Get the contigs present in the bam file
        """
        contigs = []
        if self.header is None:
            self._get_header()
        for header_fields in self.header:
            if header_fields[0] == "@SQ":
                contigs.append(header_fields[1].split(":")[1])
        return contigs

    def _get_read_groups(self):
        """
        Get the read groups from the incoming bam
        """
        read_groups = []
        if self.header is None:
            self._get_header()
        for header_fields in self.header:
            if header_fields[0] == "@RG":
                read_groups.append(' '.join(header_fields[1:]))
        return read_groups

    def get_read_pairs(self, contig=None):
        """
        process the bam to find the paired reads and output them to a FASTQ
        """
        # The magic 256 number is to filter out secondary alignments
        command_line = shlex.split("samtools view -F 256")
        command_line.append(self.bam)
        if contig is not None:
            command_line.append(contig)
        samtools = subprocess.Popen(command_line,
                                    stdout=subprocess.PIPE)
        reads = {}
        unmapped_reads = []
        for read in samtools.stdout:
            line = SamLine(read)
            pair = reads.pop(line.name, None)
            if pair is None:
                reads[line.name] = line
            else:
                for i in range(1000):
                    try:
                        print pair.sam
                        print line.sam
                    except IOError as e:
                        if e.errno == errno.EPIPE:
                            time.sleep(i)
                            continue
                        else:
                            raise
                    break
        unmapped_reads.sort()
        last_line = line
        for line in unmapped_reads:
            if last_line == line:
                print last_line.sam
                print line.sam


def main():
    """
    Main including the argument parser
    """
    parser = argparse.ArgumentParser(
        description="convert bam to unordered SAM for mapping")
    parser.add_argument('bam')
    args = parser.parse_args()
    in_bam = Bam(args.bam)

    in_bam.get_read_pairs()


if __name__ == '__main__':
    main()

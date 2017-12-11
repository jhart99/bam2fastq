#!/usr/bin/env python
"""
bam2fastq
"""

import time
import shlex
import argparse
import subprocess
import gzip
import os


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
        fields = read.split("\t")
        self.name = fields[0]
        if read_group is None:
            for attrib in fields[11:]:
                if attrib[0:2] == "RG":
                    self.read_group = attrib[6:]
                    break
            # self.read_group = read[0:5] + '.' + read[16]
        else:
            self.read_group = read_group
        self.flag = bin(int(fields[1]))
        self.unmapped = self.flag[-3] == "1"
        self.next_unmapped = self.flag[-4] == "1"
        self.reverse = self.flag[-5] == "1"
        self.secondary = self.flag[-9] == "1"
        # self.seq = fields[9]
        # self.qual = fields[10]

    def __str__(self):
        return self.read

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


class OutFiles(object):
    """
    A class to store the files which are being written
    """
    def __init__(self):
        self.files = {}
        self.path = "."

    def write_fastq(self, read1, read2):
        """
        writes a fastq line to the appropriate target for the read_group
        """
        target_files = self.files.get(read1.read_group)
        if target_files is None:
            target_files = (
                open(self.path + "/" + read1.read_group + ".1.fq", 'w'),
                open(self.path + "/" + read2.read_group + ".2.fq", 'w'))
            self.files[read1.read_group] = target_files
        target_files[0].write(read1.fastq)
        target_files[1].write(read2.fastq)

    def close(self):
        """
        close all of the files at the end
        """
        for (file1, file2) in self.files.values():
            file1.close()
            file2.close()


class OutFilesGzip(OutFiles):
    """
    A collection of files to write with gzip compression
    """
    def write_fastq(self, read1, read2):
        target_files = self.files.get(read1.read_group)
        if target_files is None:
            target_files = (
                gzip.open(
                    self.path + "/" + read1.read_group + ".1.fq.gz", 'w'),
                gzip.open(
                    self.path + "/" + read2.read_group + ".2.fq.gz", 'w'))
            self.files[read1.read_group] = target_files
        target_files[0].write(read1.fastq)
        target_files[1].write(read2.fastq)


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
        command_line = shlex.split("samtools view")
        command_line.append(self.bam)
        if contig is not None:
            command_line.append(contig)
        samtools = subprocess.Popen(command_line,
                                    stdout=subprocess.PIPE)
        reads = {}
        begin_time = time.clock()
        lines_read = 0
        lines_written = 0
        secondary_alignments = 0
        unmapped_reads = []
        for read in samtools.stdout:
            lines_read += 1
            if len(self.read_groups) > 1:
                line = SamLine(read)
            else:
                line = SamLine(read, self.rgs[0])
            if line.secondary:
                # don't want to include secondary alignments when
                # recreating the bam
                secondary_alignments += 1
                continue
            pair = reads.pop(line.name, None)
            if pair is None:
                reads[line.name] = line
            else:
                lines_written += 1
                self.files.write_fastq(pair, line)
            if lines_read % 100000 == 0:
                print lines_read, 2*lines_written, time.clock()-begin_time,
                print len(reads), secondary_alignments
        unmapped_reads.sort()
        last_line = line
        for line in unmapped_reads:
            if last_line == line:
                lines_written += 1
                self.files.write_fastq(last_line, line)
        self.files.close()


def main():
    """
    Main including the argument parser
    """
    parser = argparse.ArgumentParser(
        description="convert bam to multiple FASTQ")
    parser.add_argument('bam')
    parser.add_argument('outDir', default=".")
    parser.add_argument('--gzip', dest='fileType', action='store_const',
                        const=OutFilesGzip, default=OutFiles,
                        help='Output gzipped fastq')
    args = parser.parse_args()
    in_bam = Bam(args.bam)
    in_bam.files = args.fileType()
    in_bam.files.path = args.outDir
    read_group_file = "".join(os.path.basename(args.bam) + ".rg.txt")
    read_group_file = os.path.join(args.outDir, read_group_file)
    read_group_file_handle = open(read_group_file, 'w')
    for read_group in in_bam.read_groups:
        read_group_file_handle.write(read_group)
    read_group_file_handle.close()

    in_bam.get_read_pairs()


if __name__ == '__main__':
    main()

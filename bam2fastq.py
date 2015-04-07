import cProfile
import time
import fileinput
import shlex
import argparse
import subprocess

def bitfield(n):
    # returns bit array as little endian.  I add 4096 so that it will always 
    # have all of the bit in sam flag field.
    return [digit=='1' for digit in bin(int(n)+4096)][::-1]
#    return [1 if digit=='1' else 0 for digit in bin(int(n)+4096)][::-1]

def reverseComplement(seq):
    flip = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    return ''.join([flip[x] for x in seq[::-1]])

class samLine(object):
    def __init__(self, read):
        self.read = read
        fields = read.split("\t")
        self.name = fields[0]
        self.rg = read[0:5] + '.' + read[16]
        #self.flag = bitfield(fields[1])
        self.flag = bin(int(fields[1]))
        self.unmapped = self.flag[-3] == "1"
        self.nextUnmapped = self.flag[-4] == "1"
        self.reverse = self.flag[-5] == "1"
        self.secondary = self.flag[-9] == "1"
        #self.seq = fields[9]
        #self.qual = fields[10]
    def __str__(self):
        return self.read
    @property
    def fastq(self):
        fields = self.read.split("\t")
        seq = fields[9]
        qual = fields[10]
        if self.reverse:
            return "@%s\n%s\n+\n%s" % (self.name, reverseComplement(seq),
                    qual[::-1])
        else:
            return "@%s\n%s\n+\n%s" % (self.name, seq, qual)
    def __lt__(self, other):
        return self.name < other.name
    def __eq__(self, other):
        return self.name == other.name

class outFiles(object):
    def __init__(self):
        self.files = {}
    def writeFastq(self, read1, read2):
        _outFiles = self.files.get(read1.rg)
        if _outFiles is None:
            _outFiles = (open(read1.rg +".1.fq", 'w'),
                    open(read2.rg +".2.fq", 'w'))
            self.files[read1.rg] = _outFiles
        _outFiles[0].write(read1.fastq)
        _outFiles[1].write(read2.fastq)
    def close(self):
        for (file1, file2) in self.files.values():
            file1.close()
            file2.close()

class bam(object):
    def __init__(self, inBam):
        self.bam = inBam
        self.contigs = self._getContigs()
    def _getContigs(self):
        contigs = []
        commandLine = shlex.split("samtools view -H")
        commandLine.append(self.bam)
        samtools = subprocess.Popen(commandLine,
                stdout=subprocess.PIPE)
        for line in samtools.stdout:
            headerFields = line.split("\t")
            if headerFields[0] == "@SQ":
                contigs.append(headerFields[1].split(":")[1])
        return contigs
    def getReadPairs(self, contig=None):
        commandLine = shlex.split("samtools view")
        commandLine.append(self.bam)
        if contig is not None:
            commandLine.append(contig)
        samtools = subprocess.Popen(commandLine,
                stdout=subprocess.PIPE)
        reads = {}
        beginTime = time.clock()
        linesRead = 0
        linesWritten = 0
        unmappedReads = []
        for read in samtools.stdout:
            linesRead +=1
            line = samLine(read)
            pair = reads.pop(line.name, None)
            if pair is None:
                reads[line.name] = line
            else:
                linesWritten += 1
                self.files.writeFastq(pair, line)
            if linesRead % 100000 == 0:
                print linesRead, linesWritten, time.clock()-beginTime, len(reads)
            if linesRead % 2000000 == 0:
                break
        unmappedReads.sort()
        lastLine = line
        for line in unmappedReads:
            if lastLine == line:
                linesWritten += 1
                self.files.writeFastq(lastLine, line)
        self.files.close()


def main():
    parser = argparse.ArgumentParser(description="convert bam to multiple FASTQ")
    parser.add_argument('bam')
    args = parser.parse_args()

    


    inBam = bam(args.bam)
    inBam.files = outFiles()

    inBam.getReadPairs()
    

if __name__ == '__main__':
    main()
    #cProfile.run('main()')

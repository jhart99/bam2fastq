import cProfile
import time
import fileinput
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

def main():
    files = {}
    reads = {}
    beginTime = time.clock()
    linesRead = 0
    linesWritten = 0
    unmappedReads = []
    for read in fileinput.input():
        linesRead +=1
        line = samLine(read)
#        if line.unmapped or line.nextUnmapped:
#            unmappedReads.append(line)
        pair = reads.pop(line.name, None)
        if pair is None:
            reads[line.name] = line
        else:
            linesWritten += 1
            outFiles = files.get(line.rg)
            if outFiles is None:
                outFiles = (open(line.rg +".1.fq", 'w'),
                            open(line.rg +".2.fq", 'w'))
                files[line.rg] = outFiles
            outFiles[0].write(pair.fastq)
            outFiles[1].write(line.fastq)
        if linesRead % 100000 == 0:
            print linesRead, linesWritten, time.clock()-beginTime, len(reads), len(unmappedReads)
#        if linesRead % 2000000 == 0:
#            break
    unmappedReads.sort()
    lastLine = line
    for line in unmappedReads:
        if lastLine == line:
            linesWritten += 1
            outFiles = files.get(line.rg)
            if outFiles is None:
                outFiles = (open(line.rg +".1.fq", 'w'),
                            open(line.rg +".2.fq", 'w'))
                files[line.rg] = outFiles
            outFiles[0].write(lastLine.fastq)
            outFiles[1].write(line.fastq)
    for (file1, file2) in files.values():
        file1.close()
        file2.close()


#        print read[0:5] + '.' + read[16]
#        print read.split("\t", 1)[0]
#
#        print readread[0:5] + '.' + read[16]
#


if __name__ == '__main__':
    main()
    #cProfile.run('main()')

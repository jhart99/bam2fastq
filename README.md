# bam2fastq.py
## Description
This small utility takes in a bam file and creates a pair of FASTQ files per
read group.
## Motivation
I created this utility specifically to allow me to realign
sequencing data from the [Cancer Cell Line
Encyclopedia](https://cghub.ucsc.edu/datasets/ccle.html). For that data, only
bams are available and I needed to align them to a different build of the human
genome to facilitate certain analyses. I tried several different approaches
including Picard SamToFastq, Tophat's bam2fastx and Samtools -> unix shell
tools. Each had problems which I found discouraging. For each of these
approaches I found that producing the fastq files took longer than STAR to
align them.
## Use
python bam2fastq.py [bam]
## Results
2 Fastq files per read group named [RG].1.fq and [RG].2.fq
## Requirements
This utility requires samtools to be installed and in the path. No unusual
python modules are used. Memory requirements are small. In practice a maximum
of 50 mb was seen.  The program can take advantage of 2 cores, but should work
with 1 core. Samtools and this program run concurrently, but in practice the
python uses a full core while samtools uses 30-50% of a core. In my hands, on
my desktop which uses an Intel i5-3470 processor I see about 1 million reads
processed per 10 sec.

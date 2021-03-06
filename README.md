# bam2fastq
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

# New: bam2sam.py
## Description
This little script allows you to generate unmapped SAM files from
mapped BAM files.  This will let you directly remap mapped BAM files
to a new reference without a FASTQ intermediate.  To use this with
STAR, add --readFileType SAM PE or SAM SE as appropriate and
"--readFilesCommand bam2sam.py" or "--readFilesCommand bam2sam.py
--se".  This will remove all of the mapping information, flip the
reads if they are inverted when aligned, and strip out any mapping
flags other than read groups.  In practice, I've seen mapping speeds
comparable with STAR on FASTQ files.  Using this for on the fly
remapping will consume 1.5 cores and up to 200MB of memory in addition
to STAR's requirements.
[![Analytics](https://ga-beacon.appspot.com/UA-110461825-1/bam2fastq?pixel)](https://github.com/jhart99/bam2fastq)

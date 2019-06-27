# pileup consensus
Script that takes samtools pileup output and write consensus fasta file as output.

usage: pileup_consensus.py [-h] [-v] [-n BASE] [-r READS] [input] [output]

Produce consensus fatsa from samtools pileup output

positional arguments:
  input                 pileup file 
  output                Output file if no file result will be directed to
                        stander output

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -n BASE, --base BASE  Number of bases in line. (default: 80)
  -r READS, --reads READS
                        Minimum number of reads to accept nucleotide as a
                        variant. (default: 2)

pileup_consensus.py version 0.01. use command -h for info.

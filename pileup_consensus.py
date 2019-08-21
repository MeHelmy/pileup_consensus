
#!/usr/bin/env python3

"""
This script takes samtools pileup output and write out fasta file.
"""
from __future__ import print_function
import argparse
import sys, re



def get_args():
    parser = argparse.ArgumentParser(epilog="%(prog)s version 0.01. use command -h for info.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Produce consensus fatsa from samtools pileup output',
                                     add_help=True, )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')

    parser.add_argument('input', nargs='?', help="pileup file",
                             type=argparse.FileType('r'),
                             default=sys.stdin)
    parser.add_argument('output', nargs='?', help="Output file if no file result will be directed to stander output",
                                 type=argparse.FileType('w+'),
                                 default=sys.stdout)
    parser.add_argument('-n', '--base', help="Number of bases in line.", type=int, action='store', default=80)
    parser.add_argument('-r', '--reads', help="Minimum number of reads to accept nucleotide as a  variant.", type=int, action='store', default=2)


    parser.set_defaults(func=allel_count)

    # if not argument print help.
    if len(sys.argv) == 1 and sys.stdin.isatty():
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def allel_count(args):
    output = args.output
    input = args.input
    cov = {}
    first_time = True
    with output as data_out, input as data_in:
        # Info about samtools mpileup format https://en.wikipedia.org/wiki/Pileup_format
        chr = ""
        for line in data_in:
            contig, chrpos, ref, coverage, bases, *_ = line.split()
            # chr = chr if chr == contig else contig
            ref_upper = ref.upper()
            ref_lower = ref.lower()
            allowed_chars = set("ACGT")
            good_ref =  set(ref_upper).issubset(allowed_chars)
            cov['A'] = cov['C'] = cov['G'] = cov['T'] = cov['D'] = 0
            cov['a'] = cov['c'] = cov['g'] = cov['t'] = 0
            slen = len(bases)
            pos = 0
            while( pos < slen ):
                base = bases[pos]
                pos += 1
                if base == '^':
                    pos += 1
                elif base == '.' and good_ref:
                    cov[ref_upper] += 1
                elif base == ',' and good_ref:
                    cov[ref_lower] += 1
                elif base in 'ACGTNacgtn':
                    cov[base] += 1
                elif base == '+' or base == '-':
                    indel = re.search('[0-9]+[AGCTNagctn]+', bases[pos:]).group()
                    # adding how many deletion ot insertion we have
                    cov['D'] += int(re.search('[0-9]+', bases[pos:]).group())
                    # move by position to next location
                    pos += len(indel)
                    # cov['D'] += int(re.search('[0-9]+', bases[pos:]).group())
                elif base == '*':
                    cov['D'] += 1


            base_values = {}
            base_values['A'] = cov['A'] + cov['a']
            base_values['C'] = cov['C'] + cov['c']
            base_values['G'] = cov['G'] + cov['g']
            base_values['T'] = cov['T'] + cov['t']
            base_values['D'] = cov['D']

            max_value = max(base_values, key=base_values.get)
            position_value = ""


            if max_value in "ACGT" and max(cov[max_value], cov[max_value.lower()]) / 2 < min(cov[max_value], cov[max_value.lower()]) and max(cov[max_value], cov[max_value.lower()]) > args.reads:
                position_value = max_value
            elif max_value != "D":
                position_value = ref_upper


            if first_time:
                first_time = False
                old_chr = contig
                data_out.write(">{}\n{}".format(contig, position_value))
                n = 1
            elif not first_time and old_chr == contig:
                if n == args.base:
                    data_out.write("\n{}".format(position_value))
                    n = 1
                else:
                    data_out.write(position_value)
                    n += 1
            else:
                data_out.write("\n>{}\n{}".format(contig, position_value))
                old_chr = contig
                n = 1

def main():
    args = get_args()



if __name__ == "__main__":
    main()

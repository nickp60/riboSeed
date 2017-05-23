#!/usr/bin/env python3.5
"""
"""
import argparse
import os
import sys
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


def get_args():  # pragma: no cover
    """get the arguments as a main parser with subparsers
    for named required arguments and optional arguments
    """
    parser = argparse.ArgumentParser(
        description="This is used to combine regions extracted with " +
        "riboSnag, creating a single sequence")
    parser.add_argument("indir", help="Directory with fasta's to concatenate")
    parser.add_argument("ext", help="Extension of files to concatenate")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output",
                               help="output directory;"
                               "default: %(default)s",
                               default=os.getcwd(),
                               type=str, dest="output")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-n", "--name",
                          help="name for output fasta",
                          default='concatenated_seq', type=str)

    args = parser.parse_args()
    return args


def concat_genome(input_dir, ext, outpath, verbose=False):
    """for each fasta, read in, add to existing string, and when finished,
    write out as single-entry fasta
    """
    fastas = sorted(glob.glob(str(input_dir + ext)))
    if len(fastas) == 0:
        if verbose:
            print("No files found!")
        return 1
    if verbose:
        print(str("combining the following files matching extension " +
                  "{0}:{1}".format(ext, " ".join(fastas))))
    new_seq = ""
    for filen in fastas:
        print("Adding %s to combined sequence" % filen)
        with open(filen, 'r') as i_file:
            seq_rec = list(SeqIO.parse(i_file, 'fasta'))[0]
            new_seq = new_seq + str(seq_rec.seq)
            if verbose:
                print(str("Len of sequence:{0}\nLen of concatenated " +
                          "sequence:{1}").format(len(seq_rec),
                                                 len(new_seq)))
    try:

        with open(outpath, 'w') as o_file:
            success = SeqIO.write(
                SeqRecord(
                    seq=Seq(new_seq, IUPAC.IUPACAmbiguousDNA()),
                    id="concatenated_genome"), o_file, 'fasta')
            if success:
                print("wrote out concatenated file!")
                return 0
    except Exception as e:
        if verbose:
            print(e)
        return 1


if __name__ == "__main__":
    args = get_args()
    print(args)
    try:
        os.makedirs(args.output)
    except OSError:
        print("Output directory already exists.  Exiting")
        sys.exit(1)

    res = concat_genome(input_dir=args.indir, ext=args.ext,
                        outpath=os.path.join(args.output,
                                             str(args.name + ".fasta")),
                        verbose=True)
    if res == 1:
        print("exiting")
    else:
        print("Have a fantastic day!")

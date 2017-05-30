#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
"""

import argparse
import sys
import time
import logging
import os
import traceback
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# need this line for unittesting
sys.path.append(os.path.join('..', 'riboSeed'))
from pyutilsnrw.utils3_5 import set_up_logging

# --------------------------- classes --------------------------- #


# --------------------------- methods --------------------------- #


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="This uses the JC model of genome evolution " +
        "to test the effectiveness of riboSeed on divergent " +
        "reference sequences.",
        add_help=False)  # to allow for custom help
    parser.add_argument("fasta", action="store",
                        help="(multi)fasta file containing the sequences to " +
                        "be mutated")

    # taking a hint from http://stackoverflow.com/questions/24180527
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("-f", "--frequency", dest='frequency',
                          action="store",
                          default=.01, type=float,
                          help="Probability of mutated bases" +
                          "default: %(default)s")
    optional.add_argument("-e", "--end_length", dest='end_length',
                          action="store",
                          default=None, type=int,
                          help="if value given, only mutated the ends " +
                          "of the sequences and ignore the middle" +
                          "default: %(default)s")
    optional.add_argument("-s", "--seed", dest='seed',
                          action="store",
                          default=27, type=int,
                          help="cause reproduciblity; " +
                          "default: %(default)s")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    args = parser.parse_args()
    return args


def last_exception():
    """ Returns last exception as a string, or use in logging.
    stolen verbatim from pyani
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


def ageSequence(rec, outfile, freq, end_length, logger=None):
    assert logger is not None, "must use logging"
    logger.info("frequncy of mutation: %f", freq)
    change_counter = 0
    distr_counter = 0
    hist = []
    newseqlist = list(rec.seq)
    # print(rec.seq[1:10])
    # print(newseqlist[1:10])
    alph = ["A", "T", "C", "G"]
    seqlen = len(rec.seq)
    if end_length is None:
        ignore_region = []
    elif not seqlen - (2 * end_length) > 1:
        raise ValueError("Edge width cannot be greater than half the " +
                         "length of the sequence ")
    else:
        ignore_region = set([idx for sublist in
                             [range(end_length, seqlen - end_length)]
                             for idx in sublist])
    for i, base in enumerate(rec.seq):
        if i % 5000 == 0:
            logger.debug("processing base %d of %d", i, seqlen)
        if i in ignore_region:
            hist.append(0)
        else:
            distr_counter = distr_counter + 1
            # this should be the geometric distribution
            if random.random() < (1 - ((1 - freq) ** distr_counter)):
                hist.append(1)
                distr_counter = 0
                change_counter = change_counter + 1
                choices = alph[:]
                choices.remove(base)
                newseqlist[i] = random.choice(choices)
            else:
                hist.append(0)
    # print(hist[1:20])
    logger.info("Changed %d of %d bases", change_counter, i)
    newrec = SeqRecord(
        id=rec.id,
        # description="riboSim mutation frequency" + str(freq),
        seq=Seq("".join(newseqlist),
                IUPAC.IUPACAmbiguousDNA()))

    with open(outfile, "a") as o:
        SeqIO.write(newrec, o, "fasta")

    assert len(newseqlist) == len(rec.seq), "something bad happened!"


if __name__ == "__main__":  # pragma: no cover
    args = get_args()
    # allow user to give relative paths
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory already exists; exiting...")
        sys.exit(1)
    t0 = time.time()
    log_path = os.path.join(output_root,
                            str("riboSim_log_{0}.txt".format(
                                time.strftime("%Y%m%d%H%M"))))
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=log_path,
                            name=__name__)

    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))

    output_file = os.path.join(
        output_root,
        os.path.basename(args.fasta))

    random.seed(args.seed)
    with open(args.fasta, "r") as infile:
        for rec in SeqIO.parse(infile, "fasta"):
            ageSequence(rec, freq=args.frequency, end_length=args.end_length,
                        outfile=output_file, logger=logger)
    # Report that we've finished
    logger.info("Done: %s", time.asctime())
    logger.info("Time taken: %.3fs" % (time.time() - t0))

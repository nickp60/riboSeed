#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""
"""

import argparse
import sys
import time
import logging
import os
import traceback
from numpy import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from .shared_methods import set_up_logging

# --------------------------- classes --------------------------- #


# --------------------------- methods --------------------------- #


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="This uses the JC model of genome evolution " +
        "to test the effectiveness of riboSeed on divergent " +
        "reference sequences.",
        add_help=False)  # to allow for custom help
    parser.prog = "ribo sim"
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
                          default=None, type=int,
                          help="cause reproduciblity; " +
                          "default: %(default)s")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    args = parser.parse_args(sys.argv[2:])
    return args


def last_exception():
    """ Returns last exception as a string, or use in logging.
    stolen verbatim from pyani
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


def substitute_base(strlist, position, alph):
    """  string should be converted to list prior for mutability"""
    oldbase = strlist[position]
    choices = alph[:]
    choices.remove(oldbase)
    strlist[position] = random.choice(choices)


def ageSequence(rec, outfile, freq, end_length, seed, logger=None):
    assert logger is not None, "must use logging"
    logger.info("frequncy of mutation: %f", freq)
    newseqlist = list(rec.seq)
    alph = ["A", "T", "C", "G"]
    seqlen = len(rec.seq)
    if end_length is None or end_length is 0:
        ignore_region = []
    elif not seqlen - (2 * end_length) > 1:
        raise ValueError("Edge width cannot be greater than half the " +
                         "length of the sequence ")
    else:
        ignore_region = set([idx for sublist in
                             [range(end_length, seqlen - end_length)]
                             for idx in sublist])
    newseqlist = list(rec.seq)
    seq_len = len(newseqlist)
    random.seed(seed)
    # subst_idxs = random.sample(range(0, seq_len), int(round(seq_len * freq)))
    idxs = list(range(0, seq_len))
    random.shuffle(idxs)
    subst_idxs = idxs[0: int(round(seq_len * freq))]
    # ignore the indexes in the regions we are leaving unchanaged
    executed_subst_idxs = [x for x in subst_idxs if x not in ignore_region]
    for i in executed_subst_idxs:
        if i in ignore_region:
            pass
        else:
            substitute_base(strlist=newseqlist, position=i, alph=alph)
    logger.info("Changed %d of %d bases", len(executed_subst_idxs), seq_len)
    newrec = SeqRecord(
        id=rec.id,
        # description="riboSim mutation frequency" + str(freq),
        seq=Seq("".join(newseqlist),
                IUPAC.IUPACAmbiguousDNA()))

    with open(outfile, "a") as o:
        SeqIO.write(newrec, o, "fasta")

    assert len(newseqlist) == len(rec.seq), \
        "something bad happened! unequal lengths of input and output sequences"


def main(args, logger=None):
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

    with open(args.fasta, "r") as infile:
        for rec in SeqIO.parse(infile, "fasta"):
            ageSequence(rec, freq=args.frequency, end_length=args.end_length,
                        outfile=output_file, seed=args.seed, logger=logger)
    # Report that we've finished
    logger.info("Done: %s", time.asctime())
    logger.info("Time taken: %.3fs" % (time.time() - t0))

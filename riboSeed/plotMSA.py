#!/usr/bin/env python
"""
Minor version changes:
 - cleaned up
Input:
- genbank file
- dictionary
- specific features : 16S, 5S
- upstream, downstream widths

Output:
-dir containing DNA fastas in their

"""
import os
import subprocess
import datetime
import time
import argparse
import sys
import math
import re
import shutil
import itertools
import multiprocessing

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from collections import defaultdict  # for calculating kmer frequency
from itertools import product  # for getting all possible kmers
# from heatmapcluster import heatmapcluster

#from pyutilsnrw import utils3_5
sys.path.append(os.path.join('..', 'riboSeed'))

from pyutilsnrw.utils3_5 import set_up_logging, check_installed_tools,\
    combine_contigs, check_version_from_cmd

from riboSnag import calc_entropy_msa, annotate_msa_conensus, \
    plot_scatter_with_anno

def get_args():  # pragma: no cover
    """get the arguments as a main parser with subparsers
    for named required arguments and optional arguments
    """
    parser = argparse.ArgumentParser(description="riboSnag lite.  Use this " +
                                     "to generate the plots made in riboSnag",
                                     add_help=False)
    parser.add_argument("msa", help="msa in fasta format")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output",
                               help="output directory; default: %(default)s",
                               default=os.getcwd(),
                               type=str, dest="output")

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-n", "--name",
                          help="rename the contigs with this prefix" +
                          # "default: %(default)s",
                          "default: date (YYYYMMDD)",
                          default=None, dest="name",
                          action="store", type=str)
    optional.add_argument("-v", "--verbosity",
                          dest='verbosity', action="store",
                          default=2, type=int,
                          help="1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("--barrnap_exe", dest="barrnap_exe",
                          action="store", default="barrnap",
                          help="Path to barrnap executable; " +
                          "default: %(default)s")
    optional.add_argument("--kingdom", dest="kingdom",
                          action="store", default="bac",
                          choices=["mito", "euk", "arc", "bac"],
                          help="kingdom for barrnap; " +
                          "default: %(default)s")
    optional.add_argument("-s", "--seq_name", dest='seq_name',
                          action="store",
                          help="name of genome; "
                          "default: inferred from file name, as many cases " +
                          "involve multiple contigs, etc, making inference " +
                          "from record intractable;" +
                          "default: %(default)s",
                          default="plotMSA",
                          type=str)
    # had to make this explicitly to call it a faux optional arg
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    output_root = os.path.abspath(os.path.expanduser(args.output))
    # Create output directory only if it does not exist
    try:
        os.makedirs(args.output)
    except FileExistsError:
        # '#' is needed in case streaming output eventually
        print("#Selected output directory %s exists" %
              args.output)
        print("exiting")
        sys.exit(1)
    log_path = os.path.join(output_root,
                            str("{0}_makeMSA_log.txt".format(
                                time.strftime("%Y%m%d%H%M"))))
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=log_path,
                            name=__name__)

    logger.debug("Usage:\n{0}\n".format(str(" ".join([x for x in sys.argv]))))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, v)
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    # test whether executables are there
    executables = [args.barrnap_exe]
    test_ex = [check_installed_tools(x, logger=logger) for x in executables]
    if all(test_ex):
        logger.debug("All needed system executables found!")
        logger.debug(str([shutil.which(i) for i in executables]))

    # make MSA and calculate entropy
    seq_entropy, names, tseq = calc_entropy_msa(args.msa)
    gff, consensus_cov, annos = annotate_msa_conensus(
        tseq_array=tseq,
        pattern='product=(.+?)$',
        seq_file=os.path.join(
            output_root,
            "test_consensus.fasta"),
        barrnap_exe=args.barrnap_exe,
        countdashcov=False,   # include -'s in coverage
        collapseNs=False,  # include -'s in consensus
        excludedash=False,
        kingdom=args.kingdom,
        logger=logger)
    label_name = args.seq_name
    title = str("Shannon Entropy by Position\n" +
                label_name)

    return_code = plot_scatter_with_anno(
        data=seq_entropy,
        consensus_cov=consensus_cov,
        names=["Position", "Entropy"],
        title=title,
        anno_list=annos,
        output_prefix=os.path.join(
            output_root,
            "entropy_plot"))

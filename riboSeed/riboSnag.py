#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

import os
import subprocess
import datetime
import argparse
import sys
import math
import re
import shutil
import itertools
import multiprocessing

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import \
        FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib import gridspec
    import matplotlib.patches as patches
    # matplotlib.rc('font', family='sans-serif')
    # matplotlib.rcParams['text.usetex'] = True
    # matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
    # matplotlib.rcParams['font.family'] = 'sans-serif'
    # matplotlib.rcParams['font.sans-serif'] = 'cm'
    PLOT = True
except Exception as e:  # likely an ImportError, but not taking chances
    print(e)
    print("\nlooks like you have some issue with matplotlib.  " +
          "Classic matplotlib, amirite? Plotting is disabled\n")
    PLOT = False

import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from collections import defaultdict  # for calculating kmer frequency
from itertools import product  # for getting all possible kmers
# from heatmapcluster import heatmapcluster


from .classes import LociCluster, Locus

from .shared_methods import add_gb_seqrecords_to_cluster_list, \
    extract_coords_from_locus, pad_genbank_sequence, \
    parse_clustered_loci_file, set_up_logging, \
    combine_contigs, check_version_from_cmd

def get_args(test_args=None):  # pragma: no cover
    """get the arguments as a main parser with subparsers
    for named required arguments and optional arguments
    """
    parser = argparse.ArgumentParser(prog="ribo snag",
                                     description="Use to extract regions " +
                                     "of interest based on supplied Locus " +
                                     "tags and evaluate the extracted regions",
                                     add_help=False)
    parser.prog = "ribo snag"
    parser.add_argument("genbank_genome", help="Genbank file (WITH SEQUENCE)")
    parser.add_argument("clustered_loci", help="output from riboSelect")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output",
                               help="output directory; default: %(default)s",
                               default=os.getcwd(),
                               type=str, dest="output")

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-n", "--name",
                          help="rename the contigs with this prefix; " +
                          # "default: %(default)s",
                          "default: date (YYYYMMDD)",
                          default=None, dest="name",
                          action="store", type=str)
    optional.add_argument("-l", "--flanking_length",
                          help="length of flanking regions, in bp; " +
                          "default: %(default)s",
                          default=1000, type=int, dest="flanking")
    optional.add_argument("--msa_kmers",
                          help="calculate kmer similarity based on aligned " +
                          "sequences instead of raw sequences;" +
                          "default: %(default)s",
                          default=False, action="store_true", dest="msa_kmers")
    optional.add_argument("--skip_kmers",
                          help="Just plot entropy if MSA" +
                          "default: %(default)s",
                          default=False, action="store_true", dest="skip_kmers")
    optional.add_argument("--skip_blast",
                          help="Skip running BLAST Comparisons" +
                          "default: %(default)s",
                          default=False, action="store_true", dest="skip_blast")
    # optional.add_argument("-c", "--circular",
    #                       help="if the genome is known to be circular, and " +
    #                       "an region of interest (including flanking bits) " +
    #                       "extends past chromosome end, this extends the " +
    #                       "sequence past chromosome origin forward by 5kb; " +
    #                       "default: %(default)s",
    #                       default=False, dest="circular", action="store_true")
    optional.add_argument("--linear",
                          help="if the genome is not circular, and " +
                          "an region of interest (including flanking bits) " +
                          "extends past chromosome end, this extends the " +
                          "sequence past chromosome origin forward by 5kb; " +
                          "default: %(default)s",
                          default=False, dest="linear", action="store_true")
    optional.add_argument("-p", "--padding", dest='padding', action="store",
                          default=5000, type=int,
                          help="if treating as circular, this controls the " +
                          "length of sequence added to the 5' and 3' ends " +
                          "to allow for selecting regions that cross the " +
                          "chromosome's origin; default: %(default)s")
    optional.add_argument("-v", "--verbosity",
                          dest='verbosity', action="store",
                          default=2, type=int,
                          help="1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("--title",
                          help="String for plot title;" +
                          " uses matplotlib math processing for italics " +
                          "(you know, the LaTeX $..$ syntax): " +
                          "https://matplotlib.org/users/mathtext.html " +
                          "default: inferred from --seq_name",
                          action='store',
                          default=None,
                          dest="title")
    optional.add_argument('--pubplot',
                          action="store_true",
                          # want slightly cleaner figures? Try this for
                          # bigger fonts, adjusted label spacing, and
                          help=argparse.SUPPRESS)
    optional.add_argument("--clobber",
                          help="overwrite previous output files; " +
                          "default: %(default)s", action='store_true',
                          default=False, dest="clobber")
    optional.add_argument("--no_revcomp",
                          help="default returns reverse complimented seq " +
                          "if majority of regions on reverse strand. if  " +
                          "--no_revcomp, this is overwridden; " +
                          "default: %(default)s",
                          action='store_true',
                          default=False, dest="no_revcomp")
    optional.add_argument("-j", "--just_extract",
                          help="Dont bother making an MSA, calculating " +
                          "Shannon Entropy, BLASTing, generating plots etc; " +
                          " just extract the regions ; " +
                          "default: %(default)s",
                          action='store_true',
                          default=False, dest="just_extract")
    optional.add_argument("--msa_tool", dest="msa_tool",
                          choices=["mafft", "prank"],
                          action="store", default="mafft",
                          help="Path to PRANK executable; " +
                          "default: %(default)s")
    optional.add_argument("--prank_exe", dest="prank_exe",
                          action="store", default="prank",
                          help="Path to PRANK executable; " +
                          "default: %(default)s")
    optional.add_argument("--mafft_exe", dest="mafft_exe",
                          action="store", default="mafft",
                          help="Path to MAFFT executable; " +
                          "default: %(default)s")
    optional.add_argument("--barrnap_exe", dest="barrnap_exe",
                          action="store", default="barrnap",
                          help="Path to barrnap executable; " +
                          "default: %(default)s")
    optional.add_argument("--makeblastdb_exe", dest="makeblastdb_exe",
                          action="store", default="makeblastdb",
                          help="Path to makeblastdb executable; " +
                          "default: %(default)s")
    optional.add_argument("--kingdom", dest="kingdom",
                          action="store", default="bac",
                          choices=["mito", "euk", "arc", "bac"],
                          help="kingdom for barrnap; " +
                          "default: %(default)s")
    optional.add_argument("-s", "--seq_name", dest='seq_name',
                          action="store",
                          help="name of genome; "
                          "default: inferred from file name, as many cases" +
                          "involve multiple contigs, etc, making inference  " +
                          "from record intractable",
                          type=str)
    # had to make this explicitly to call it a faux optional arg
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    if test_args is None:
        args = parser.parse_args(sys.argv[2:])
    else:
        args = parser.parse_args(test_args)
    return args


def check_loci_file_not_genbank(filepath):
    """ raise error if "cluster file" ends in .gb, etc

    this covers common case where user submits genbank and cluster file
    in the wrong order.

    Args:
        filepath (str): path to what should be a loci file
    Returns:
        None
    Raises:
       FileNotFoundError: if the file smells like a genbank file

    """
    if not (os.path.isfile(filepath) and os.path.getsize(filepath) > 0):
        raise ValueError("Cluster File not found!")
    if filepath.endswith(("gb", "genbank", "gbk")):
        raise FileNotFoundError("Hmm, this cluster file looks like genbank;" +
                                " it ends in {0}".format(os.path.splitext(
                                    filepath)[1]))


def stitch_together_target_regions(cluster,
                                   flanking=500,
                                   logger=None, circular=False,
                                   revcomp=False):
    """
    given a single LociCluster object containing Locus objects (usually)
    of length 3 (16,5,and 23 rRNAs),
    return the object with ammended sequence info
    revamped 20161004
    """
    assert logger is not None, "Must have logger for this function"

    # TODO : make this safer. coord list is constructed sequentially but this
    # is a backup. Throws sort of a cryptic error. but as I said, its a backup
    if sorted([x.start_coord for x in cluster.loci_list]) != \
       [x.start_coord for x in cluster.loci_list]:
        logger.warning("Coords are not in increasing order; " +
                       "you've been warned")
    start_list = sorted([x.start_coord for x in cluster.loci_list])
    logger.debug("Start_list: {0}".format(start_list))

    logger.debug("stitching together the following coords:")
    for i in cluster.loci_list:
        logger.debug(str(i.__dict__))
    #  This works as long as coords are never in reverse order
    cluster.global_start_coord = min([x.start_coord for
                                      x in cluster.loci_list]) - flanking
    # if start is negative, just use 1, the beginning of the sequence
    if cluster.global_start_coord < 1:
        logger.warning("Caution! Cannot retrieve full flanking region, as " +
                       "the 5' flanking region extends past start of " +
                       "sequence. If this is a problem, try using a smaller " +
                       "--flanking region, and/or if  appropriate, run " +
                       "without with --linear.")
        cluster.global_start_coord = 1
    cluster.global_end_coord = max([x.end_coord for
                                    x in cluster.loci_list]) + flanking
    if cluster.global_end_coord > len(cluster.seq_record):
        logger.warning("Caution! Cannot retrieve full flanking region, as " +
                       "the 5' flanking region extends past start of " +
                       "sequence. If this is a problem, try using a smaller " +
                       "--flanking region, and/or if  appropriate, run " +
                       "without with --linear.")
        cluster.global_end_coord = len(cluster.seq_record)

    logger.debug("global start and end: %s %s", cluster.global_start_coord,
                 cluster.global_end_coord)
    #  the minus one makes things go from 1 based to zero based
    seq_with_ns = str(cluster.seq_record.seq[cluster.global_start_coord - 1:
                                             cluster.global_end_coord])
    seq_len = len(seq_with_ns[:])
    logger.debug(seq_len)
    # again, plus 1 corrects for 0 index.
    # len("AAAA") = 4 vs AAAA[-1] - AAAA[0] = 3
    logger.debug("\nexp length %i \nact length %i",
                cluster.global_end_coord - cluster.global_start_coord + 1,
                seq_len)

    ## Change coords in ID When using padding
    if not circular:
        seq_id = str(cluster.sequence_id + "_" + str(cluster.global_start_coord) +
                     ".." + str(cluster.global_end_coord))
    else:  # correct for padding
        seq_id = str(cluster.sequence_id + "_" + str(cluster.global_start_coord -
                                                     cluster.padding) +
                     ".." + str(cluster.global_end_coord - cluster.padding))

    # if most are on - strand, return sequence reverse complement
    strands = [x.strand for x in cluster.loci_list]
    if revcomp and \
       (sum([x == -1 for x in strands]) > sum([x == 1 for x in strands])):
        logger.info("returning the reverse compliment of the sequence")
        cluster.extractedSeqRecord = SeqRecord(
            Seq(seq_with_ns,
                IUPAC.IUPACAmbiguousDNA()).reverse_complement(),
            id=str(seq_id + "_RC"))
    else:
        cluster.extractedSeqRecord = SeqRecord(
            Seq(seq_with_ns,
                IUPAC.IUPACAmbiguousDNA()),
            id=seq_id)
    ### last minuete check
    cluster.extractedSeqRecord.description = "from riboSnag"
    for prop, value in vars(cluster).items():
        if value is None:
            logger.debug("%s has a value of None!", prop)
    return cluster


def prepare_prank_cmd(outdir, combined_fastas, prank_exe,
                      add_args="", outfile_name="best_MSA",
                      logger=None):
    """returns command line for constructing MSA with
    PRANK and the path to results file
    """
    assert logger is not None, "Must use logger"
    if not os.path.exists(outdir):
        raise FileNotFoundError("output directory not found!")
    prank_cmd = "{0} {1} -d={2} -o={3}".format(
        prank_exe, add_args, combined_fastas,
        os.path.join(outdir, outfile_name))
    logger.debug("PRANK command: \n %s", prank_cmd)
    return (prank_cmd, os.path.join(outdir, str(outfile_name + ".best.fas")))


def prepare_mafft_cmd(outdir, combined_fastas, mafft_exe,
                      add_args="", outfile_name="best_MSA",
                      logger=None):
    """returns command line for constructing MSA with
    mafft and the path to results file.
    We call it best.fas to match the mafft result
    """
    assert logger is not None, "Must use logger"
    if not os.path.exists(outdir):
        raise FileNotFoundError("output directory not found!")
    mafft_cmd = "{0} {1} {2} > {3}".format(
        mafft_exe, add_args, combined_fastas,
        os.path.join(outdir, outfile_name + ".best.fas"))
    logger.debug("MAFFT command: \n %s", mafft_cmd)
    return (mafft_cmd, os.path.join(outdir, outfile_name + ".best.fas"))


def calc_Shannon_entropy(matrix):
    """ $j$ has entropy $H(j)$ such that
    $H(j) = -sum_{i=(A,C,T,G)} p_i(j) log p_i(j)$
    """
    entropies = []
    for instance in matrix:
        unique = set(instance)
        proportions = {}
        for i in unique:
            proportions[i] = sum([x == i for x in instance]) / len(instance)
        entropy = -sum([prob * (math.log(prob, math.e)) for
                        prob in proportions.values()])
        entropies.append(entropy)
    return entropies


def calc_entropy_msa(msa_path):
    """givn a path to an MSA in FASTA format, this gets the
    entropy of each position in batches, so as to not gobble memory
    return list
    """
    batch_size = 1000  # read seequences in chunks this long
    lengths = []
    with open(msa_path) as fh:
        msa_seqs = list(SeqIO.parse(fh, 'fasta'))
    seq_names = []
    for rec in msa_seqs:
        lengths.append(len(rec))
        seq_names.append(rec.id)
    if not all([i == lengths[0] for i in lengths]):
        raise ValueError("Sequences must all be the same length!")
    entropies = []
    tseq = []
    # save memory by reading in chunks
    for batch in range(0, (math.ceil(lengths[0] / batch_size))):
        # get all sequences into an array
        seq_array = []
        for nseq, record in enumerate(msa_seqs):
            seq_array.append(
                [x for x in record.seq[(batch * batch_size):
                                       ((batch + 1) * batch_size)]])
        # transpose
        tseq_array = list(map(list, zip(*seq_array)))
        tseq.extend(tseq_array)
        entropies.extend(calc_Shannon_entropy(tseq_array))
    # check length of sequence is the same as length of the entropies
    assert len(entropies) == lengths[0]
    return (entropies, seq_names, tseq)


def annotate_msa_conensus(tseq_array, seq_file, barrnap_exe,
                          kingdom="bact",
                          pattern='product=(.+?)$',
                          countdashcov=True,   # include -'s in coverage
                          collapseNs=False,  # include -'s in consensus
                          excludedash=False,
                          logger=None):
    """ returns annotations (as gfflist),the consensus sequence as a
    list[base, cov], and named coords  as a list
    TODO: The 'next_best' thing fails is an N is most frequent. Defaults to a T
    """
    assert logger is not None, "Must use logger"
    if excludedash:
        logger.warning("CAUTION: excludedash selected. There is a known " +
                       "bug in the 'next_best' thing fails if an " +
                       "N is most frequent. Defaults to a T")
    consensus = []
    nseqs = len(tseq_array[0])
    logger.info("calc coverage for each of the %i positions", len(tseq_array))
    for position in tseq_array:
        if all([x == position[0] for x in position]):
            if position[0] == '-':
                if collapseNs:
                    continue
                elif not countdashcov:
                    consensus.append([position[0], 0])
            else:
                consensus.append([position[0], nseqs])
        else:
            max_count = 0  # starting count
            nextbest_count = 0
            best_nuc = None
            nextbest_nuc = None
            for nuc in set(position):
                count = sum([nuc == z for z in position])
                # if max count, swap with max to nextbest and update max
                if count > max_count:
                    nextbest_count = max_count
                    max_count = count  # update count if better
                    nextbest_nuc = best_nuc
                    best_nuc = nuc
                else:
                    pass
            # contol whether gaps are allowed in consensus
            if (
                    all([x != '-' for x in position]) and
                    best_nuc == '-' and
                    excludedash):
                # if we dont want n's, choose nextbest
                if nextbest_nuc is None:
                    nextbest_nuc = 't'  # I hate myself for this
                consensus.append([nextbest_nuc, nextbest_count])
            elif best_nuc == '-' and not countdashcov:
                consensus.append([best_nuc, 0])
            else:
                consensus.append([best_nuc, max_count])
    # if any are '-', replace with n's for barrnap
    seq = str(''.join([x[0] for x in consensus])).replace('-', 'n')

    # annotate seq
    with open(seq_file, 'w') as output:
        SeqIO.write(SeqRecord(
            Seq(seq, IUPAC.IUPACAmbiguousDNA())), output, "fasta"),
    barrnap_cmd = "{0} --kingdom {1} {2}".format(barrnap_exe,
                                                 kingdom, seq_file)
    barrnap_gff = subprocess.run(barrnap_cmd,
                                 shell=sys.platform != "win32",
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 check=True)
    results_list = [x.split('\t') for x in
                    barrnap_gff.stdout.decode("utf-8").split("\n")]

    ###  make [name, [start_coord, end_coord]] list
    named_coords = []
    for i in results_list:
        if i[0].startswith("#") or not len(i) == 9:
            logger.debug("skipping gff line: %s", i)
            continue
        m = re.search(pattern, i[8])
        if m:
            found = m.group(1)
        named_coords.append([found, [int(i[3]), int(i[4])]])
    if len(named_coords) == 0:
        raise ValueError(str("Error extracting coords from barrnap gff" +
                             "with pattern %s!" %  pattern))
    return (results_list, consensus, named_coords)


def plot_scatter_with_anno(data,
                           consensus_cov,
                           anno_list,
                           names=["Position", "Entropy"],
                           title="$Strain Name$",
                           output_prefix="entropy_plot.png", pubplot=True):
    """Given annotation coords [feature, [start, end]],
    consensus cov list ['base', coverage_depth],
    entropy values (list) and consensus sequence
    (same length for consensus_cov and data, no funny business),
    plot out the entropies for each position,
    plot the annotations, and return 0
    """
    if len(consensus_cov) != len(data):
        raise ValueError("data and consensus different lengths!")
    df = pd.DataFrame({names[0]: range(1, len(data) + 1),
                       names[1]: data})  # columns=names)
    df_con = pd.DataFrame(consensus_cov, columns=["base",
                                                  "depth"])
    cov_max_depth = max(df_con['depth'])
    fig = Figure()
    FigureCanvas(fig)
    #fig.SubplotParams(hspace=0.0)
    # fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
    #                                gridspec_kw={'height_ratios': [4, 1]})
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax1.get_shared_x_axes().join(ax1, ax2)
    colors = ['#ff4c05', '#FFFB07', '#04FF08', '#06B9FF', '#6505FF', '#FF012F',
              '#ff4c05', '#FFFB07', '#04FF08', '#06B9FF', '#6505FF', '#FF012F']
    # ax1.set_title("Shannon Entropy by Position\n" +
    #               title, y=1.08, fontsize=20)
    if pubplot:
        fontscale = 1.5
    else:
        fontscale = 1
    ax1.set_title(title, y=1.08, fontsize=18 * fontscale)
    xmin, xmax = 0, len(data)
    ymin, ymax = -0.1, (max(data) * 1.2)
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
    if pubplot:
        # make minimal y labels for coeverage
        ax2.set_yticks([0, max(df_con["depth"])])
    yjust = -.1
    # add the highlighted bits sowing where te coding regions are
    for index, anno in enumerate(anno_list):
        rect1 = patches.Rectangle(
            (anno[1][0],  # starting x
             ymin),  # starting y
            anno[1][1] - anno[1][0],  # rel x end
            ymax - ymin,  # rel y end
            facecolor=matplotlib.colors.ColorConverter().to_rgba(
                colors[index], alpha=0.2),
            edgecolor=matplotlib.colors.ColorConverter().to_rgba(
                colors[index], alpha=0.2))
        rect2 = patches.Rectangle(
            (anno[1][0],  # starting x
             1),  # starting y
            anno[1][1] - anno[1][0],  # rel x end
            cov_max_depth,  # dont -1 beacuse start at 1
            facecolor=matplotlib.colors.ColorConverter().to_rgba(
                colors[index], alpha=0.2),
            edgecolor=matplotlib.colors.ColorConverter().to_rgba(
                colors[index], alpha=0.2))
        ax1.add_patch(rect1)
        ax2.add_patch(rect2)
        if not pubplot:
            ax1.text((anno[1][0] + anno[1][1]) / 2,    # x location
                     ymax - 0.48 - yjust,              # y location
                     anno[0][0:20],                     # text first 20 char
                     ha='center', color='red', weight='bold', fontsize=11)
            yjust = yjust * - 1
    ax1.scatter(x=df["Position"], y=df["Entropy"],
                marker='o', color='black', s=2)
    # add smoothing for kicks
    df["fit"] = savitzky_golay(df["Entropy"].values, 351, 3)  # window size 51, polynomial order 3
    ax1.scatter(x=df["Position"], y=df["fit"], color='red', s=1)
    #
    ax1.set_ylabel('Shannon Entropy', fontsize= 13 * fontscale)
    ax1.get_yaxis().set_label_coords(-.05, 0.7)
    ax2.set_xlim([xmin, xmax])
    ax2.invert_yaxis()  # we want the bars pointing down
    ax2.set_ylabel('Consensus Coverage', fontsize= 13 * fontscale)
    ax2.set_xlabel('Position (bp)', fontsize= 13 * fontscale)
    ax2.get_yaxis().set_label_coords(-.05, 0.5)
    # ax2.set_ylim([1, cov_max_depth + 1]) #, 1])
    ax2.bar(df_con.index, df_con["depth"],
            width=1, color='darkgrey', linewidth=0, edgecolor='darkgrey')
    # ax2.step(df_con.index, df_con["depth"],
    #          where='mid', color='darkgrey')
    for ax in [ax1, ax2]:
        ax.spines['right'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax1.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')
    ax1.xaxis.set_ticks_position('top')
    ax1.tick_params(axis='y', colors='dimgrey', labelsize = 10 * fontscale)
    ax2.tick_params(axis='y', colors='dimgrey', labelsize = 10 * fontscale)
    ax1.tick_params(axis='x', colors='dimgrey', labelsize = 10 * fontscale)
    ax2.tick_params(axis='x', colors='dimgrey', labelsize = 10 * fontscale)
    ax1.yaxis.label.set_color('black')
    ax2.yaxis.label.set_color('black')
    ax1.xaxis.label.set_color('black')
    ax2.xaxis.label.set_color('black')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.set_size_inches(12, 7.5)
    fig.savefig(str(output_prefix + '.png'), dpi=(300))
    fig.savefig(str(output_prefix + '.pdf'), dpi=(300))
    return 0


def get_all_kmers(alph="", length=3):
    """Given an alphabet of charactars, return a list of all permuations
    not actually used, I dont think I'll need this
    """
    mers = [''.join(p) for p in product(alph, repeat=length)]
    return mers


def profile_kmer_occurances(rec_list, k, logger=None):
    """ given a list of seq records, an alphabet, and a value k,
    retrun counts dict of kmer occurances and list of seq names
    """
    counts = defaultdict(list)
    names_list = []
    # part 1: get all kmers from seqs so that we can have equal length lisst
    all_mers = []
    for rec in rec_list:
        all_mers.extend([str(rec.seq).lower()[x: x + k] for
                         x in range(0, len(rec.seq) - k + 1)])
    # unique_mers = list(set().union(all_mers))
    unique_mers = set(all_mers)
    for i in unique_mers:
        counts[i] = []  # initialixe counts dictionary with ker keys
    # part two: count 'em
    for n, rec in enumerate(rec_list):
        logger.info("counting kmer occurances for %s", rec.id)
        names_list.append(rec.id)
        logger.debug("converting to lower")
        string = str(rec.seq).lower()
        logger.debug("getting %imers from seq", k)
        string_mers = [string[x: x + k] for x in range(0, len(string) - k + 1)]
        logger.debug("counting mer occurances")
        # Add counts for those present
        for value in set(string_mers):
            counts[value].append(sum([value == mer for mer in string_mers]))
        # filling in where not found
        for missing in unique_mers - set(string_mers):
            counts[missing].append(0)
    return counts, names_list


def plot_pairwise_least_squares(counts, names_list, output_prefix):
    """given a  list of sequence names, a prefix for plot file names,
    and a list of counts from plot_kmer_occurances,
    retruns a pandas df of least squares after plotting heatmaps
    """
    res_list = []
    counts_list = []
    for v in counts.values():
        counts_list.append(v)
    # this gives each an index and gets all the pairs
    all_pairs = [[index, value] for index, value in
                 enumerate(product(range(0, len(names_list)), repeat=2))]
    for i in all_pairs:
        this_pairs_diffs = []
        for row in counts_list:
            this_pairs_diffs.append((row[i[1][0]] - row[i[1][1]]) ** 2)
        res_list.append([names_list[i[1][0]],
                         names_list[i[1][1]],
                         sum(this_pairs_diffs)])
    lsdf_wNA = pd.DataFrame(res_list, columns=["locus_1", "locus_2", "sls"])
    wlsdf = lsdf_wNA.pivot(index='locus_1', columns='locus_2', values='sls')
    fig = Figure()
    FigureCanvas(fig)
    ax = fig.add_subplot(111)
    #fi, ax = plt.subplots(1, 1)
    lsdf = lsdf_wNA.fillna(value=0)
    heatmap = ax.pcolormesh(wlsdf, cmap='Greens')
    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(wlsdf.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(wlsdf.shape[1]) + 0.5, minor=False)
    # # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    # # Set the labels
    ax.set_xticklabels([x.replace('_', '-') for x in wlsdf.columns.values],
                       minor=False, rotation=90)
    ax.set_yticklabels([x.replace('_', '-') for x in wlsdf.index],
                       minor=False)
    fig.colorbar(heatmap)  # add colorbar key
    fig.tight_layout()  # pad=0, w_pad=5, h_pad=.0)
    # fig.subplots_adjust(top=0.4, bottom=0.0, hspace=0, wspace=0.3 )
    fig.set_size_inches(8, 8)
    fig.savefig(str(output_prefix + "heatmap.png"), dpi=(300))
    fig.savefig(str(output_prefix + "heatmap.pdf"), dpi=(300))
    ####  plot clustered heatmap
#    plt.close('all')
#    plt.figure(1, figsize=(6, 6))
    # h = heatmapcluster(wlsdf.as_matrix(), row_labels=wlsdf.index,
    #                    col_labels=wlsdf.columns.values,
    #                    num_row_clusters=2, num_col_clusters=2,
    #                    label_fontsize=6,
    #                    xlabel_rotation=-75,
    #                    cmap=plt.cm.coolwarm,
    #                    show_colorbar=True,
    #                    top_dendrogram=True)
    # # plt.tight_layout()  # pad=0, w_pad=5, h_pad=.0)
    # # fig2.set_size_inches(16, 16)
    # plt.savefig(str(output_prefix + "clustered_heatmap.png"), dpi=(200))
    # plt.savefig(str(output_prefix + "clustered_heatmap.pdf"), dpi=(200))
    return lsdf_wNA


def make_msa(msa_tool, unaligned_seqs, prank_exe, mafft_exe,
             args, outdir, logger=None):
    """returns msa cmd and results path
    """
    assert logger is not None, "Must use logger"
    if msa_tool == "prank":
        if shutil.which(prank_exe):
            msa_cmd, results_path = prepare_prank_cmd(
                outdir=outdir,
                outfile_name="best_MSA",
                combined_fastas=unaligned_seqs,
                prank_exe=prank_exe,
                add_args=args,
                logger=logger)
        else:
            raise ValueError("Construction of MSA skipped because " +
                             "%s is not a valid executable!", prank_exe)
    elif msa_tool == "mafft":
        if shutil.which(mafft_exe):
            msa_cmd, results_path = prepare_mafft_cmd(
                outdir=outdir,
                outfile_name="best_MSA",
                combined_fastas=unaligned_seqs,
                mafft_exe=mafft_exe,
                add_args=args,
                logger=logger)
        else:
            raise ValueError("Construction of MSA skipped because " +
                             "%s is not a valid executable!", mafft_exe)
    return(msa_cmd, results_path)


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""
    smoothing algorithm
    http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
    Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k ** i for i in order_range] for k in
                range(- half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1: half_window + 1][:: -1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1: -1][:: -1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


def get_rec_from_generator(recordID, gen, method=None):
    """ given a record ID and and SeqIO generator return sequence of
    genbank record that has all the loci, and call method to refresh generator
    If on different sequences, return error
    """
    for record in gen:
        if recordID == record.id:
            if method is not None:
                method()
            return record
        else:
            pass
    # if none found, raise error
    raise ValueError("no record found matching record id %s!" % recordID)


def submain(clusters, gb_path, logger, verbose, no_revcomp,
         output, circular, flanking, prefix_name, args):
    get_rev_comp = no_revcomp is False  # kinda clunky
    flanking_regions_output = os.path.join(output, "flanking_regions_output")
    try:
        os.makedirs(flanking_regions_output)
    except FileExistsError:
        if args.clobber:
            pass
        else:
            sys.exit(1)
    extracted_regions = []
    logger.debug(clusters)
    for cluster in clusters:  # for each cluster of loci
        # get seq record that cluster is  from
        logger.debug("fetching genbank record for %s", cluster.sequence_id)
        try:
            genome_records_gen = SeqIO.parse(gb_path, 'genbank')
            cluster.seq_record = get_rec_from_generator(
                recordID=cluster.sequence_id,
                gen=genome_records_gen,
                method=None)
            logger.debug(cluster.seq_record)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        # make coord list
        try:
            extract_coords_from_locus(
                cluster=cluster, feature=cluster.feat_of_interest,
                logger=logger,
                verbose=False)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        logger.debug(
            "Here are the detected region,coords, strand, product, " +
            "locus tag, subfeatures and sequence id of the results:")
        logger.debug(str(cluster.__dict__))
        if circular:
            cluster_post_pad = pad_genbank_sequence(cluster=cluster,
                                                    logger=logger)
        else:
            cluster_post_pad = cluster

        #  given coords and a sequnce, extract the region as a SeqRecord
        try:
            cluster_post_stitch =\
                stitch_together_target_regions(cluster=cluster_post_pad,
                                               flanking=flanking,
                                               logger=logger,
                                               circular=circular,
                                               revcomp=get_rev_comp)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        # this is to allow sorting later
        cluster_post_stitch.extractedSeqRecord.name = \
            cluster_post_stitch.global_start_coord
        extracted_regions.append(cluster_post_stitch.extractedSeqRecord)
    # after each cluster has been extracted, write out results
    logger.debug(extracted_regions)
    # This keeps the records being produced in the same order each time
    sorted_regions = sorted(extracted_regions, key=lambda x: int(x.name),
                            reverse=False)
    file_list = []
    for index, region in enumerate(sorted_regions):
        logger.debug("index: %d", index)
        logger.debug("region: %s", region)
        if prefix_name is None:
            prefix_name = date
        # make dir for blastable partial seeds
        graded_output = os.path.join(
            output,
            "{0}_region_{1}_{2}.fasta".format(
                prefix_name, index + 1, "graded_flanking_regions"))
        ####  make this a cmdline option, maybe
        do_blast = True
        if do_blast:
            try:
                os.makedirs(graded_output)
            except FileExistsError:
                if args.clobber:
                    pass
                else:
                    sys.exit(1)
            for i in range(0, flanking + 100, 100):
                short_rec = SeqRecord(
                    #  needs the -1 to avoid the  index-by-negative-zero-error
                    Seq(str(region.seq[i: -i - 1]), IUPAC.IUPACAmbiguousDNA()),
                    id="{0}_{1}-bp_flanks".format(region.id, flanking - i))
                blast_fasta_filename = \
                    "{0}_region_{1}_{2}_{3}-bp_flanks.fasta".format(
                        prefix_name, index + 1, "riboSnag", flanking - i)
                with open(os.path.join(graded_output, blast_fasta_filename),
                          "w") as outfile_blast:
                    SeqIO.write(short_rec, outfile_blast, "fasta")
                file_list.append(
                    os.path.join(graded_output, blast_fasta_filename))
        ####
        filename = "{0}_region_{1}_{2}.fasta".format(
            prefix_name, index + 1, "riboSnag")
        filename_flanking = "{0}_region_{1}_{2}_flanking_regions.fasta".format(
            prefix_name, index + 1, "riboSnag")
        # TODO make discription work when writing seqrecord
        # TODO move date tag to fasta description?
        # i.description = str("{0}_riboSnag_{1}_flanking_{2}_within".format(
        #                           output_index, args.flanking, args.within))
        # output full region
        with open(os.path.join(output, filename), "w") as outfile:
            SeqIO.write(region, outfile, "fasta")
            outfile.write('\n')

        # Output 5' region and 3' region
        upstream_flanking = SeqRecord(
            Seq(str(region.seq[0: flanking]), IUPAC.IUPACAmbiguousDNA()),
            id="{0}_upstream_flanking_{1}:{2}".format(region.id, 0, flanking))
        downstream_flanking = SeqRecord(
            Seq(str(region.seq[-flanking: len(region.seq)]),
                IUPAC.IUPACAmbiguousDNA()),
            id="{0}_downstream_flanking_{1}:{2}".format(
                region.id, len(region.seq)-flanking, len(region.seq)))
        with open(os.path.join(flanking_regions_output,
                               filename_flanking), "w") as outfile:
            SeqIO.write(upstream_flanking, outfile, "fasta")
            outfile.write('\n')
            SeqIO.write(downstream_flanking, outfile, "fasta")
        # output 3' region


    # write out the whole file as a fasta as well...
    # append in case of multiple records
    ref_fasta = os.path.join(output, str(prefix_name + "_genome.fasta"))
    # refresh the generator
    genome_records_gen = SeqIO.parse(gb_path, 'genbank')
    with open(ref_fasta, "a") as outfa:
        for record in genome_records_gen:
            SeqIO.write(record, outfa, "fasta")
            outfa.write('\n')
    return extracted_regions, ref_fasta, file_list


def get_makeblastdb_cmd(input_file, input_type="fasta", dbtype="prot",
                        title="blastdb", out="blastdb",
                        makeblastdb_exe='', logger=None):
    """
    This runs make blast db with the given parameters
    requires logging, os, subprocess, shutil
    """
    if makeblastdb_exe == '':
        makeblastdb_exe = shutil.which("makeblastdb")
        if makeblastdb_exe is None:
            logger.error("error finding makeblastdb executable")
            raise ValueError("no executable for makeblastdb found!")
    makedbcmd = str("{0} -in {1} -input_type {2} -dbtype {3} " +
                    "-title {4} -out {5}").format(makeblastdb_exe,
                                                  input_file,
                                                  input_type,
                                                  dbtype, title, out)
    return makedbcmd


def make_blast_cmds(filename_list, blast_type, output, blastdb, date,
                    logger=None):
    """given a file, make a blast cmd, and return path to output csv
    """
    blast_cmds = []
    blast_outputs = []
    for f in filename_list:
        if logger:
            logger.debug("creating blast cmds for %s", f)
        output_path_tab = "{0}_{1}_results_{2}.tab".format(
            os.path.join(output, date),
            blast_type,
            os.path.basename(f))
        blast_cline = NcbiblastnCommandline(query=f,
                                            db=blastdb, evalue=10,
                                            outfmt=6, out=output_path_tab)
        if blast_type == 'blastn':
            add_params = str(
                " -num_threads 1 -max_target_seqs " +
                "2000 -task blastn -perc_identity 97")
        elif blast_type == 'dc_megablast':
            add_params = str(
                " -num_threads 1 -max_target_seqs " +
                "2000 -task dc_megablast -perc_identity 97")
        elif blast_type == 'tblastx':
            add_params = str(
                " -num_threads 1 -max_target_seqs 2000 " +
                "-query_gencode 11 -db_gencode 11 -perc_identity 97")
        else:
            raise ValueError("must use either blastn or tblastx")
        blast_command = str(str(blast_cline) + add_params)
        blast_cmds.append(blast_command)
        blast_outputs.append(output_path_tab)
    return(blast_cmds, blast_outputs)


def merge_outfiles(files, outfile_name, logger=None):
    """
    """
    # only grab .tab files, ie, the blast output
    logger.debug("files to merge:")
    logger.debug("\n".join(files))
    filelist = [i for i in files if "tab" in os.path.splitext(i)[-1]]
    if len(filelist) == 0:
        if logger:
            logger.error("No BLAST output files found!")
        raise FileNotFoundError

    elif len(filelist) == 1:
        if logger:
            logger.warning("only one file found! no merging needed")
        return(filelist)
    else:
        nfiles = len(filelist)
        if logger:
            logger.info("merging %i blast results to %s", nfiles, outfile_name)
        fout = open(outfile_name, "a")
        # first file:
        for line in open(filelist[0]):
            fout.write(line)
        #  now the rest, ignoring the header:
        for num in range(1, nfiles):
            f = open(filelist[num])
            for line in f:
                fout.write(line)
            f.close()  # not really needed
        fout.close()
    return(outfile_name)


def run_blast(query_list, ref, name, output, mbdb_exe='', logger=None):
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    logger.info("constructing makeblastdb command")
    db_name = os.path.join(output, str(name + "_db"), name)
    mbdb_cmd = get_makeblastdb_cmd(
        input_file=ref, input_type="fasta", dbtype="nucl",
        title=name, out=db_name,
        makeblastdb_exe=mbdb_exe, logger=logger)
    logger.debug(mbdb_cmd)
    logger.info("creating blast cmds")
    blast_outdir = os.path.join(output, "BLAST_results")
    os.makedirs(blast_outdir)
    blast_cmds, paths_to_outputs = make_blast_cmds(
        filename_list=query_list, blast_type="blastn",
        output=blast_outdir, blastdb=db_name, date=date, logger=logger)
    logger.info("Making BLAST Database")
    subprocess.run(mbdb_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    # pool = multiprocessing.Pool(processes=args.cores)
    pool = multiprocessing.Pool()
    logger.info("Running BLAST commands")
    logger.debug("Running the following commands in parallel " +
                 "(this could take a while):")
    logger.debug("\n" + "\n".join([x for x in blast_cmds]))
    results = [
        pool.apply_async(subprocess.run,
                         (cmd,),
                         {"shell": sys.platform != "win32",
                          "stdout": subprocess.PIPE,
                          "stderr": subprocess.PIPE,
                          "check": True})
        for cmd in blast_cmds]
    pool.close()
    pool.join()
    reslist = []
    reslist.append([r.get() for r in results])
    logger.debug("merging resulting blast files")
    merged_tsv = merge_outfiles(
        paths_to_outputs,
        os.path.join(blast_outdir,
                     str(date + "_results_merged.csv")),
        logger=logger)
    logger.info("Merged BLAST results can be found here for perusal or plotting:")
    logger.info(merged_tsv)


def main(args, logger=None):
    output_root = os.path.abspath(os.path.expanduser(args.output))
    # Create output directory only if it does not exist
    try:
        os.makedirs(args.output)
    except FileExistsError:
        # '#' is needed in case streaming output eventually
        print("#Selected output directory %s exists" %
              args.output)
        if not args.clobber:
            print("exiting")
            sys.exit(1)
        else:
            print("# continuing, and risking potential loss of data")
    log_path = os.path.join(output_root, "riboSnag.log")
    if logger is None:
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
    # test_ex = [check_installed_tools(x, logger=logger) for x in executables]
    if shutil.which(args.barrnap_exe):
        logger.debug("All needed system executables found!")
        logger.debug(str([shutil.which(i) for i in executables]))
    check_version_from_cmd(
        exe=args.barrnap_exe,
        cmd='', line=2,
        pattern=r"barrnap (?P<version>[^-]+)",
        where='stderr',
        logger=logger,
        coerce_two_digit=True,
        min_version="0.0.7")
    # parse cluster file
    check_loci_file_not_genbank(args.clustered_loci, )
    try:
        clusters = parse_clustered_loci_file(args.clustered_loci,
                                             gb_filepath=args.genbank_genome,
                                             output_root='',
                                             padding=args.padding,
                                             circular=not args.linear,
                                             logger=logger)
        clusters = add_gb_seqrecords_to_cluster_list(
            cluster_list=clusters,
            gb_filepath=args.genbank_genome)
    except Exception as e:
        logger.error(e)
        sys.exit(1)
    # parse genbank records
    # with open(args.genbank_genome) as fh:
    #     genome_records = list(SeqIO.parse(fh, 'genbank'))
    # genome_records_gen = SeqIO.parse(args.genbank_genome, 'genbank')
    regions = []
    logger.debug("clusters:")
    for cluster in clusters:
            logger.debug(cluster.__dict__)
    if args.name is None:  # if none given, use date for name (for out files)
        args.name = date
    regions, ref_fasta, region_files = submain(
        clusters=clusters,
        gb_path=args.genbank_genome,
        logger=logger,
        verbose=False,
        flanking=args.flanking,
        output=output_root,
        circular= not args.linear,
        prefix_name=args.name,
        no_revcomp=args.no_revcomp,
        args=args
    )

    # make MSA and calculate entropy
    # reminder, PLOT is set by matplotlib avail
    if not args.just_extract and PLOT:
        if args.clobber:
            logger.error("Cannot safely check MSA when --clobber is used!")
            sys.exit(1)

        unaligned_seqs = combine_contigs(contigs_dir=output_root,
                                         pattern="*riboSnag*",
                                         contigs_name="riboSnag_unaligned",
                                         ext=".fasta", verbose=False,
                                         logger=logger)
        msa_cmd, results_path = make_msa(msa_tool=args.msa_tool,
                                         unaligned_seqs=unaligned_seqs,
                                         prank_exe=args.prank_exe,
                                         args='',
                                         mafft_exe=args.mafft_exe,
                                         outdir=output_root,
                                         logger=logger)

        logger.info("Running %s for MSA", args.msa_tool)
        subprocess.run(msa_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)

        seq_entropy, names, tseq = calc_entropy_msa(results_path)
        if not args.skip_kmers:
            if args.msa_kmers:
                with open(results_path, 'r') as resfile:
                    kmer_seqs = list(SeqIO.parse(resfile, "fasta"))
            else:
                kmer_seqs = regions
            counts, names = profile_kmer_occurances(kmer_seqs,
                                                    # alph='atcg-',
                                                    k=5,
                                                    logger=logger)

            mca_df = plot_pairwise_least_squares(
                counts=counts, names_list=names,
                output_prefix=os.path.join(
                    output_root,
                    "sum_least_squares"))
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
        if args.title is None:
            label_name = args.seq_name if args.seq_name is not None else \
                         os.path.basename(
                             os.path.splitext(
                                 args.genbank_genome)[0])
        else:
            label_name = args.title
        return_code = plot_scatter_with_anno(
            data=seq_entropy,
            consensus_cov=consensus_cov,
            names=["Position", "Entropy"],
            title=args.title,
            anno_list=annos,
            pubplot=args.pubplot,
            output_prefix=os.path.join(
                output_root,
                "entropy_plot"))
        # plot_alignment_3d(
        #     output_prefix=os.path.join(output_root, "entropy_plot"),
        #     consensus=consensus_cov,
        #     tseq=tseq)
        if not args.skip_blast:
            run_blast(query_list=region_files, ref=ref_fasta,
                      mbdb_exe=args.makeblastdb_exe, name=args.name,
                      output=args.output, logger=logger)

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
    matplotlib.rc('font', family='sans-serif')
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = 'cm'
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

from pyutilsnrw.utils3_5 import check_installed_tools,\
    combine_contigs, check_version_from_cmd

from .classes import LociCluster, Locus

from .shared_methods import add_gb_seqrecords_to_cluster_list, \
    extract_coords_from_locus, pad_genbank_sequence, \
    parse_clustered_loci_file, set_up_logging

def get_args():  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(prog="ribo structure",
                                     description="Given a list of genomes, " +
                                     "plot the organizational structure  " +
                                     "of the ribosomal coding regions",
                                     add_help=False)
    parser.add_argument("dir", help="folder with (multi)fasta genomes.")
    args = parser.parse_args(sys.argv[2:])
    return args



def plot_scatter_with_anno(data,
                           consensus_cov,
                           anno_list,
                           names=["Position", "Entropy"],
                           title="$Strain Name$",
                           output_prefix="entropy_plot.png"):
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
    ax1.set_title("Shannon Entropy by Position\n" + title, y=1.08)
    xmin, xmax = 0, len(data)
    ymin, ymax = -0.1, (max(data) * 1.2)
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
    yjust = -.1
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
        ax1.text((anno[1][0] + anno[1][1]) / 2,    # x location
                 ymax - 0.48 - yjust,                      # y location
                 anno[0][0:20],                          # text first 20 char
                 ha='center', color='red', weight='bold', fontsize=10)
        yjust = yjust * - 1
    ax1.scatter(x=df["Position"], y=df["Entropy"],
                marker='o', color='black', s=2)
    # add smoothing for kicks
    df["fit"] = savitzky_golay(df["Entropy"].values, 351, 3)  # window size 51, polynomial order 3
    ax1.scatter(x=df["Position"], y=df["fit"], color='red', s=1)
    #
    ax1.set_ylabel('Shannon Entropy')
    ax1.get_yaxis().set_label_coords(-.05, 0.5)
    ax2.set_xlim([xmin, xmax])
    ax2.invert_yaxis()
    ax2.set_ylabel('Consensus Coverage')
    ax2.set_xlabel('Position (bp)')
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
    ax1.tick_params(axis='y', colors='dimgrey')
    ax2.tick_params(axis='y', colors='dimgrey')
    ax1.tick_params(axis='x', colors='dimgrey')
    ax2.tick_params(axis='x', colors='dimgrey')
    ax1.yaxis.label.set_color('black')
    ax2.yaxis.label.set_color('black')
    ax1.xaxis.label.set_color('black')
    ax2.xaxis.label.set_color('black')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.set_size_inches(12, 7.5)
    fig.savefig(str(output_prefix + '.png'), dpi=(200))
    fig.savefig(str(output_prefix + '.pdf'), dpi=(200))
    return 0


class RRNA(object):
    def __init__(self, record_name=None, name=None, product=None, strand=None,
                 start=None):
        self.record_name = record_name
        self.name = name
        self.product = product
        self.strand = strand
        self.start = start
    def __str__(self):
        return str("rRNA:\t{seq}\t{name}\t{product}\t{strand}\t{start}".format(
            seq=self.record_name,
            name=self.name,
            product=self.product,
            strand=self.strand,
            start=self.start))


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

def get_duplicates(a):
    seen = set()
    uniq = []
    for x in a:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    return seen

def find_pairs(list_16, list_23, max_dist):
    pairs = []
    orphans = []
    for st in list_16:
        bestpair_diff = None
        for tw in list_23:
            diff = abs(st.start - tw.start)
            if diff < max_dist:
                # if this distance if smaller than our previous best, replace
                if bestpair_diff is None or  diff < bestpair_diff[1]:
                    best_pair = [(st, tw), diff]
        pairs.append(bestpair_diff[0])
    if len(pairs) > 1:
        # if we have duplicates, see if there are any that we can shift from one to
        duplicated_23_pairs = [(x, y) for x, y in pairs if
                               y in get_duplicates([w for w, z in pairs])]

        duplicated_23_pairs = get_duplicates([y for x, y in pairs])
        for l in [duplicated_16, duplicated_23]:



    duplicated_23

    from glob import glob

    fastas = glob(args.dir + os.path.sep + "*.fa*")
    for genome in fastas:
        # run barrnap
        # get list of RRNAs
        # get pair list and orphan list
        refine_coord list

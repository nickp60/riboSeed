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
from glob import glob
import shutil
import itertools
import multiprocessing

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import \
        FigureCanvasAgg as FigureCanvas
    from matplotlib.patches import FancyBboxPatch
    from matplotlib.figure import Figure
    from matplotlib import gridspec
    import matplotlib.patches as patches
    matplotlib.rc('font', family='sans-serif')
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

# import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from collections import defaultdict  # for calculating kmer frequency
from itertools import product  # for getting all possible kmers
# from heatmapcluster import heatmapcluster

from .shared_methods import make_barrnap_cmd, set_up_logging



def get_args(test_args=None):  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(prog="ribo structure",
                                     description="Given a list of genomes, " +
                                     "plot the organizational structure  " +
                                     "of the ribosomal coding regions",
                                     add_help=False)
    parser.prog = "ribo structure"
    parser.add_argument("dir", help="folder with (multi)fasta genomes.")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-g", "--gffdir", dest='gffdir',
                          help="If used, look here for gff files " +
                          "corresponding to the sequences rather than " +
                          "generating with barrnap",
                          action="store", type=str)
    optional.add_argument("-k", "--kingdom", dest='kingdom',
                          action="store",
                          choices=["bac", "euk", "arc", "mito"],
                          help="whether to look for eukaryotic, archaeal, or" +
                          " bacterial rDNA; " +
                          "default: %(default)s", default="bac",
                          type=str)
    optional.add_argument("-t", "--id_thresh", dest='id_thresh',
                          action="store", type=float,
                          help="partial rRNA hits below this threshold will " +
                          "be ignored. " +
                          "default: %(default)s", default=0.5,
                          )
    optional.add_argument("-b", "--barrnap_exe", dest='barrnap_exe',
                          action="store",
                          help="path to barrnap executable; " +
                          "default: %(default)s", default="barrnap",
                          type=str)
    optional.add_argument("-c", "--cores", dest='cores',
                          action="store",
                          help="number of threads/cores to use; " +
                          "default: %(default)s", default=2,
                          choices=[1, 2, 4, 8, 16], type=int)
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    args = parser.parse_args(sys.argv[2:])
    return args

mycolors = {  # pragma: no cover
    "pinkish": matplotlib.colors.ColorConverter().to_rgba(
        "#ff4c05", alpha=1),
    "redish": matplotlib.colors.ColorConverter().to_rgba(
        "#ff4c05", alpha=1),
    "yellish": matplotlib.colors.ColorConverter().to_rgba(
        "#FFFB07", alpha=1),
    "greenish": matplotlib.colors.ColorConverter().to_rgba(
        "#00D424", alpha=1),
    "bluish": matplotlib.colors.ColorConverter().to_rgba(
        "#06B9FF", alpha=1),
    "greyish": matplotlib.colors.ColorConverter().to_rgba(
        "#7E7F97", alpha=1),
    "lightergreyish": matplotlib.colors.ColorConverter().to_rgba(
        "#9A9BB8", alpha=1),
    "clear": matplotlib.colors.ColorConverter().to_rgba(
        "#FF012F", alpha=0),
}


bgcols = {  # pragma: no cover
    "purle": matplotlib.colors.ColorConverter().to_rgba(
        "#EB87A3", alpha=0.5),
    "green": matplotlib.colors.ColorConverter().to_rgba(
        "#5EA662", alpha=0.5),
    "yellow": matplotlib.colors.ColorConverter().to_rgba(
        "#EBE418", alpha=0.5),
    "red": matplotlib.colors.ColorConverter().to_rgba(
        "#EB7D7D", alpha=0.5),
    "blue": matplotlib.colors.ColorConverter().to_rgba(
        "#6795A6", alpha=0.5),
}


def plot_rDNAs(
        gff_lists,
        featuremin,
        maxlen,
        breakwidth=40,
        # names=["Position", "Entropy"],
        title="Shannon Entropy by Position",
        output_prefix="entropy_plot.png", logger=None):
    """
    plot rDNA positions from a list of genomes

    given a list of gff-like structures (one per genome),
    and a list of records for each genome,
    Add a line for each.

    """
    logger.debug("sources for datasets for plotting")
    logger.debug([x[0][0] for x in gff_lists])
    max_combined_len = 10000 + maxlen
    fig = Figure()
    FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])#, sharey=ax)
    ax2.get_shared_y_axes().join(ax, ax2)
    ax.set_title(title, y=1.08)
    # set the centers as starting relative to  relheight - (2* codingdepth)
    # and out relative height: assuming max 10Mb genome, 1/10 should be nice
    ngffs = len(gff_lists)
    coding_height = 1.35
    rRNA_height = coding_height / 3
    print("coding height " + str(coding_height))

    # define our center lines of the graphs
    centers = [x + .5  for x in range(0, (len(gff_lists) * 2) + 1, 2)]
    centers[-1] = centers[-1] - 1  # bit of a hack to tidy up the spacing

    logger.info(centers)
    xmin, xmax = 0, max_combined_len
    ymin, ymax = 0, ngffs #* len(gff_lists) #* 1.5
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    # ax.grid(color='r', linestyle='-', linewidth=2)
    ax2.set_ylim([ymin, ymax])
    #  plot the color shadings
    unused_cols = ["red", "green", "yellow", "purple", "red", "blue"]
    nudge = coding_height / 2
    patch_list = []
    # add annotations
    print(len(gff_lists))
    names, descs, short_descs, lens = [], [], [], []
    for gff_i, gff in enumerate(gff_lists):
        logger.info("plotting %i of %i", gff_i + 1, len(gff_lists))
        # if i > 2:
        #     continue
        with open(gff[0][0], "r") as inf:
            recs = list(SeqIO.parse(inf, "fasta"))
        last_chrom_end = 0
        names.append(",".join([x.id for x in recs]))
        descs.append(",".join([x.description for x in recs]))
        short_name = recs[0].description
        if len(recs) > 1:
            short_name = short_name + "; " + ",".join([x.id for x in recs[1:]])
        short_descs.append(short_name)
        # lens.append(",".join([str(len(x.seq)) for x in recs]))
        lens.append(sum([len(x.seq) for x in recs]))
        coding_lengths = [len(x.seq) for x in recs]
        coding_box = patches.Rectangle(
            (last_chrom_end, # x1
             centers[gff_i] - (coding_height *.5)), # y1
            sum(coding_lengths), # x2
            abs(coding_height), # y2
            fc=mycolors['lightergreyish'],
            ec=mycolors['clear'],
            linewidth=0
        )

        ax.add_patch(coding_box)
        # add chrom breaks
        last_chrom_end = last_chrom_end + coding_lengths[0]
        if len(coding_lengths) > 1:
            for length_i, length in enumerate(coding_lengths):
                if length_i == 0:
                    continue
                buffer_box = patches.Rectangle(
                    (last_chrom_end - breakwidth,
                     centers[gff_i] - coding_height * .6),
                    breakwidth,
                    coding_height * 1.2,
                    fc="black",
                    ec=mycolors['clear']
                )
                last_chrom_end = last_chrom_end + length
                ax.add_patch(buffer_box)
        # add tracking ticks to label
        diff_this_max =  maxlen - sum(coding_lengths)
        if diff_this_max > 0:
            ax.plot([last_chrom_end + breakwidth, maxlen],
                    [centers[gff_i],centers[gff_i]], color = 'black', linestyle='--', dashes=(2,5))

        previous_end = 0
        for f, seq, prog, feat, start, end, x1, strand, x2, product in gff:
            start, end = int(start), int(end)
            if feat != "rRNA" and i == 0:
                #Exclude this feature
                continue
            feat_len = int(end) - int(start)
            # ------------------------ .5 total
            # |                      |
            # ------------------------
            # |                      |
            # ------------------------
            # |                      |
            # ------------------------   0
            if "16S" in product:
                color = mycolors['redish']
                lencode = featuremin * 40
                yjust = (2 * coding_height)/3
            elif "23S" in product:
                color = mycolors['yellish']
                lencode = featuremin * 60
                yjust = coding_height/3
            elif "5S" in product:
                color = mycolors['greenish']
                lencode = featuremin * 20
                yjust = 0
            else:
                continue
            anno_box = patches.Rectangle(
                (start, # x1
                 centers[gff_i] - (coding_height * .5 ) + yjust ),  # y1
                lencode, # x2
                abs(rRNA_height),
                fc=color,
                ec="black",
                alpha=.75,
                linewidth=.5
            )
            ax.add_patch(anno_box)
            # previous_end = this_end
    # add a legend to this mess
    red_patch = patches.Patch(color=mycolors['redish'], label='16S')
    yellish_patch = patches.Patch(color=mycolors['yellish'], label='23S')
    green_patch = patches.Patch(color=mycolors['greenish'], label='5S')

    ax.legend(bbox_to_anchor=(.1, -0.1, .9, -.102), loc=3,
              ncol=3, mode="expand", borderaxespad=0.,
              handles=[red_patch, yellish_patch, green_patch],
              frameon=False)
    # # for righthad labels
    loc = -.1
    for desc_i, desc in enumerate(short_descs):
        ax2.text(loc,    # x location
                 centers[desc_i] - (nudge * .4),                      # y location
                 desc,                          # text first 20 char
                 ha='left', color='black', weight='normal', fontsize=14)
    # label right with lenth
    ax.set_yticklabels(lens, fontsize=14)
    # ax.set_yticklabels(centers, fontsize=14)
    ax.set_ylabel('Total Genome Length (including plasmids)', fontsize=14)
    ax.get_yaxis().set_label_coords(-.15, .5) # left/right, up/down
    ax.yaxis.label.set_color('dimgrey')
    # ax.get_yaxis().set_label_coords(.05, .1)
    # remove all labels for axis 2
    ax2.get_yaxis().set_visible(False)
    ax2.get_xaxis().set_visible(False)

    #set the ticks/spines
    ax.spines['top'].set_visible(True)
    ax2.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('top')
    ax2.xaxis.set_ticks_position('none')
    # ax.tick_params(axis='y', colors='dimgrey')
    ax.tick_params(axis='x', colors='dimgrey')
    ax.xaxis.get_major_formatter().set_scientific(False)
    # set axis labels
    ax.xaxis.label.set_color('black')

    # remove all the other spines
    for x in [ax, ax2]:
        x.set_yticks(np.array(centers))
        x.spines["left"].set_visible(False)
        x.spines["right"].set_visible(False)
        x.spines["bottom"].set_visible(False)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.set_size_inches(22, (.75* len(names)))
    logger.info("writing plots to file")
    # fig.savefig(str(output_prefix + '.png'), bbox_inches='tight', dpi=(200), transparent=True)
    fig.savefig(str(output_prefix + '.pdf'), bbox_inches='tight', dpi=(200))
    return 0





class RRNA(object):
    newid = itertools.count()

    def __init__(self, record_name=None, kind=None, product=None, strand=None,
                 start=None):
        self.name = next(RRNA.newid)
        self.kind = kind
        self.record_name = record_name
        self.strand = strand
        self.start = start
    def __str__(self):
        return str("rRNA:\t{seq}\t{name}\t{kind}\t{strand}\t{start}".format(
            seq=self.record_name,
            kind=self.kind,
            name=self.name,
            strand=self.strand,
            start=self.start))
    def __repr__(self):
        return str(self.name)



def get_duplicates(a):
    seen = set()
    uniq = []
    for x in a:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    return seen



def find_pairs(list_16, list_23, max_dist):
    each_possibilities = []
    for _16s in list_16:
        print(_16s.name)
        # make an empt list for each of the possible neighbors
        each_possibilities.append((_16s.name, []))
        for _23s in list_23:
            print(_23s)
            diff = abs(_16s.start - _23s.start)
            print(diff)
            if diff < max_dist:
                # add it to the sublist in the most recent item
                each_possibilities[-1][1].append(_23s)
    print(each_possibilities)
    n_rDNAs = len(each_possibilities)
    print(n_rDNAs)
    final = []
    unavailable_23s = []
    still_needs_reducing = True
    max_i = 50
    i = 0
    while len(final) != n_rDNAs:
        i = i + 1
        print("hshgv")
        for _16s, poss in each_possibilities:
            poss = [x for x in poss if x not in unavailable_23s]
            if len(poss) == 0:
                final.append((_16s, None))
                each_possibilities = [(x, y) for x, y in each_possibilities if x != _16s]
            elif len(poss) == 1:
                final.append((_16s, poss[0]))
                unavailable_23s.append(poss[0])
                each_possibilities = [(x, y) for x, y in each_possibilities if x != _16s]
            else:
                all_single = False
        if i > max_i:
            raise ValueError("Unable to reduce neighbors!")
    return final


def main(args, logger=None):
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
    log_path = os.path.join(output_root, "riboStructure.log")
    if logger is None:
        logger = set_up_logging(verbosity=args.verbosity,
                                outfile=log_path,
                                name=__name__)
    logger.debug("Usage:\n{0}\n".format(str(" ".join([x for x in sys.argv]))))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, v)

    new_gffs = []
    source_list = []
    plot_gff_lists = []
    input_fastas = sorted(glob(args.dir + os.path.sep + "*.f*a*"))
    if args.gffdir is not None:
        new_gffs = sorted(glob(args.gffdir + os.path.sep + "*.gff"))
        if len(new_gffs) != len(input_fastas):
            raise ValueError("number of fastas must match gffs!")
        for i, gff in enumerate(new_gffs):
            if os.path.basename(os.path.splitext(gff)[0]) != \
               os.path.basename(os.path.splitext(input_fastas[i])[0]):
                raise ValueError("fasta basnames must match gffs!")
    else: # run with barrnap
        for f in input_fastas:
            source = os.path.basename(os.path.splitext(f)[0])
            barrnap_cmd = make_barrnap_cmd(
                infasta=f,
                outgff=os.path.join(output_root, "{0}.gff".format(source)),
                exe=args.barrnap_exe,
                threads=args.cores,
                thresh=args.id_thresh,
                kingdom=args.kingdom)
            logger.info("running barrnap cmd: %s", barrnap_cmd)
            subprocess.run(barrnap_cmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
            new_gff = os.path.join(output_root, "{0}.gff".format(source))
            new_gffs.append(new_gff)
    for ix, f in enumerate(input_fastas):
    # for i, gfffile in enumerate(new_gffs):
        gff_list = []
        with open(new_gffs[ix], 'r') as g:
            for idx, line in enumerate(g):
                if idx == 0:
                    #gets rid of header line
                    pass
                else:
                    l = [f] + line.strip().split("\t")
                    gff_list.append(l)
        rRNA_list = []
        for kind in ["16S", "23S", "5S"]:
            for line in gff_list:
                if kind in line[9]:
                    rRNA_list.append(
                        RRNA(record_name=line[1],
                             kind=kind,
                             start=int(line[4]),
                             strand=line[7]))
        # _16s_list = [x for x in rRNA_list if x.kind == "16S"]
        # _23s_list = [x for x in rRNA_list if x.kind == "23S"]
        # print(find_pairs(list_16=_16s_list, list_23=_23s_list, max_dist=17000))
        # for bit in gff_list:
        #     print(bit)
        if len(gff_list) > 0:
            plot_gff_lists.append(gff_list)
            source_list.append(f)

    #####################################################################
    print("plotting %i strains" % len(plot_gff_lists))
    sequence_lens = []
    for source in source_list:
        with open(source, "r") as inf:
            sequence_lens.append(sum(
                [len(x.seq) for x in SeqIO.parse(inf, "fasta")]))
    logger.debug("max: %s", max(sequence_lens))
    # zip and sort by the sequence length, low to high
    sorted_plot_gff_lists = [x for y, x in sorted(zip(sequence_lens, plot_gff_lists))]
    # for x, y in sequence_lens, plot_mauve_compare
    plot_rDNAs(
        sorted_plot_gff_lists,
        breakwidth=4000,
        featuremin=1000,
        maxlen=max(sequence_lens),
        # names=["Position", "Entropy"],
        title="Relative rRNA coding sequences",
        output_prefix=os.path.join(output_root, "rDNA_relative_locations"),
        logger=logger)

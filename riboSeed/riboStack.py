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
import glob
import re
import math
import statistics as stats
from bisect import bisect
import shutil
import subprocess
import logging
import os
import traceback

# need this line for unittesting
sys.path.append(os.path.join('..', 'riboSeed'))
from .shared_methods import set_up_logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# --------------------------- classes --------------------------- #


# --------------------------- methods --------------------------- #


def n_samples_type(x):
    x = int(x)
    if x < 2:
        raise argparse.ArgumentTypeError("Minimum number of samples is 2")
    return x

def get_args(test_args=None):  # pragma: no cover

    parser = argparse.ArgumentParser(
        description="This facilitates the mapping of reads to a reference " +
        "and comparison of coverage depths in rDNA regions to assess " +
        "disparity in rDNA counts between the reference and your reads",
        add_help=False)
    parser.prog = "ribo stack"
    parser.add_argument("riboScan_dir", action="store",
                        help="We need the gff and fasta files from your " +
                        "riboScan run.")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-b", "--bam", dest='bam', action="store",
                          help="BAM file; tested with BWA output; " +
                          "default: %(default)s",
                          type=str)
    optional.add_argument("-r", "--riboSeed_dir", dest='riboSeed_dir',
                          action="store",
                          help="look for BAM file in this riboSeed output " +
                          "directory ",
                          required="-b" not in sys.argv, # require if no BAM
                          type=str)
    optional.add_argument("-n", "--n_samples", dest='n_samples',
                          action="store",
                          help="Number of regions to compare rDNA depth to; " +
                          "must be greater than 1; " +
                          "default: %(default)s",
                          type=n_samples_type, default=10)
    optional.add_argument("-e", "--experiment_name", dest='experiment_name',
                          action="store",
                          help="prefix for results files; " +
                          "default: %(default)s",
                          default="riboStack", type=str)
    optional.add_argument("-i", "--infer", dest='infer',
                          action="store_true",
                          help="If --infer, ignore the name and length " +
                          "of reference sequence, using the BAM " +
                          "header instead " +
                          "default: %(default)s", default=False)
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    if test_args is None:
        args = parser.parse_args(sys.argv[2:])
    else:
        args = parser.parse_args(test_args)
    return args


def last_exception():
    """ Returns last exception as a string, or use in logging.
    stolen verbatim from pyani
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


def makeRegions(outdir, gff, dest, name="", logger=None):
    coord_list = []
    with open(gff, "r") as gf:
        for idx, line in enumerate(gf):
            if idx is not 0:
                fields = line.strip().split("\t")
                coord_list.append([fields[0], fields[3], fields[4]])
    logger.debug("coord list:")
    logger.debug(coord_list)

    with open(dest, "w") as out:
        for coords in coord_list:
            line = "\t".join(coords) if name is "" else \
                   "{0}\t{1}\t{2}".format(
                       name, coords[1], coords[2])
            out.write(line + "\n")

    # # skip first line, get columns with chrom, start, end, separated by tabs
    # cmd = "awk '{if (NR!=1) {print $1,$4,$5}}' OFS='\\t' " + "{0} > {1}".format(
    #     gff,
    #     dest)
    # if logger:
    #     logger.debug("cmd to extract regions %s", cmd)
    # subprocess.run(cmd,
    #                shell=sys.platform != "win32",
    #                stdout=subprocess.PIPE,
    #                stderr=subprocess.PIPE,
    #                check=True)

def get_bam_from_riboSeed(d, iteration=0, logger=None):
    # I have both mapping_for_iteration and mapping_iteration
    # because I hate myself, apparently...
    if not os.path.isdir(d):
        logger.error("riboSeed directory not found: %s", d)
        raise(FileNotFoundError)
    pattern = os.path.join(
        d,
        "*_mapping_for_iteration_{}".format(iteration),
        "*_mapping_iteration_{}_sorted.bam".format(iteration))
    options = glob.glob(pattern)
    if len(options) == 0:
        logger.error("No BAM file found in riboSeed directory with glob %s",
                     pattern)
        raise(FileNotFoundError)
    if len(options) > 1:
        logger.warning("multiple bam files found; using first: %s",
                       options[0])
    return options[0]


def makeBedtoolsShuffleCmd(region, destdir, genome, bedtools_exe, n=10):
    cmd_list = []
    results_list = []
    for i in range(0, n):
        bedcmd = "{0} shuffle -i {1} -g {2} > {3}_{4}".format(
            bedtools_exe, region, genome,
            os.path.join(destdir, "sample_region"), i + 1)
        cmd_list.append(bedcmd)
        results_list.append(
            os.path.join(destdir, "sample_region_" + str(i + 1)))
    return(cmd_list, results_list)


def samtoolsGetDepths(samtools_exe, bam, ref_reg_file,
                      sample_file_list, outdir):
    """ return list [ref_file, [list, of, sample, files]]
    """
    ref_out = os.path.join(outdir, "ref_out_depth")
    cmds = []
    ref_cmd = "{0} view -b -F 256 {1} | {0} depth -b {2} - > {3}".format(
        samtools_exe,
        bam,
        ref_reg_file,
        ref_out)
    cmds.append(ref_cmd)
    sample_outs = []
    for idx, f in enumerate(sample_file_list):
        sample_out = os.path.join(outdir, "sample_out_depth_" + str(idx + 1))
        cmd = "{0} view -b -F 256 {1} | {0} depth -b {2} - > {3}".format(
            samtools_exe,
            bam,
            f,
            sample_out)
        cmds.append(cmd)
        sample_outs.append(sample_out)
    return [ref_out, sample_outs, cmds]


def getRecLengths(fasta, name=""):
    """ returns a list of [rec.id, record lengths]
    """
    with open(fasta, "r") as fa:
        flist = [[x.id, len(x)] for x in list(SeqIO.parse(fa, "fasta"))]
    if name is not "":
        flist[0][0] = name
    return flist


def printPlot(data, line=None, ymax=30, xmax=60, tick=.2,
              title="test", fill=False, logger=None):
    """ ascii not what your program can do for you...
    """
    assert logger is not None, "must use logging"
    data = sorted(data, reverse=True)
    xaxis = "|"
    yaxis = "_"
    avbin = []
    scaledy = []
    # ylab = str(max(data))
    ylab = str(len(data))
    sli = math.ceil(len(data) / ymax)
    if sli == 0 or len(data) == 0:
        return 1
    #  get rough averages for a window
    for i in range(0, ymax + 1):
        avbin.append(int(sum(data[i * sli: (i + 1) * sli]) / sli))

    avmax = max(avbin)
    logger.debug("scaling to max of %i", avmax)
    for j in avbin:
        scaledy.append(int((j / avmax) * xmax))
    if line is not None:
        scaled_line = int((line / avmax) * xmax)
        lineidx = len(scaledy) - bisect(sorted(scaledy), scaled_line)
    plotlines = []
    fillchar = " " if not fill else "X"
    for idx, j in enumerate(scaledy):
        if idx == 0:
            plotlines.append(" " * int(xmax * .25) + title)
            plotlines.append(" " * 9 + "0" + " " * (xmax - 2) + ylab)
            plotlines.append(" " * 10 + xaxis + (yaxis * xmax) + "|")
        if line is not None:
            if lineidx == idx:
                plotlines.append(" " * 10 + xaxis +
                                 "*" * (xmax - 4) + "  " + str(line))
                fillchar = " "
        if idx % int(tick * ymax) == 0:
            # plot ticks at increments
            ticlab = str(int((j / xmax) * avmax))
            plotlines.append(ticlab.rjust(10, " ") + xaxis + fillchar *
                             j + "O")
        else:
            plotlines.append(" " * 10 + xaxis + fillchar * j + "O")
    logger.info("\n" + "\n".join(plotlines))
    return(plotlines)


def main(args, logger=None):
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory already exists; exiting...")
        sys.exit(1)
    t0 = time.time()
    log_path = os.path.join(output_root, "riboStack.log")

    if logger is None:
        logger = set_up_logging(verbosity=args.verbosity,
                                outfile=log_path,
                                name=__name__)
    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))

    # check executables
    bedtools_exe = shutil.which("bedtools")
    samtools_exe = shutil.which("samtools")
    for exe in [samtools_exe, bedtools_exe]:
        if not exe:
            logger.error("must have bedtools and samtools installed " +
                         "and available in $PATH")
            sys.exit(1)

    gff = os.path.join(args.riboScan_dir, "scannedScaffolds.gff")
    fasta = os.path.join(args.riboScan_dir, "scannedScaffolds.fa")
    logger.info("GFF file: %s", gff)
    logger.info("fasta file: %s", fasta)
    if args.riboSeed_dir is not None:
        if args.bam is not None:
            logger.warning("Both BAM and riboSeed dir were provided; " +
                           "using bam from riboSeed dir")
        args.bam = get_bam_from_riboSeed(
            args.riboSeed_dir, iteration=0, logger=logger)
    if not os.path.exists(args.bam):
        logger.error("BAM file not found: %s", args.bam)
        sys.exit(1)
    for f in [gff, fasta]:
        if not os.path.exists(f):
            logger.error("file %s not found!", f)
            sys.exit(1)
    logger.info("making region list from gff file")
    rDNA_regions = os.path.join(output_root, "rDNA_regions")
    name = ""
    if args.infer:
        cmd = "{0} view -H {1} |grep '^@SQ'".format(
            samtools_exe, args.bam)
        res = subprocess.run(cmd,
                             shell=sys.platform != "win32",
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             check=True)
        res = res.stdout.decode("utf-8").split("\n")[0]
        res = res.split("\t")
        name, length = res[1].split(":")[1], int(res[2].split(":")[1])
        logger.debug("nameL %s", name)

    makeRegions(outdir=output_root, gff=gff, name=name,
                dest=rDNA_regions, logger=logger)
    os.makedirs(os.path.join(output_root, "shuffleRegions"))
    flist = getRecLengths(fasta, name=name)
    logger.info("genome coords")
    logger.info(flist)
    with open(os.path.join(output_root, "bedtools_genome"), "w") as bg:
        for entry in flist:
            line = "\t".join([str(x) for x in entry]) + "\n"
            logger.info(line)
            bg.write(line)
    bedtools_cmds, bed_results_list = makeBedtoolsShuffleCmd(
        region=rDNA_regions,
        bedtools_exe=bedtools_exe,
        destdir=os.path.join(output_root, "shuffleRegions"),
        n=args.n_samples, genome=os.path.join(output_root, "bedtools_genome"))

    ref_depth_path, sample_depths_paths, samtools_cmds = samtoolsGetDepths(
        samtools_exe=samtools_exe,
        sample_file_list=bed_results_list,
        bam=args.bam,
        ref_reg_file=rDNA_regions,
        outdir=output_root )

    for cmd_list in [bedtools_cmds, samtools_cmds]:
        for cmd in cmd_list:
            logger.debug(cmd)
            subprocess.run(
                cmd,
                shell=sys.platform != "win32",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True)

    ref_depths = []
    with open(ref_depth_path, "r") as r:
        for line in r:
            ref_depths.append(int(line.strip().split("\t")[2]))
    sample_depths_list = []
    for sample in sample_depths_paths:
        depths = []
        with open(sample, "r") as r:
            for line in r:
                try:
                    depths.append(int(line.strip().split("\t")[2]))
                except:
                    logger.warning("could not parse the following line:")
                    logger.warning(line)
                    depths.append(0)
        sample_depths_list.append(depths)
    logger.debug("first 100 reference depths")
    logger.debug(ref_depths[0:100])
    printPlot(
        data=ref_depths, line=round(stats.mean(ref_depths), 2), ymax=10, xmax=60,
        tick=.2, title="reference", fill=False, logger=logger)
    sample_means = []
    logger.debug("first 100 sample depths")
    for sample in sample_depths_list:
        logger.debug(sample[1:100])
        if len(sample) != 0:  # dont average in zeros
            sample_means.append(stats.mean(sample))
    printPlot(
        data=sample_means, line=round(stats.mean(sample_means), 2), ymax=10,
        xmax=60, tick=.2, title="samples 1-10", fill=False, logger=logger)

    print("Average depth in rDNA regions:\t%0.2f" % stats.mean(ref_depths))
    print(str("Average depth in %d sets of randomly sampled " +
              "non-rDNA regions:\t%0.2f") % (args.n_samples, stats.mean(sample_means)))
    with open(os.path.join(output_root, "riboStack_results.txt"), "w") as outf:
        outf.write("{0}\t{1}\t{2}\t{3}\n".format(
            args.experiment_name, "rDNA_mean_depth", stats.mean(ref_depths), stats.stdev(ref_depths)))
        outf.write("{0}\t{1}\t{2}\t{3}\n".format(
            args.experiment_name, "rDNA_mean_depth", stats.mean(sample_means), stats.stdev(sample_means)))

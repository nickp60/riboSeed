#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
"""

import argparse
import sys
import time
import glob
import re
import math
from bisect import bisect
import shutil
import subprocess
import logging
import os
import traceback

# need this line for unittesting
sys.path.append(os.path.join('..', 'riboSeed'))
from pyutilsnrw.utils3_5 import set_up_logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --------------------------- classes --------------------------- #


# --------------------------- methods --------------------------- #


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="This facilitates the mapping of reads to a reference and comparison of coverage depths in rDNA regions to assess disparity in rDNA counts between the reference and your reads",
        add_help=False)  # to allow for custom help
    parser.add_argument("riboScan_dir", action="store",
                        help="We need the gff and fasta files from your riboScan run.")

    # taking a hint from http://stackoverflow.com/questions/24180527
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-b", "--bam", dest='bam', action="store",
                               help="BAM file; tested with BWA output; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-n", "--name", dest='name',
                          action="store",
                          help="name to give the contig files; "
                          "default: infered from file",
                          type=str)
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
    args = parser.parse_args()
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

def bedtoolsShuffle(region, destdir, genome, bedtools_exe, n=10, logger=None):
    for i in range(0, n):
        bedcmd = "{0} shuffle -i {1} -g {2} > {3}_{4}".format(
            bedtools_exe, region, genome,
            os.path.join(destdir, "sample_region"), i)
        if logger:
            logger.debug("cmd to shuffle regions: %s", bedcmd)
        subprocess.run(bedcmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)

def samtoolsGetDepths(samtools_exe, ref, bam, ref_reg_file, sample_reg_dir, outdir ):
    """ return list [ref_file, [list, of, sample, files]]
    """

    sample_files = glob.glob(os.path.join(sample_reg_dir, "*"))
    ref_out = os.path.join(outdir, "ref_out_depth")
    cmds = []
    ref_cmd = "{0} view -b -F 256 {1} | {0} depth -b {2} - > {3}".format(
        samtools_exe,
        bam,
        ref_reg_file,
        ref_out)
    cmds.append(ref_cmd)
    sample_outs = []
    for idx, f in enumerate(sample_files):
        sample_out = os.path.join(outdir, "sample_out_depth_" + str(idx))
        cmd = "{0} view -b -F 256 {1} | {0} depth -b {2} - > {3}".format(
            samtools_exe,
            bam,
            f,
            sample_out)
        cmds.append(cmd)
        sample_outs.append(sample_out)
    for cmd in cmds:
        if logger:
            logger.debug("cmd to get depths regions: %s", cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    return [ref_out, sample_outs]


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
    if sli == 0 or len(data) == 0 :
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


def mean(x):
    return float(sum(x)) / max(len(x), 1)


if __name__ == "__main__":
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
                            str("riboStack_log_{0}.txt".format(
                                time.strftime("%Y%m%d%H%M"))))
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=log_path,
                            name=__name__)
    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))

    # check exes
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
    logger.info("Shuffling with bedtools")
    os.makedirs(os.path.join(output_root, "shuffleRegions"))
    flist = getRecLengths(fasta, name=name)
    logger.info("genome coords")
    logger.info(flist)
    with open(os.path.join(output_root, "bedtools_genome"), "w") as bg:
        for entry in flist:
            line = "\t".join([str(x) for x in entry]) + "\n"
            logger.info(line)
            bg.write(line)
    bedtoolsShuffle(region=rDNA_regions,
                    bedtools_exe=bedtools_exe,
                    destdir=os.path.join(output_root, "shuffleRegions"),
                    n=10, genome=os.path.join(output_root, "bedtools_genome"),
                    logger=logger)



    ref_depth_path, sample_depths_paths = samtoolsGetDepths(
        samtools_exe=samtools_exe,
        ref=fasta,
        bam=args.bam,
        ref_reg_file=rDNA_regions,
        sample_reg_dir=os.path.join(output_root, "shuffleRegions"),
        outdir=output_root )

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
    logger.debug("regerence depths")
    logger.debug(ref_depths)
    printPlot(data=ref_depths, line=mean(ref_depths), ymax=10, xmax=60, tick=.2,
              title="reference", fill=False, logger=logger)
    sample_means = []
    logger.debug("sample depths")
    for sample in sample_depths_list:
        logger.debug(sample)
        sample_means.append(mean(sample))
    printPlot(data=sample_means, line=round(mean(sample_means), 2), ymax=10, xmax=60, tick=.2,
              title="samples 1-10", fill=False, logger=logger)

    print("Average depth in rDNA regions:\t%d" % mean(ref_depths))
    print("Average depth in 10 sets of randomly sampled non-rDNA regions:\t%d" % mean(sample_means))

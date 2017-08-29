#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
Created on Sun Jul 24 19:33:37 2016

See README.md for more info and usage
"""

import argparse
import sys
import time
import os
import shutil
import pkg_resources

try:  # development mode
    from _version import __version__
except ImportError:  # ie, if an installed pkg from pip or other using setup.py
    __version__ = pkg_resources.require("riboSeed")[0].version

def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="Given cluster file of rDNA regions from riboSelect and " +
        "either paired-end or single-end reads, assembles the mapped reads " +
        "into pseduocontig 'seeds', and uses those with the reads to run" +
        "de fere novo and de novo assembly with SPAdes",
        add_help=True)
    parser.add_argument("-o", "--outdir", dest='outdir', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=False)
    args = parser.parse_args()
    return args


def config_exes():
    # find all needed system requiremnts
    req_programs = [("BARRNAP_EXE", "barrnap"),
                    ("SEQRET_EXE", "seqret"),
                    ("SPADES_EXE", "spades.py"),
                    ("BWA_EXE", "bwa"),
                    ("SAMTOOLS_EXE", "samtools"),
                    ("BALST_EXE", "blastn"),
                    ("MAUVE_ALIGNER", "mauveAligner")]
    config_lines = ["# Required programs"]
    for k, v in req_programs:
        config_lines.append("# executable for " + v)
        if shutil.which(v):
            config_lines.append(k + " = '" + shutil.which(v) + "'\n")
        else:
            config_lines.append(k + " = None\n" )

    # # find all optional sys requirements
    # opt_programs = ["quast.py", "smalt"]
    # # set options based on system requirements
    # pass
    return config_lines

def config_run():
    # set library type
    # se
    pass

def config_bug():
    pass


def make_config_header():
    hashbang = ["#!/usr/bin/env python3",
                "#-*- coding: utf-8 -*-"]
    copy_header = [
        "# Copyright 2017, National University of Ireland and The James Hutton Insitute",
        "# Author: Nicholas Waters",
        "#",
        "# This code is part of the riboSeed package, and is governed by its licence.",
        "# Please see the LICENSE file that should have been included as part of",
        "# this package.\n\n"]
    t0 = time.asctime()
    timestamp = "this config file was generated " + str(t0) + "\n\n"
    hashbang.extend(copy_header)
    hashbang.append(timestamp)
    return(hashbang)


def write_config(header, outfile):
    with open(outfile, "w") as of:
        for hline in header:
            of.write(hline)
            of.write("\n")


def append_config(lines, outfile):
    with open(outfile, "a") as of:
        for line in lines:
            of.write(line)
            of.write("\n")


def main(args):
    lines = config_exes()
    header = make_config_header()
    outfile = os.path.join(os.path.abspath(os.path.expanduser(args.outdir)),
                           str(time.strftime("%Y-%m-%dT%H:%M") +
                               "_riboSeed_config.py"))

    write_config(header=header, outfile=outfile)
    append_config(lines=lines, outfile=outfile)

if __name__ == "__main__":
    args = get_args()
    main(args)

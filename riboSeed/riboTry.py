#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

import pkg_resources
import sys
import os
import shutil
import subprocess
import argparse

from .shared_methods import set_up_logging

helpstring = """
Welcome to the ribo try! Here we test the integration of several parts of the
riboSeed pipeline.  First, `ribo run` is performed on the included test
dataset.  Then, essentially the same thing is done, but calling the
individual subcommands (`ribo scan`, `ribo select`, etc)

If all goes well, no errors should occur, and you should essentially have
two "identical" riboSeed assemblies (although due to random assignments
of mapping duplicates, the nature of error correction, etc, I can't
guarantee that you will get the exact same result

Have fun!
"""


def get_args(test_args=None):  # pragma: no cover
    parser = argparse.ArgumentParser(
        prog="ribo try",
        description=helpstring,
        add_help=False)  # to allow for custom help
    parser.prog = "ribo try"
    parser.add_argument("-o", "--output", dest='output', action="store",
                        help="output directory; " +
                        "default: %(default)s",
                        default=os.path.join(
                            os.getcwd(), "riboSeed_sample_results"),
                        type=str)
    parser.add_argument("-v", "--verbosity", dest='verbosity',
                        action="store",
                        default=2, type=int, choices=[1, 2, 3, 4, 5],
                        help="Logger writes debug to file in output dir; " +
                        "this sets verbosity level sent to stderr. " +
                        " 1 = debug(), 2 = info(), 3 = warning(), " +
                        "4 = error() and 5 = critical(); " +
                        "default: %(default)s")
    parser.add_argument("-c", "--cores", dest='cores', action="store",
                        default=2, type=int,
                        help="cores to be used" +
                        "; default: %(default)s")
    parser.add_argument("-t", "--threads", dest='threads',
                        action="store",
                        default=1, type=int,
                        choices=[1, 2, 4],
                        help="if your cores are hyperthreaded, set number" +
                        " threads to the number of threads per processer." +
                        "If unsure, see 'cat /proc/cpuinfo' under 'cpu " +
                        "cores', or 'lscpu' under 'Thread(s) per core'." +
                        ": %(default)s")
    parser.add_argument("-m", "--memory", dest='memory', action="store",
                        default=8, type=int,
                        help="system memory available" +
                        "; default: %(default)s")
    parser.add_argument("-h", "--help",
                        action="help", default=argparse.SUPPRESS,
                        help="Displays this help message")
    args = parser.parse_args(sys.argv[2:])
    return args


def main(args, logger=None):
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory %s already exists; exiting..." % output_root)
        sys.exit(1)

    log_path = os.path.join(output_root, "riboTry.log")
    if logger is None:
        logger = set_up_logging(verbosity=args.verbosity,
                                outfile=log_path,
                                name=__name__)
    logger.info("Testing your installation of riboSeed on some test data")
    # here we locate the test data we packaged with riboSeed -
    # some reads and a reference
    resource_package = pkg_resources.Requirement.parse("riboSeed")
    logger.debug(resource_package)
    # this looks like I should be using os.path.join, but the package resource
    # stuff needs unix-style path seps

    resource_path_fasta = '/'.join(('riboSeed',
                                    'integration_data', 'concatenated_seq.fasta'))
    resource_path_reffasta = '/'.join(('riboSeed',
                                       'integration_data', 'NC_000913.3.fasta'))
    resource_path_1 = '/'.join(('riboSeed',
                                'integration_data', 'test_reads1.fq'))
    resource_path_2 = '/'.join(('riboSeed',
                                'integration_data', 'test_reads2.fq'))

    logger.debug(resource_path_fasta)
    fasta = pkg_resources.resource_filename(resource_package, resource_path_fasta)
    reffasta = pkg_resources.resource_filename(resource_package,
                                               resource_path_reffasta)
    fastq1 = pkg_resources.resource_filename(resource_package, resource_path_1)
    fastq2 = pkg_resources.resource_filename(resource_package, resource_path_2)
    # fasta_path = pkg_resources.resource_string("/", resource_path)
    logger.debug(fasta)
    logger.debug(reffasta)
    logger.debug(fastq1)
    logger.debug(fastq2)

    for i in ["blastn", "spades.py", "bwa", "mafft",
              "samtools", "barrnap"]:
        assert shutil.which(i) is not None, \
            "{0} executable not found in PATH!".format(i)

    ribo_run_cmd = str(
        "ribo run -r {0} -o {1} -F {2} -R {3} --serialize -v 1 " +
        "--subassembler skesa  " +
        "--stages stack score spec --cores {4} --threads {5} --memory {6}"
    ).format(
        fasta,
        os.path.join(output_root, "run"),
        fastq1,
        fastq2,
        args.cores,
        args.threads,
        args.memory)
    logger.info("running " + ribo_run_cmd)
    logger.info("This usually take about ~4-5 minutes to run all the modules")
    subprocess.run([ribo_run_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    logger.info("finished running integration test with ribo run!")

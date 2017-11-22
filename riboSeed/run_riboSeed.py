#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

import argparse
import yaml
import pkg_resources
import time
import sys
import os
import shutil
import subprocess
import multiprocessing
import re

from .shared_methods import set_up_logging
from argparse import Namespace

from ._version import __version__
from . import riboScan as rscan
from . import riboSelect as rsel
from . import riboSeed as rseed
from . import riboScore as rscore
# note that we import sketch later in case of matplotlib gotchas
from . import make_riboSeed_config as mrc


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        prog="ribo run",
        description="Run the riboSeed pipeline of scan, select, seed, " +
        "sketch, and score.  Uses a config file to wrangle all the args not " +
        "available via these commandline args.\n\n\n" +
        "This can either be run by providing (as minimum) a reference, " +
        "some reads, and an output directory; or, if you have a completed " +
        "config file, you can run it with just that.",
        formatter_class=argparse.MetavarTypeHelpFormatter,
        add_help=False)  # to allow for custom help
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-r", "--reference_fasta",
                          dest="REFERENCE_FASTA", action="store",
                          help="path to a (multi)fasta or a directory " +
                          "containing one or more chromosomal " +
                          "sequences in fasta format. Required, unless " +
                          "using a config file", type=str,
                          metavar="reference.fasta",
                          default=None)
    optional.add_argument("-c", "--config", dest='RUN_CONFIG', action="store",
                          help="config file; if none given, create one; " +
                          "default: %(default)s", default=os.getcwd(),
                          type=str, required="-r" not in sys.argv,
                          metavar="config_file")
    optional.add_argument("-o", "--output", dest='RUN_OUTPUT', action="store",
                          help="output directory; " +
                          "default: %(default)s",
                          default=os.path.join(
                              os.getcwd(),
                              str(time.strftime("%Y-%m-%dT%H%M") +
                                  "_riboSeed_pipeline_results"), ""),
                          type=str,
                          metavar="/output/dir/")
    optional.add_argument("-n", "--experiment_name",
                          dest='RUN_EXPERIMENT_NAME',
                          action="store",
                          help="prefix for results files; " +
                          "default: inferred",
                          default=None, type=str,
                          metavar="experiment_name")
    # riboScan args
    optional.add_argument("-K", "--Kingdom", dest='RUN_KINGDOM',
                          action="store",
                          choices=["bac", "euk", "arc", "mito"],
                          help="whether to look for eukaryotic, archaeal, or" +
                          " bacterial rDNA; " +
                          "default: %(default)s", default="bac",
                          type=str)
    optional.add_argument("-S", "--specific_features",
                          dest="RUN_SPECIFIC_FEATURES",
                          help="colon:separated -- specific features"
                          "; default: %(default)s",
                          default='16S:23S:5S', type=str,
                          metavar='16S:23S:5S')
    #riboSelect args
    optional.add_argument("--clusters",
                          help="number of rDNA clusters;"
                          "if submitting multiple records, must be a "
                          "colon:separated list whose length matches number "
                          "of genbank records.  Default is inferred from "
                          "specific feature with fewest hits", default='',
                          type=str, dest="RUN_CLUSTERS")
    # riboSeed args
    optional.add_argument("-C", "--cluster_file",
                          dest="RUN_CLUSTER_FILE",
                          help="clustered_loci file output from riboSelect;"
                          "this is created by default from run_riboSeed, "
                          "but if you don't agree with the operon structure "
                          "predicted by riboSelect, you can use your " +
                          "alternate clustered_loci file. " +
                          "default: %(default)s",
                          default=None,
                          type=str)
    optional.add_argument("-F", "--fastq1", dest='RUN_FASTQ1', action="store",
                          help="path to forward fastq file, can be compressed",
                          type=str, default=None,
                          metavar="reads_F.fq")
    optional.add_argument("-R", "--fastq2", dest='RUN_FASTQ2', action="store",
                          help="path to reverse fastq file, can be compressed",
                          type=str, default=None,
                          metavar="reads_R.fq")
    optional.add_argument("-S1", "--fastq_single1", dest='RUN_FASTQS1',
                          action="store",
                          help="path to single fastq file", type=str,
                          default=None,
                          metavar="reads_S.fq")
    optional.add_argument("-s", "--score_min", dest='RUN_SCORE_MIN',
                          action="store",
                          default=None, type=int,
                          help="If using smalt, this sets the '-m' param; " +
                          "default with smalt is inferred from " +
                          "read length. If using BWA, reads mapping with AS" +
                          "score lower than this will be rejected" +
                          "; default with BWA is half of read length")
    optional.add_argument("--ref_as_contig", dest='RUN_REF_AS_CONTIG',
                          action="store", type=str,
                          default="infer",
                          choices=["ignore", "infer", "trusted", "untrusted"],
                          help="ignore: reference will not be used in " +
                          "subassembly. trusted: SPAdes will use the seed" +
                          " sequences as a --trusted-contig; untrusted: " +
                          "SPAdes will treat as --untrusted-contig. " +
                          "infer: if mapping percentage " +
                          "over 80%%, 'trusted'; else 'untrusted'." +
                          " See SPAdes docs for details.  default: infer")
    optional.add_argument("--linear",
                          help="if genome is known to not be circular and " +
                          "a region of interest (including flanking bits) " +
                          "extends past chromosome end, this extends the " +
                          "seqence past chromosome origin forward by " +
                          "--padding; " +
                          "default: %(default)s",
                          default=False, dest="RUN_LINEAR",
                          action="store_true")
    optional.add_argument("-j", "--just_seed", dest='RUN_JUST_SEED',
                          action="store_true",
                          default=False,
                          help="Don't do an assembly, just generate the long" +
                          " read 'seeds'; default: %(default)s")
    optional.add_argument("--score", dest='RUN_SCORE',
                          action="store_true",
                          default=False,
                          help="run riboScore too! " +
                          "default: %(default)s")
    # optional.add_argument("--sketch", dest='RUN_SKETCH',
    #                       action="store_true",
    #                       default=False,
    #                       help="run riboSketch too! " +
    #                       "default: %(default)s")
    optional.add_argument("-l", "--flanking_length",
                          help="length of flanking regions, in bp; " +
                          "default: %(default)s",
                          default=1000, type=int, dest="RUN_FLANKING")
    optional.add_argument("-k", "--kmers", dest='RUN_KMERS', action="store",
                          default="21,33,55,77,99,127", type=str,
                          help="kmers used for final assembly" +
                          ", separated by commas such as" +
                          "21,33,55,77,99,127. Can be set to 'auto', where " +
                          "SPAdes chooses.  We ensure kmers are not " +
                          "too big or too close to read length" +
                          "; default: %(default)s",
                          metavar="21,33,55,77,99,127")
    optional.add_argument("--force_kmers", dest="RUN_FORCE_KMERS",
                          action="store_true",
                          default=False,
                          help="skip checking to see if kmerchoice is " +
                          "appropriate to read length. Sometimes kmers " +
                          "longer than reads can help in the final assembly," +
                          " as the long reads generated by riboSeed contain " +
                          "kmers longer than the read length")
    optional.add_argument("-p", "--pre_kmers", dest='RUN_PRE_KMERS',
                          action="store",
                          default="21,33,55,77,99", type=str,
                          help="kmers used during seeding assemblies, " +
                          "separated bt commas" +
                          "; default: %(default)s", metavar="21,33,55,77,99")
    optional.add_argument("-d", "--min_flank_depth",
                          help="a subassembly won't be performed if this " +
                          "minimum depth is not achieved on both the 3' and" +
                          "5' end of the pseudocontig. " +
                          "default: %(default)s",
                          default=0, dest="RUN_MIN_FLANKING_DEPTH", type=int)
    optional.add_argument("--clean_temps", dest='RUN_CLEAN_TEMPS',
                          default=False, action="store_true",
                          help="if --clean_temps, mapping files will be " +
                          "removed once they are no no longer needed during" +
                          " the mapping iterations to save space; " +
                          "default: %(default)s")
    optional.add_argument("-i", "--iterations", dest='RUN_ITERATIONS',
                          action="store",
                          default=3, type=int,
                          help="if iterations>1, multiple seedings will " +
                          "occur after subassembly of seed regions; " +
                          "if setting --target_len, seedings will continue " +
                          "until --iterations are completed or --target_len"
                          " is matched or exceeded; " +
                          "default: %(default)s")
    optional.add_argument("-v", "--verbosity", dest='RUN_VERBOSITY',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("--cores", dest='RUN_CORES', action="store",
                          default=None, type=int,
                          help="cores used" +
                          "; default: %(default)s")
    optional.add_argument("--memory", dest='RUN_MEMORY', action="store",
                          default=8, type=int,
                          help="cores for multiprocessing" +
                          "; default: %(default)s")
    optional.add_argument("--damn_the_torpedos", dest='damn_the_torpedos',
                          action="store_true",
                          default=False,
                          help="Ignore certain  errors, full speed ahead!")
    optional.add_argument("-t", "--threads", dest='RUN_THREADS',
                          action="store",
                          default=1, type=int,
                          choices=[1, 2, 4],
                          help="if your cores are hyperthreaded, set number" +
                          " threads to the number of threads per processer." +
                          "If unsure, see 'cat /proc/cpuinfo' under 'cpu " +
                          "cores', or 'lscpu' under 'Thread(s) per core'." +
                          ": %(default)s")
    optional.add_argument("-z", "--serialize", dest='serialize',
                          action="store_true",
                          default=False,
                          help="if --serialize, runs seeding and assembly " +
                          "without multiprocessing. We recommend this for " +
                          "machines with less than 8GB RAM: %(default)s")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    optional.add_argument('--version', action='version',
                          version='%(prog)s {version}'.format(
                              version=__version__))
    args = parser.parse_args(sys.argv[2:])
    return args


def detect_or_create_config(config_file, output_root, theseargs,
                            newname=None, logger=None):
    """ return path to config file, or die trying.
    if a config file exists, just return the path.  If not,
    create it using "make_riboSeed_config.py"
    """
    assert logger is not None, "must use logging"
    if not os.path.isfile(config_file):
        logger.debug("creating config file")
        make_config_args = Namespace(
            outdir=output_root,
            name=newname)
        config_file = mrc.main(make_config_args)
        add_these_params_to_config(config_file=config_file,
                                   args=theseargs)
    else:
        logger.info("using provided config file! ignoring any other args" +
                    "provided via commandline")
    return config_file


def add_these_params_to_config(config_file, args):
    """ write the args used by run_riboSeed to config
    """
    with open(config_file, 'a') as stream:
        yaml.dump(vars(args), stream, default_flow_style=False)


def parse_config(config_file, logger=None):
    """ Read the config file, make a namespace object out of it
    """
    with open(config_file, 'r') as stream:
        try:
            yamldata = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            logger.error(exc)
            logger.error("error parsing config file!")
            sys.exit(1)
    logger.debug(yamldata)
    newns = argparse.Namespace()
    newns.__dict__ = yamldata
    return(newns)


def new_log_for_diff(logfile_path):
    """ make new log file without DDDD-MM-DD HH:MM:SS in line
    sometimes you want timestamps, sometimes you just wanna diff stuff
    """
    pattern  = re.compile(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}")
    with open(os.path.join(os.path.dirname(logfile_path),
                           "run_riboSeed_notime.log"), "w") as outlog:
        with open(logfile_path, "r") as inlog:
            for line in inlog:
                splitline = re.sub(pattern, "", line)
                outlog.write(splitline)


def main(args):
    # set up output directory and logging
    output_root = os.path.abspath(os.path.expanduser(args.RUN_OUTPUT))
    try:
        os.makedirs(output_root)
    except FileExistsError:
        print("Selected output directory %s exists" %
              output_root)
        sys.exit(1)
    log_path = os.path.join(output_root, "run_riboSeed.log")
    logger = set_up_logging(verbosity=args.RUN_VERBOSITY,
                            outfile=log_path,
                            name=__name__)
    logger.info("Usage:\n%s\n", " ".join([x for x in sys.argv]))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, str(v))
    # if not using a config, groom some of the arguments
    if args.REFERENCE_FASTA is not None:
        # detect system attributes if not given expressly
        if args.RUN_CORES is None:
            args.RUN_CORES = multiprocessing.cpu_count()
            logger.info("Using %i cores", args.RUN_CORES)

    # if starting a fresh run, create a config file.  If not, read it in
    args.RUN_CONFIG = detect_or_create_config(
        config_file=args.RUN_CONFIG,
        theseargs=args,
        output_root=output_root, logger=logger)
    conf = parse_config(args.RUN_CONFIG, logger=logger)
    # if no name given, infer from the name of the contigs file, from either
    # args or the config file
    if conf.RUN_EXPERIMENT_NAME is not None:
        experiment_name = conf.RUN_EXPERIMENT_NAME
    elif conf.REFERENCE_FASTA is not  None:
        experiment_name = os.path.basename(
            os.path.splitext(conf.REFERENCE_FASTA)[0])
    else:  # infer from riboScan arg
        assert conf.SCAN_CONTIG_NAME is not None, \
            "no reference fasta found in args or config! Exiting"
        experiment_name = os.path.basename(
            os.path.splitext(conf.SCAN_CONTIG_NAME)[0])

    logger.debug(experiment_name)

    if conf.RUN_CLUSTER_FILE is None:
        cluster_txt_file = os.path.join(
            output_root, "select", "riboSelect_grouped_loci.txt")
    else:
        cluster_txt_file = conf.RUN_CLUSTER_FILE

    scan_args = Namespace(
        contigs=conf.REFERENCE_FASTA,
        output=os.path.join(output_root, "scan"),
        kingdom=conf.RUN_KINGDOM,
        id_thresh=conf.SCAN_ID_THRESH,
        name=conf.SCAN_CONTIG_NAME,
        cores=conf.RUN_CORES,
        barrnap_exe=conf.BARRNAP_EXE,
        seqret_exe=conf.SEQRET_EXE,
        min_length=conf.SCAN_MIN_LENGTH,
        verbosity=conf.SCAN_VERBOSITY)
    select_args = Namespace(
        genbank_genome=os.path.join(
            output_root, "scan", "scannedScaffolds.gb"),
        output=os.path.join(output_root, "select"),
        feature=conf.SELECT_FEATURE,
        specific_features=conf.RUN_SPECIFIC_FEATURES,
        clobber=False,
        clusters=conf.RUN_CLUSTERS,
        verbosity=conf.SELECT_VERBOSITY,
        debug=False)
    seed_args = Namespace(
        clustered_loci_txt=cluster_txt_file,
        reference_genbank=os.path.join(
            output_root, "scan", "scannedScaffolds.gb"),
        output=os.path.join(output_root, "seed"),
        fastq1=conf.RUN_FASTQ1,
        fastq2=conf.RUN_FASTQ2,
        fastqS1=conf.RUN_FASTQS1,
        just_seed=conf.RUN_JUST_SEED,
        min_flank_depth=conf.RUN_MIN_FLANKING_DEPTH,
        exp_name=experiment_name,
        clean_temps=conf.RUN_CLEAN_TEMPS,
        flanking=conf.RUN_FLANKING,
        method=conf.SEED_MAP_METHOD,
        iterations=conf.RUN_ITERATIONS,
        cores=conf.RUN_CORES,
        threads=conf.RUN_THREADS,
        memory=conf.RUN_MEMORY,
        kmers=conf.RUN_KMERS,
        pre_kmers=conf.RUN_PRE_KMERS,
        force_kmers=conf.RUN_FORCE_KMERS,
        score_min=conf.RUN_SCORE_MIN,
        min_assembly_len=conf.SEED_MIN_ASSEMBLY_LENGTH,
        include_short_contigs=conf.SEED_INCLUDE_SHORTS,
        subtract=conf.SEED_SUBTRACT,
        linear=conf.RUN_LINEAR,
        skip_control=conf.SEED_SKIP_CONTROL,
        target_len=conf.SEED_TARGET_LEN,
        smalt_scoring=conf.SEED_SMALT_SCORING,
        serialize=conf.serialize,
        ref_as_contig=conf.RUN_REF_AS_CONTIG,
        mapper_args=conf.SEED_MAPPER_ARGS,
        initial_consensus=conf.SEED_INITIAL_CONSENSUS,
        spades_exe=conf.SPADES_EXE,
        samtools_exe=conf.SAMTOOLS_EXE,
        smalt_exe=conf.SMALT_EXE,
        bwa_exe=conf.BWA_EXE,
        bcftools_exe=conf.BCFTOOLS_EXE,
        quast_exe=conf.QUAST_EXE,
        verbosity=conf.SEED_VERBOSITY)
    sketch_args = Namespace(
        indir=os.path.join(output_root, "seed","mauve"),
        outdir=os.path.join(output_root, "sketch"),
        assembly_ext=conf.SKETCH_ASSEMBLY_EXT,
        ref_ext=conf.SKETCH_REF_EXT,
        names="{0},{1},{2}".format(
            os.path.basename(conf.REFERENCE_FASTA),
            os.path.basename(experiment_name) + " de fere novo",
            os.path.basename(experiment_name) + " de novo"),
        replot=False,
        mauve_jar=conf.MAUVE_JAR,
        verbosity=conf.SKETCH_VERBOSITY)
    score_args = Namespace(
        indir=os.path.join(output_root, "seed","mauve"),
        output=os.path.join(output_root, "score"),
        flanking=conf.RUN_FLANKING,
        ref_ext=conf.SKETCH_REF_EXT,
        assembly_ext=conf.SKETCH_ASSEMBLY_EXT,
        blast_full=False,
        verbosity=conf.SCORE_VERBOSITY)


    logger.info("\nrunning riboScan\n")
    rscan.main(scan_args, logger)
    if conf.RUN_CLUSTER_FILE is not None:
        logger.info(
            "Skipping riboSelect, using user-supplied cluster file: %s",
            conf.RUN_CLUSTER_FILE)
    else:
        logger.info("\nrunning riboSelect\n")
        rsel.main(select_args, logger)
    logger.info("\nrunning riboSeed\n")
    rseed.main(seed_args, logger)
    # if conf.RUN_SKETCH:
    #     # import this here in case there are issues with mpl''s X windows
    #     from . import riboSketch as rsketch
    #     if conf.MAUVE_JAR is not None:
    #         logger.info("\nrunning riboSketch\n")
    #         rsketch.main(sketch_args, logger=logger)
    #     else:
    #         logger.info(
    #             "Skipping riboSketch: no Mauve.jar found. To fix, " +
    #             "add the path to Mauve.jar in the config file from this " +
    #             "run, and re-run with -c path/to/config.py")
    if conf.RUN_SCORE:
        if conf.BLAST_EXE is not None:
            logger.info("\nrunning riboScore\n")
            rscore.main(score_args, logger=logger)
        else:
            logger.info("Skipping riboScore, as no blastn executable was " +
                        "found in path.")
    new_log_for_diff(logfile_path=log_path)

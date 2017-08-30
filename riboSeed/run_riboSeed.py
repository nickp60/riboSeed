#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import argparse
import pkg_resources
import time
import sys
import os
import shutil
import subprocess
import multiprocessing

from pyutilsnrw.utils3_5 import set_up_logging
from argparse import Namespace
# from make_riboSeed_config import make_config_header, write_config, \
#     append_config

try:  # development mode
    from _version import __version__
except ImportError:  # ie, if an installed pkg from pip or other using setup.py
    __version__ = pkg_resources.require("riboSeed")[0].version

sys.path.append(os.path.join('..', 'riboSeed'))

import riboScan as rscan
import riboSelect as rsel
import riboSeed as rseed


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="Given a directory of one or more chromosomes as fasta " +
        "files, this facilitates reannotation of rDNA regions with Barrnap " +
        " and outputs all sequences as a single, annotated genbank file",
        add_help=False)  # to allow for custom help
    parser.add_argument("contigs", action="store",
                        help="either a (multi)fasta or a directory " +
                        "containing one or more chromosomal " +
                        "sequences in fasta format")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-c", "--config", dest='config', action="store",
                          help="config file; if none given, create one; " +
                          "default: %(default)s", default=os.getcwd(),
                          type=str, required=False)
    optional.add_argument("-o", "--output", dest='output', action="store",
                          help="output directory; " +
                          "default: %(default)s",
                          default=os.path.join(
                              os.getcwd(),
                              str(time.strftime("%Y-%m-%dT%H:%M") +
                                  "_riboSeed_pipeline_results"), ""),
                          type=str)
    optional.add_argument("-n", "--experiment_name", dest='exp_name',
                          action="store",
                          help="prefix for results files; " +
                          "default: %(default)s",
                          default="riboSeed", type=str)
    # riboScan args
    optional.add_argument("-K", "--Kingdom", dest='kingdom',
                          action="store",
                          choices=["bac", "euk", "arc", "mito"],
                          help="whether to look for eukaryotic, archaeal, or" +
                          " bacterial rDNA; " +
                          "default: %(default)s", default="bac",
                          type=str)
    optional.add_argument("-s", "--specific_features",
                          help="colon:separated -- specific features"
                          "; default: %(default)s",
                          default='16S:23S:5S', type=str)
    #riboSelect args
    optional.add_argument("--clusters",
                          help="number of rDNA clusters;"
                          "if submitting multiple records, must be a "
                          "colon:separated list whose length matches number "
                          "of genbank records.  Default is inferred from "
                          "specific feature with fewest hits", default='',
                          type=str, dest="clusters")
    # riboSeed args
    optional.add_argument("-F", "--fastq1", dest='fastq1', action="store",
                          help="forward fastq reads, can be compressed",
                          type=str, default=None)
    optional.add_argument("-R", "--fastq2", dest='fastq2', action="store",
                          help="reverse fastq reads, can be compressed",
                          type=str, default=None)
    optional.add_argument("-S1", "--fastq_single1", dest='fastqS1',
                          action="store",
                          help="single fastq reads", type=str, default=None)
    optional.add_argument("--linear",
                          help="if genome is known to not be circular and " +
                          "a region of interest (including flanking bits) " +
                          "extends past chromosome end, this extends the " +
                          "seqence past chromosome origin forward by " +
                          "--padding; " +
                          "default: %(default)s",
                          default=False, dest="linear", action="store_true")
    optional.add_argument("-j", "--just_seed", dest='just_seed',
                          action="store_true",
                          default=False,
                          help="Don't do an assembly, just generate the long" +
                          " read 'seeds'; default: %(default)s")
    optional.add_argument("-l", "--flanking_length",
                          help="length of flanking regions, in bp; " +
                          "default: %(default)s",
                          default=1000, type=int, dest="flanking")
    optional.add_argument("-k", "--kmers", dest='kmers', action="store",
                          default="21,33,55,77,99,127", type=str,
                          help="kmers used for final assembly" +
                          ", separated by commas such as" +
                          "21,33,55,77,99,127. Can be set to 'auto', where " +
                          "SPAdes chooses.  We ensure kmers are not " +
                          "too big or too close to read length" +
                          "; default: %(default)s")
    optional.add_argument("-p", "--pre_kmers", dest='pre_kmers',
                          action="store",
                          default="21,33,55,77,99", type=str,
                          help="kmers used during seeding assemblies, " +
                          "separated bt commas" +
                          "; default: %(default)s")
    optional.add_argument("-d", "--min_flank_depth",
                          help="a subassembly won't be performed if this " +
                          "minimum depth is not achieved on both the 3' and" +
                          "5' end of the pseudocontig. " +
                          "default: %(default)s",
                          default=0, dest="min_flank_depth", type=float)
    optional.add_argument("--clean_temps", dest='clean_temps',
                          default=False, action="store_true",
                          help="if --clean_temps, mapping files will be " +
                          "removed once they are no no longer needed during" +
                          " the mapping iterations to save space; " +
                          "default: %(default)s")
    optional.add_argument("-i", "--iterations", dest='iterations',
                          action="store",
                          default=3, type=int,
                          help="if iterations>1, multiple seedings will " +
                          "occur after subassembly of seed regions; " +
                          "if setting --target_len, seedings will continue " +
                          "until --iterations are completed or --target_len"
                          " is matched or exceeded; " +
                          "default: %(default)s")
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("--cores", dest='cores', action="store",
                          default=None, type=int,
                          help="cores used" +
                          "; default: %(default)s")
    optional.add_argument("--memory", dest='memory', action="store",
                          default=8, type=int,
                          help="cores for multiprocessing" +
                          "; default: %(default)s")
    optional.add_argument("-t", "--threads", dest='threads',
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
    args = parser.parse_args()
    return args


def run_make_config(output_root, logger=None):
    assert logger is not None, "must use logging"
    cmd = "{0} {1} -o {2}".format(
        sys.executable,
        shutil.which("make_riboSeed_config.py"),
        output_root)
    logger.debug(cmd)
    result = subprocess.run(
        cmd,
        shell=sys.platform != "win32",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True)
    logger.debug(result)
    conf = result.stdout.decode("utf-8").split("\n")[0].split("\t")
    assert len(conf) == 1, "make_riboSeed_config.py writes too much to stdout!"
    return(conf[0])


def parse_config(config_file, output_root, logger=None):
    assert logger is not None, "must use logging"
    if not os.path.isfile(config_file):
        config_file = run_make_config(output_root, logger)
    else:
        pass
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        os.path.splitext(os.path.basename(config_file))[0],
        config_file)
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return(config)


def main(args):
    output_root = os.path.abspath(os.path.expanduser(args.output))
    # Create output directory only if it does not exist
    try:
        os.makedirs(output_root)
    except FileExistsError:
         # leading comment char'#' added for stream output usage
        print("#Selected output directory %s exists" %
              output_root)
        sys.exit(1)
    log_path = os.path.join(output_root, "run_riboSeed.log")
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=log_path,
                            name=__name__)
    if args.cores is None:
        args.cores = multiprocessing.cpu_count()
        logger.info("Using %i cores", args.cores)
    conf = parse_config(args.config, output_root, logger=logger)
    scan_args = Namespace(
        contigs=args.contigs,
        output=os.path.join(output_root, "scan"),
        kingdom=args.kingdom,
        id_thresh=conf.SCAN_ID_THRESH,
        name=conf.SCAN_CONTIG_NAME,
        cores=args.cores,
        barrnap_exe=conf.BARRNAP_EXE,
        seqret_exe=conf.SEQRET_EXE,
        min_length=conf.SCAN_MIN_LENGTH,
        verbosity=conf.SCAN_VERBOSITY)
    select_args = Namespace(
        genbank_genome=os.path.join(
            output_root, "scan", "scannedScaffolds.gb"),
        output=os.path.join(output_root, "select"),
        feature=conf.SELECT_FEATURE,
        specific_features=conf.SELECT_SPECIFIC_FEATURES,
        clobber=False,
        clusters=args.clusters,
        verbosity=conf.SELECT_VERBOSITY,
        debug=False)
    seed_args = Namespace(
        clustered_loci_txt=os.path.join(
            output_root, "select", "riboSelect_grouped_loci.txt"),
        reference_genbank=os.path.join(
            output_root, "scan", "scannedScaffolds.gb"),
        output=os.path.join(output_root, "seed"),
        fastq1=args.fastq1,
        fastq2=args.fastq2,
        fastqS1=args.fastqS1,
        just_seed=args.just_seed,
        min_flank_depth=args.min_flank_depth,
        exp_name=args.exp_name,
        clean_temps=args.clean_temps,
        flanking=args.flanking,
        method=conf.SEED_MAP_METHOD,
        iterations=args.iterations,
        cores=args.cores,
        threads=args.threads,
        memory=args.memory,
        kmers=args.kmers,
        pre_kmers=args.pre_kmers,
        score_min=conf.SEED_SCORE_MIN,
        min_assembly_len=conf.SEED_MIN_ASSEMBLY_LENGTH,
        include_short_contigs=conf.SEED_INCLUDE_SHORTS,
        subtract=conf.SEED_SUBTRACT,
        linear=args.linear,
        skip_control=conf.SEED_SKIP_CONTROL,
        target_len=conf.SEED_TARGET_LEN,
        smalt_scoring=conf.SEED_SMALT_SCORING,
        serialize=args.serialize,
        ref_as_contig=conf.SEED_REF_AS_CONTIG,
        mapper_args=conf.SEED_MAPPER_ARGS,
        spades_exe=conf.SPADES_EXE,
        samtools_exe=conf.SAMTOOLS_EXE,
        smalt_exe=conf.SMALT_EXE,
        bwa_exe=conf.BWA_EXE,
        quast_exe=conf.QUAST_EXE,
        verbosity=conf.SEED_VERBOSITY)
    logger.info("\nrunning riboScan\n")
    rscan.main(scan_args, logger)
    logger.info("\nrunning riboSelect\n")
    rsel.main(select_args, logger)
    logger.info("\nrunning riboSeed\n")
    rseed.main(seed_args, logger)

if __name__ == "__main__":
    args = get_args()
    main(args)

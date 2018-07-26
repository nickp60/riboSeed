#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""
Created on Sun Jul 24 19:33:37 2016

See README.md for more info and usage
"""

import argparse
import sys
import time
import random
import os
import re
import shutil
import multiprocessing
import subprocess
import traceback
import pysam
import math
import pkg_resources

from riboSeed import __version__
from .classes import SeedGenome, LociMapping, Exes, NgsLib
from bisect import bisect
from itertools import chain
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from distutils.version import StrictVersion

# plotting with mpl is depreciated till I can figure out why
# it wont work in screen sessions ( see bioconda issue  #6451)

# try:
#     import numpy as np
#     import matplotlib as mpl
#     from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
#     from matplotlib.figure import Figure
#     import matplotlib.patches as patches
#     PLOT = True
# except Exception as e:  # most likely an ImportError, but Im not taking chances
#     print(e)
#     print("\nlooks like you have some issue with matplotlib.  " +
#           "Classic matplotlib, amirite? Plotting is disabled\n")
#     PLOT = False


from .shared_methods import parse_clustered_loci_file, pad_genbank_sequence, \
    combine_contigs, file_len, get_number_mapped, \
    keep_only_first_contig, get_fasta_lengths, \
    extract_coords_from_locus, add_gb_seqrecords_to_cluster_list, \
    set_up_logging, check_version_from_cmd

# GLOBALS
SAMTOOLS_MIN_VERSION = '1.3.1'


def get_args(test_args=None):  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(prog="ribo seed",
        description="Given cluster file of rDNA regions from riboSelect and " +
        "either paired-end or single-end reads, assembles the mapped reads " +
        "into pseduocontig 'seeds', and uses those with the reads to run" +
        "de fere novo and de novo assembly with SPAdes",
        add_help=False)  # to allow for custom help
    parser.add_argument("clustered_loci_txt", action="store",
                        help="output from riboSelect")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-r", "--reference_genbank",
                               dest='reference_genbank',
                               action="store", default='', type=str,
                               help="genbank reference, used to estimate " +
                               "insert sizes, and compare with QUAST",
                               required=True)
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    # had to make this faux "optional" parse so that the named required ones
    # above get listed aboce the ther args when displaying the help message
    # read libraries
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-F", "--fastq1", dest='fastq1', action="store",
                          help="forward fastq reads, can be compressed",
                          type=str, default=None)
    optional.add_argument("-R", "--fastq2", dest='fastq2', action="store",
                          help="reverse fastq reads, can be compressed",
                          type=str, default=None)
    optional.add_argument("-S1", "--fastq_single1", dest='fastqS1',
                          action="store",
                          help="single fastq reads", type=str, default=None)
    # parameters for run
    optional.add_argument("-l", "--flanking_length",
                          help="length of flanking regions, in bp; " +
                          "default: %(default)s",
                          default=1000, type=int, dest="flanking")
    optional.add_argument("-j", "--just_seed", dest='just_seed',
                          action="store_true",
                          default=False,
                          help="Don't do an assembly, just generate the long" +
                          " read 'seeds'; default: %(default)s")
    optional.add_argument("-e", "--experiment_name", dest='experiment_name',
                          action="store",
                          help="prefix for results files; " +
                          "default: %(default)s",
                          default="riboSeed", type=str)
    optional.add_argument("--mapper", dest='mapper',
                          action="store", choices=["smalt", "bwa"],
                          help="available mappers: smalt and bwa; " +
                          "default: %(default)s",
                          default='bwa', type=str)
    optional.add_argument("-k", "--kmers", dest='kmers', action="store",
                          default="21,33,55,77,99,127", type=str,
                          help="kmers used for final assembly" +
                          ", separated by commas such as" +
                          "21,33,55,77,99,127 . Can be set to 'auto', where " +
                          "SPAdes chooses.  We ensure kmers are not " +
                          "too big or too close to read length" +
                          "; default: %(default)s")
    optional.add_argument("-p", "--pre_kmers", dest='pre_kmers',
                          action="store",
                          default="21,33,55,77,99", type=str,
                          help="kmers used during seeding assemblies, " +
                          "separated bt commas" +
                          "; default: %(default)s")
    optional.add_argument("--force_kmers", dest="force_kmers",
                          action="store_true",
                          default=False,
                          help="skip checking to see if kmerchoice is " +
                          "appropriate to read length. Sometimes kmers " +
                          "longer than reads can help in the final assembly," +
                          " as the long reads generated by riboSeed contain " +
                          "kmers longer than the read length")
    optional.add_argument("-s", "--score_min", dest='score_min',
                          action="store",
                          default=None, type=int,
                          help="If using smalt, this sets the '-m' param; " +
                          "default with smalt is inferred from " +
                          "read length. If using BWA, reads mapping with AS" +
                          "score lower than this will be rejected" +
                          "; default with BWA is half of read length")
    optional.add_argument("-a", "--min_assembly_len", dest='min_assembly_len',
                          action="store",
                          default=6000, type=int,
                          help="if initial SPAdes assembly largest contig " +
                          "is not at least as long as --min_assembly_len, " +
                          "reject. Set this to the length of the seed " +
                          "sequence; if it is not achieved, seeding across " +
                          "regions will likely fail; default: %(default)s")
    optional.add_argument("--include_shorts", dest='include_short_contigs',
                          action="store_true",
                          default=False,
                          help="if assembled contig is smaller than  " +
                          "--min_assembly_len, contig will still be included" +
                          " in assembly; default: inferred")
    optional.add_argument("--damn_the_torpedos", dest='damn_the_torpedos',
                          action="store_true",
                          default=False,
                          help="Ignore certain errors, full speed ahead!")
    optional.add_argument("--subtract", dest='subtract',
                          action="store_true",
                          default=False,
                          help="if --subtract reads already used in previous" +
                          "round of subassembly will not be included in " +
                          "subsequent rounds.  This can lead to problems " +
                          "with multiple mapping and inflated coverage.")
    optional.add_argument("--linear",
                          help="if genome is known to not be circular and " +
                          "a region of interest (including flanking bits) " +
                          "extends past chromosome end, this extends the " +
                          "seqence past chromosome origin forward by " +
                          "--padding; " +
                          "default: %(default)s",
                          default=False, dest="linear", action="store_true")
    optional.add_argument("-d", "--min_flank_depth",
                          help="a subassembly will not be performed if this " +
                          "minimum depth is not achieved on both the 3' and" +
                          "5' end of the pseudocontig. " +
                          "default: %(default)s",
                          default=0, dest="min_flank_depth", type=float)
    optional.add_argument("--ref_as_contig", dest='ref_as_contig',
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
    optional.add_argument("--additional_libs",
                          help="include these libraries in  final assembly " +
                          "in addition to the reads supplied as -F and -R. " +
                          "They must be supplied according to SPAdes arg " +
                          "naming scheme. Use at own risk."
                          "default: %(default)s",
                          dest="additional_libs", type=str)
    optional.add_argument("--clean_temps", dest='clean_temps',
                          default=False, action="store_true",
                          help="if --clean_temps, mapping files will be " +
                          "removed once they are no no longer needed during " +
                          "the mapping iterations to save space; " +
                          "default: %(default)s")
    optional.add_argument("--skip_control", dest='skip_control',
                          action="store_true",
                          default=False,
                          help="if --skip_control, no de novo " +
                          "assembly will be done; default: %(default)s")
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
    optional.add_argument("--target_len", dest='target_len', action="store",
                          default=None, type=float,
                          help="if set, iterations will continue until " +
                          "contigs reach this length, or max iterations (" +
                          "set by --iterations) have been completed. Set as " +
                          "fraction of original seed length by giving a " +
                          "decimal between 0 and 5, or set as an absolute " +
                          "number of base pairs by giving an integer greater" +
                          " than 50. Not used by default")
    optional.add_argument("-z", "--serialize", dest='serialize',
                          action="store_true",
                          default=False,
                          help="if --serialize, runs seeding and assembly " +
                          "without multiprocessing. This is recommended for " +
                          "machines with less than 8GB RAM: %(default)s")
    optional.add_argument("--consensus", dest='initial_consensus',
                          action="store_true",
                          default=False,
                          help="if --initial_consensus, " +
                          "generate a mpileup-based consesnsus instead of " +
                          "doing a proper spades subassembly")
    optional.add_argument("--smalt_scoring", dest='smalt_scoring',
                          action="store",
                          default="match=1,subst=-4,gapopen=-4,gapext=-3",
                          help="if mapping with SMALT, " +
                          "submit custom smalt scoring via smalt -S " +
                          "scorespec option; default: %(default)s")
    optional.add_argument("--mapper_args", dest='mapper_args',
                          action="store",
                          default="-L 0,0 -U 0 -a",
                          help="submit custom parameters to mapper. " +
                          "And by mapper, I mean bwa, cause we dont support " +
                          "this option for SMALT, sorry. " +
                          "This requires knowledge of your chosen mapper's " +
                          "optional arguments. Proceed with caution!  " +
                          "default: %(default)s")
    # # TODO  Make these check a config file
    optional.add_argument("--spades_exe", dest="spades_exe",
                          action="store", default="spades.py",
                          help="Path to SPAdes executable; " +
                          "default: %(default)s")
    optional.add_argument("--samtools_exe", dest="samtools_exe",
                          action="store", default="samtools",
                          help="Path to samtools executable; " +
                          "default: %(default)s")
    optional.add_argument("--smalt_exe", dest="smalt_exe",
                          action="store", default="smalt",
                          help="Path to smalt executable;" +
                          " default: %(default)s")
    optional.add_argument("--bwa_exe", dest="bwa_exe",
                          action="store", default="bwa",
                          help="Path to BWA executable;" +
                          " default: %(default)s")
    optional.add_argument("--quast_exe", dest="quast_exe",
                          action="store", default="quast",
                          help="Path to quast executable; " +
                          "default: %(default)s")
    optional.add_argument("--bcftools_exe", dest="bcftools_exe",
                          action="store", default="bcftools",
                          help="Path to bcftools executable; " +
                          "default: %(default)s")
    optional.add_argument("-c", "--cores", dest='cores', action="store",
                          default=None, type=int,
                          help="cores to be used" +
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
    optional.add_argument("-m", "--memory", dest='memory', action="store",
                          default=8, type=int,
                          help="system memory available" +
                          "; default: %(default)s")
    optional.add_argument('--version', action='version',
                          version='riboSeed {version}'.format(
                              version=__version__))
    # # had to make this explicitly to call it a faux optional arg
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


def get_rec_from_generator(recordID, gen, method=None):
    """ given a record ID and and SeqIO generator return sequence of
    genbank record that has the loci, and call method to refresh generator
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


def get_smalt_full_install_cmds(smalt_exe, logger=None):  # pragma: no cover
    """ TODO replace this with swg tests for bambamc installation
    In the meantime, this looks for the files included with riboSeed
    (a bam file, reference, index, and fastq file), and generates the cmds
    to run a little test mapping
    """
    smalttestdir = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                "sample_data",
                                "smalt_test", "")
    assert logger is not None, "Must Use Logging"
    logger.debug("looking for smalt test dir: {0}".format(
        smalttestdir))
    if not os.path.exists(smalttestdir):
        raise FileNotFoundError(
            "Cannot find smalt_test dir containing " +
            "files to verify bambamc install! It should here: \n%s",
            smalttestdir)
    ref = os.path.join(smalttestdir, "ref_to_test_bambamc.fasta")
    index = os.path.join(smalttestdir, "test_index")
    test_bam = os.path.join(smalttestdir, "test_mapping.bam")
    test_reads = os.path.join(smalttestdir, "reads_to_test_bambamc.fastq")
    testindexcmd = str("{0} index {1} {2}".format(smalt_exe, index, ref))
    testmapcmd = str("{0} map -f bam -o {1} {2} {3}".format(smalt_exe,
                                                            test_bam,
                                                            index,
                                                            test_reads))
    return([testindexcmd, testmapcmd])


def test_smalt_bam_install(
        cmds, logger=None):  # pragma: no cover, cause pragma, no care
    """ using test data tha tcomes with package, ensure that
    the bambamc library was properly installed with SMALT instaltation
    """
    assert logger is not None, "must use logger"
    logger.info("testing instalation of SMALT and bambamc")
    smalttestdir = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                "sample_data",
                                "smalt_test", "")
    test_index = os.path.join(smalttestdir, "test_index")
    test_bam = os.path.join(smalttestdir, "test_mapping.bam")

    for i in cmds:
        try:
            logger.debug(i)
            subprocess.run([i],
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        except:
            logger.error(
                "Error running test to check bambamc lib is " +
                "installed! See github.com/gt1/bambamc " +
                "and the smalt install guide for more details." +
                "https://sourceforge.net/projects/smalt/files/")
            sys.exit(1)

    # remove the temp files
    os.remove(test_bam)
    os.remove(str(test_index + ".sma"))
    os.remove(str(test_index + ".smi"))


def check_fastqs_len_equal(file1, file2):
    """ using file_len from pyutilsnrw, check that the fastqs contain
    the same number of lines, ie tat the pairing looks proper.
    """
    if file_len(file1) != file_len(file2):
        raise ValueError(
            "Input Fastq's are of unequal length! Try " +
            "fixing with this script: " +
            "github.com/enormandeau/Scripts/fastqCombinePairedEnd.py")


def nonify_empty_lib_files(ngsLib, logger=None):
    # sometimes, if no singletons are found, we get an empty file.
    #  this shoudl weed out any empty read files before mapping, etc
    logger.info("checking for empty read files")
    EMPTIES = 0
    for f in ["readF", "readR", "readS0"]:
        # ignore if lib is None, as those wont be used anyway
        if getattr(ngsLib, f) is None:
            logger.debug("%s is set to None, and will be ignored", f)
            EMPTIES = EMPTIES + 1
            continue
        if not os.path.exists(getattr(ngsLib, f)):
            logger.warning("read file %s not found and can not be used " +
                           "for mapping!", f)
            EMPTIES = EMPTIES + 1
            # set to None so mapper will ignore
            setattr(ngsLib, f, None)
            continue
        # if lib is not none but file is of size 0
        logger.debug("size of %s: %f", getattr(ngsLib, f),
                     os.path.getsize(getattr(ngsLib, f)))
        if not os.path.getsize(getattr(ngsLib, f)) > 0:
            logger.warning("read file %s is empty and will not be used " +
                           "for mapping!", f)
            EMPTIES = EMPTIES + 1
            # set to None so mapper will ignore
            setattr(ngsLib, f, None)
    if EMPTIES == 3:
        raise ValueError("None of the read files hold data!")


def map_to_genome_ref_smalt(mapping_ob, ngsLib, cores,
                            samtools_exe, smalt_exe,
                            genome_fasta,
                            score_minimum=None,
                            scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
                            step=3, k=5, logger=None):  # pragma: no cover
    """run smalt based on pased args
    #TODO rework this to read libtype of ngslib object
    requires at least paired end input, but can handle an additional library
    of singleton reads. Will not work on just singletons
    """

    logger.info("Mapping reads to reference genome with SMALT")
    # check min score
    assert score_minimum is not None, "must sassign score outside map function!"
    score_min = score_minimum
    logger.debug(str("using a score min of " +
                     "{0}").format(score_min))
    # index the reference
    cmdindex = str("{0} index -k {1} -s {2} {3} {3}").format(
        smalt_exe, k, step, genome_fasta)
    # map paired end reads to reference index
    smaltcommands = [cmdindex]
    if "pe" in ngsLib.libtype:
        cmdmap = str('{0} map -l pe -S {1} ' +
                     '-m {2} -n {3} -g {4} -f bam -o {5} {6} {7} ' +
                     '{8}').format(smalt_exe, scoring,
                                   score_min, cores, ngsLib.smalt_dist_path,
                                   mapping_ob.pe_map_bam, genome_fasta,
                                   ngsLib.readF,
                                   ngsLib.readR)
        smaltcommands.append(cmdmap)
    else:
        with open(mapping_ob.pe_map_bam, 'w') as tempfile:
            tempfile.write("@HD riboseed_dummy_file")
        pass
    # if singletons are present, map those too.  Index is already made
    if ngsLib.readS0 is not None:  # and not ignore_singletons:
        # because erros are thrown if there is no file, this
        cmdmapS = str(
            "{0} map -S {1} -m {2} -n {3} -g {4} -f bam -o {5} " +
            "{6} {7}").format(smalt_exe, scoring, score_min, cores,
                              ngsLib.smalt_dist_path, mapping_ob.s_map_bam,
                              genome_fasta, ngsLib.readS0)
        with open(mapping_ob.s_map_bam, 'w') as tempfile:
            tempfile.write("@HD riboseed_dummy_file")
        # merge together the singleton and pe reads
        cmdmergeS = '{0} merge -f {3} {1} {2}'.format(
            samtools_exe, mapping_ob.pe_map_bam,
            mapping_ob.s_map_bam, mapping_ob.mapped_bam)
        smaltcommands.extend([cmdmapS, cmdmergeS])
    else:
        # if not already none, set to None when ignoring singleton
        ngsLib.readS0 = None
        # 'merge', but reallt just converts
        cmdmerge = str("{0} view -bh {1} >" +
                       "{2}").format(samtools_exe, mapping_ob.pe_map_bam, mapping_ob.mapped_bam)
        smaltcommands.extend([cmdmerge])
    logger.info("running SMALT:")
    logger.debug("with the following SMALT commands:")
    for i in smaltcommands:
        logger.debug(i)
        subprocess.run(i, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
    # report simgpleton reads mapped
    if ngsLib.readS0 is not None:
        logger.info(str("Singleton mapped reads: " +
                        get_number_mapped(mapping_ob.s_map_bam,
                                          samtools_exe=samtools_exe)))
    # report paired reads mapped
    if "pe" in ngsLib.libtype:
        logger.info(str("PE mapped reads: " +
                        get_number_mapped(mapping_ob.pe_map_bam,
                                          samtools_exe=samtools_exe)))
    combined_map_string = get_number_mapped(mapping_ob.mapped_bam_unfiltered,
                                            samtools_exe=samtools_exe)
    logger.info(str("Combined mapped reads: " + combined_map_string))
    # extract overall percentage as a float
    map_percentage = float(combined_map_string.split("(")[1].split("%")[0])

    # apparently there have been no errors, so mapping success!
    ngsLib.mapping_success = True
    return map_percentage


def index_sort_BAM(inbam):
    try:
        pysam.index(inbam)
    except pysam.utils.SamtoolsError:
        sorted_bam = os.path.join(os.path.dirname(inbam), "temp_sorted.bam")
        pysam.sort("-o", sorted_bam, inbam)
        inbam = sorted_bam
        try:
            pysam.index(inbam)
        except pysam.utils.SamtoolsError:
            raise ValueError("Your bam file is corrupted! No good")
    return inbam



def filter_bam_AS(inbam, outsam, score, logger=None):
    """ This is needed because bwa cannot filter based n alignment score
    for paired reads.
    https://sourceforge.net/p/bio-bwa/mailman/message/31968535/
    Given a  bam file from bwa (has "AS" tags), write out
    reads with AS higher than --score to outsam
    read count from https://www.biostars.org/p/1890/
    """
    notag = 0
    written = 0
    score_list = []
    inbam = index_sort_BAM(inbam)
    bam = pysam.AlignmentFile(inbam, "rb")
    osam = pysam.Samfile(outsam, 'wh', template=bam)
    for read in bam.fetch():
        if read.has_tag('AS'):
            score_list.append(read.get_tag('AS'))
            if read.get_tag('AS') >= score:
                osam.write(read)
                written = written + 1
            else:
                pass
        else:
            notag = notag + 1
            pass
    bam.close()
    logger.debug("Reads after filtering: %i", written)
    # if no reads pass the filtering score
    if written == 0:
        raise ValueError("No reads pass the filtering score! This commonly happens with " +
                         "short (<65bp) reads. Try rerunning with a decreased --score_min" +
                         " argument|")
    if notag != 0:
        logger.debug("Reads lacking alignment score: %i", notag)
    return score_list


def get_bam_AS(inbam, logger=None):
    """ Return the mappign scores for downstream QC plotting.
    """
    assert logger is not None, "must use logging"
    score_list = []
    count = 0
    try:
        pysam.index(inbam)
    except pysam.utils.SamtoolsError:
        raise ValueError("It looks like your bam file is unsorted! No good")
    bam = pysam.AlignmentFile(inbam, "rb")
    for read in bam.fetch():
        count = count + 1
        if read.has_tag('AS'):
            score_list.append(read.get_tag('AS'))
        else:
            pass
    bam.close()
    if len(score_list) != count:
        logger.warning("%i reads did not have AS tags",
                       count - len(score_list))
    return score_list


def convert_sam_to_bam(samtools_exe, bam, sam, reverse=False, logger=None):
    """
    becasue pysam doesnt like to write bams in an iterator, which makes sense

    """
    assert logger is not None, "must use logging"
    logger.debug("Converting with the following command:")
    if not reverse:
        cmd = "{0} view -o {1} -bS {2}".format(samtools_exe, bam, sam)
    else:
        cmd = "{0} view -o {2} -h {1}".format(samtools_exe, bam, sam)

    logger.debug(cmd)
    subprocess.run([cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def make_bwa_map_cmds(mapping_ob, ngsLib, cores,
                      samtools_exe, bwa_exe, genome_fasta,
                      add_args='-L 0,0 -U 0 -a', logger=None):
    """ make bwa sys commands. maps PE and S reads separately,
    then combines them into a X_mapped.bam file
    return a list of commands to run.
    """
    # index the reference
    cmdindex = str("{0} index {1}").format(
        bwa_exe, genome_fasta)
    # map paired end reads to reference index.
    bwacommands = [cmdindex]
    if "pe" in ngsLib.libtype:
        cmdmap = str('{0} mem -t {1} {2} -k 15 ' +
                     '{3} {4} {5} | {6} view -bh - | ' +
                     '{6} sort -o ' +
                     '{7} - ').format(bwa_exe,  # 0
                                      cores,  # 1
                                      add_args,  # 2
                                      genome_fasta,  # 3
                                      ngsLib.readF,  # 4
                                      ngsLib.readR,  # 5
                                      samtools_exe,  # 6
                                      mapping_ob.pe_map_bam)  # 7)
        bwacommands.append(cmdmap)
    else:
        assert ngsLib.readS0 is not None, \
            str("No readS0 attribute found, cannot run mapping with " +
                "any reads in .readS0 or .readF and .readR")

    # if singletons are present, map those too.  Index is already made
    if ngsLib.readS0 is not None:  # and not ignore_singletons:
        cmdmapS = str(
            '{0} mem -t {1} {2} -k 15 ' +
            '{3} {4} | {5} view -bh - | ' +
            '{5} sort -o {6} - ').format(bwa_exe,  # 0
                                         cores,  # 1
                                         add_args,  # 2
                                         genome_fasta,  # 3
                                         ngsLib.readS0,  # 4
                                         samtools_exe,  # 5
                                         mapping_ob.s_map_bam)  # 5)
        # merge together the singleton and pe reads, if there are any
        if "s_1" == ngsLib.libtype:
            cmdmergeS = str(
                "{0} view -bh {1} > {2}"
            ).format(samtools_exe, mapping_ob.s_map_bam, mapping_ob.mapped_bam_unfiltered)
        else:
            assert ngsLib.libtype == "pe_s", "error parsing libtype"
            cmdmergeS = '{0} merge -f {3} {1} {2}'.format(
                samtools_exe, mapping_ob.pe_map_bam,
                mapping_ob.s_map_bam, mapping_ob.mapped_bam_unfiltered)
        bwacommands.extend([cmdmapS, cmdmergeS])
    else:
        # if not already none, set to None when ignoring singleton
        ngsLib.readS0 = None
        cmdmerge = str("{0} view -bh {1} > " +
                       "{2}").format(samtools_exe, mapping_ob.pe_map_bam,
                                     mapping_ob.mapped_bam_unfiltered)
        bwacommands.extend([cmdmerge])
    return bwacommands


def map_to_genome_ref_bwa(mapping_ob, ngsLib, cores,
                          samtools_exe, bwa_exe, genome_fasta,
                          score_minimum=None,
                          add_args='-L 0,0 -U 0 -a', logger=None):
    """ Map to bam.  maps PE and S reads separately,
    then combines them into a X_mapped.bam file
    TODO:: break up into execution and comamnd generation
    """
    logger.info("Mapping reads to reference genome with BWA")
    bwacommands = make_bwa_map_cmds(mapping_ob, ngsLib, cores,
                                    samtools_exe, bwa_exe, genome_fasta,
                                    add_args=add_args, logger=None)
    logger.info("running BWA:")
    logger.debug("with the following BWA commands:")

    for i in bwacommands:
        logger.debug(i)
        subprocess.run(i, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
    # report simgpleton reads mapped
    if ngsLib.readS0 is not None:
        logger.info(str("Singleton mapped reads: " +
                        get_number_mapped(mapping_ob.s_map_bam,
                                          samtools_exe=samtools_exe)))
    # report paired reads mapped
    if "pe" in ngsLib.libtype:
        logger.info(str("PE mapped reads: " +
                        get_number_mapped(mapping_ob.pe_map_bam,
                                          samtools_exe=samtools_exe)))

    combined_map_string = get_number_mapped(mapping_ob.mapped_bam_unfiltered,
                                            samtools_exe=samtools_exe)
    logger.info(str("Combined mapped reads: " + combined_map_string))
    # extract overall percentage as a float
    try:
        map_perc_string = combined_map_string.split("(")[1].split("%")[0]
        map_percentage = float(map_perc_string)
    except ValueError:
        logger.error("Error mapping reads with bwa; line: %s:", map_perc_string)
        if map_perc_string == "N/A : N/A)":
            raise ValueError(
                "Error during mapping reads to reference; please examine " +
                "the sam/bam files. This could be to using a extremely " +
                "divergent reference or (more likely) using a library that "+
                " is not FR oriented. Consider running as a single library")
    # check min score
    if score_minimum is not None:
        score_min = score_minimum
    else:
        logger.debug(
            "no bwa mapping score min provided; default is 1/2 read " +
            "length or 50, whichever is greater.")
        # hard minimum of 50
        score_min = max(int(round(float(ngsLib.readlen) / 2.0)), 50)
    logger.debug("using a score minimum of %i", score_min)

    logger.debug("filtering mapped reads with an AS score minimum of %i",
                 score_min)
    score_list = filter_bam_AS(inbam=mapping_ob.mapped_bam_unfiltered,
                               outsam=mapping_ob.mapped_sam,
                               score=score_min, logger=logger)
    convert_sam_to_bam(
        bam=mapping_ob.mapped_bam,
        sam=mapping_ob.mapped_sam,
        samtools_exe=samtools_exe,
        reverse=False,
        logger=logger)
    logger.info(str("Mapped reads after filtering: " +
                    get_number_mapped(mapping_ob.mapped_bam,
                                      samtools_exe=samtools_exe)))
    # apparently there have been no errors, so mapping success!
    ngsLib.mapping_success = True
    return (map_percentage, score_list, score_min)


def convert_bam_to_fastqs_cmd(mapping_ob, ref_fasta, samtools_exe,
                              which='mapped', source_ext="_sam",
                              single=False, logger=None):
    """generate a cmd to convert a bam file to fastq, using samtools
    """
    assert which in ['mapped', 'unmapped'], \
        "only valid options are mapped and unmapped"
    read_path_dict = {'readF': None, 'readR': None, 'readS': None}
    logger.debug("preparing to convert extracted reads to make these files:")
    for key, value in read_path_dict.items():
        read_path_dict[key] = str(os.path.splitext(
            mapping_ob.mapped_bam)[0] + "_" + which + key + '.fastq')

    assert None not in read_path_dict.values(), \
        "Could not properly construct fastq names!"
    # if converting mapped reads, get them from the bam file
    if which == 'mapped':
        source_ext = '_bam'
    # else, leave the defaultsource ext (sam)
    else:
        pass

    if not single:
        samfastq = "{0} fastq {1} -1 {2} -2 {3} -s {4}".format(
            samtools_exe,
            getattr(mapping_ob, str(which + source_ext)),
            read_path_dict['readF'],
            read_path_dict['readR'],
            read_path_dict['readS'])
        for key in ['readF', 'readR', 'readS']:
            logger.debug(read_path_dict[key])

    else:
        # This option outputs all the reads in a single fastq
        # its needed for low coverage mappings when the F and R
        # file may end up empty.  Since default behaviour is to
        # treat F and R as single libraries anyway, this works
        samfastq = "{0} fastq {1} > {2} ".format(
            samtools_exe,
            getattr(mapping_ob, str(which + source_ext)),
            read_path_dict['readS'])
        logger.debug(read_path_dict['readS'])
        # Flag the others for ignoral
        # read_path_dict['readF'] = None
        # read_path_dict['readR'] = None
    return(samfastq, NgsLib(name=which, master=False,
                            logger=logger,
                            readF=read_path_dict['readF'],
                            readR=read_path_dict['readR'],
                            readS0=read_path_dict['readS'],
                            ref_fasta=ref_fasta))


# def check_version_from_cmd2(
#         exe,
#         cmd, line,
#         pattern=r"^__version__ = '(?P<version>[^']+)'$",
#         where='stderr',
#         min_version="0.0.0", logger=None,
#         coerce_two_digit=False):
#     """the guts have been stolen from pyani; returns version
#     from an system call that should return a version string.
#     Hacky, but better than nothing.
#     line arg is 1-indexed
#     .strip() is called on match to remove whitspaces

#     This will be removed when pyutilsnrw v0.1.1 is released
#     """
#     # exe_path = shutil.which(exe)
#     # if exe_path is None:
#     #     raise ValueError("executable %s not found!" % exe)
#     from distutils.version import StrictVersion
#     result = subprocess.run("{0} {1}".format(exe, cmd),
#                              # is this a securiy risk?
#                             shell=sys.platform != "win32",
#                             stdout=subprocess.PIPE,
#                             stderr=subprocess.PIPE,
#                             check=False)
#     logger.debug(result)
#     try:
#         if where == 'stderr':
#             printout = result.stderr.decode("utf-8").split("\n")
#         elif where == 'stdout':
#             printout = result.stdout.decode("utf-8").split("\n")
#         else:
#             raise ValueError("where option can only be 'stderr' or 'stdout'")
#     except Exception as e:
#         raise e
#     if logger:
#         logger.debug(printout)
#     this_version = None
#     try:
#         m = re.search(pattern, printout[line - 1])
#     except IndexError as e:
#         raise e
#     if m:
#         this_version = m.group('version').strip()
#     if logger:
#         logger.debug("this_version: %s", this_version)
#     if coerce_two_digit:
#         this_version = "0.{0}".format(this_version)
#         if logger:
#             logger.debug("coerced this_version: %s", this_version)
#     if this_version is None:
#         raise ValueError("No version was captured with pattern" +
#                          "{0}".format(pattern))
#     try:
#         if StrictVersion(this_version) < StrictVersion(min_version):
#             raise ValueError("{0} version {1} must be greater than {2}".format(
#                 cmd, this_version, min_version))
#     except Exception as e:
#         raise e
#     return this_version


def fiddle_with_spades_exe(spades_exe, logger=None):
    """  so heres the deal.  SPAdes 3.9  can be run with python3.5 and below.
     version 3.10 can be run with 3.6 and below.  If the version of spades
     and the version of execultion doest jive, this will try to correct it.
    return (wait for it, this is harebrained) the python executable needed!
    """
    assert logger is not None, "must use logging"
    assert sys.version_info[0] != 2 and sys.version_info[1] >= 5, \
        "how did we get here? cannot use riboSeed with anything less than 3.5"
    # spades will throw an error if run with the wrong version of python.  So
    # we will assume the error means that we are running with an incompatible
    # version of python
    # SPAdes 3.11.0 : all
    # SPAdes 3.10.0 : python3.5 and below
    # SPAdes 3.9.0 : python3.4 and below
    logger.debug("Making sure python's version is compatible with SPAdes")
    logger.debug("Python executable: %s", sys.executable)
    logger.debug("SPAdes executable: %s", spades_exe)
    try:
        spades_verison = check_version_from_cmd(
            exe=sys.executable + " " + spades_exe,
            cmd='--version', line=1, where='stderr',
            pattern=r"\s*v(?P<version>[^(]+)",
            min_version="3.9.0", logger=logger)
        SPADES_VERSION_SUCCESS = True
    except Exception as e:
        SPADES_VERSION_SUCCESS = False
        logger.debug("failed initial attempt to get spades version")
        logger.debug(e)
    if not SPADES_VERSION_SUCCESS:
        # if python 3.5 or 3.6, we should try to use python3.5 explicitly
        # if the user has it
        try:
            spades_verison = check_version_from_cmd(
                exe=shutil.which("python3.5") + " " + spades_exe,
                cmd='--version', line=1, where='stderr',
                pattern=r"\s*v(?P<version>[^(]+)",
                min_version="3.9.0", logger=logger)
            SPADES_VERSION_SUCCESS = True
            return shutil.which("python3.5")
        except Exception as e:
            SPADES_VERSION_SUCCESS = False
            logger.error(e)
            logger.error("There is an apparent mismatch between python" +
                         " and spades. Check to see if your spades " +
                         "version is compatible with your python version")
            sys.exit(1)

    return sys.executable


def generate_spades_cmd(
        mapping_ob, ngs_ob, ref_as_contig, python_exe, as_paired=True,
        addLibs="", prelim=False, k="21,33,55,77,99", spades_exe="spades.py",
        single_lib=False, logger=None, check_libs=False, check_exe=True):
    """return spades command so we can multiprocess the assemblies
    wrapper for common spades setting for long illumina reads
    ref_as_contig should be either None, 'trusted', or 'untrusted'
    prelim flag is True, only assembly is run, and without coverage corrections
    """
    assert logger is not None, "Must Use Logging"
    # make kmer comamnd empty if using "auto", or set to something like
    # -k 55,77,99
    if k is not 'auto':
        kmers = "-k " + k
    else:
        kmers = ""
    #  prepare reference, if being used
    if ref_as_contig is not None:
        alt_contig = "--{0}-contigs {1}".format(
            ref_as_contig, mapping_ob.ref_fasta)
    else:
        alt_contig = ''
    libs = []
    if single_lib:
        singles = "--pe1-s {0}".format(ngs_ob.readS0)
        pairs = ""
        libs.append(ngs_ob.readS0)
    elif as_paired and ngs_ob.readS0 is not None:  # for lib with both
        singles = "--pe1-s {0}".format(ngs_ob.readS0)
        pairs = "--pe1-1 {0} --pe1-2 {1} ".format(
            ngs_ob.readF, ngs_ob.readR)
        libs.append(ngs_ob.readS0)
        libs.append(ngs_ob.readF)
        libs.append(ngs_ob.readR)
    elif as_paired and ngs_ob.readS0 is None:  # for lib with just PE
        singles = ""
        pairs = "--pe1-1 {0} --pe1-2 {1}".format(
            ngs_ob.readF, ngs_ob.readR)
        libs.append(ngs_ob.readF)
        libs.append(ngs_ob.readR)
    # for libraries treating paired ends as two single-end libs
    elif not as_paired and ngs_ob.readS0 is None:
        singles = ''
        pairs = "--pe1-s {0} --pe2-s {1}".format(
            ngs_ob.readF, ngs_ob.readR)
        libs.append(ngs_ob.readF)
        libs.append(ngs_ob.readR)
    else:  # for 3 single end libraries
        singles = "--pe3-s {0} ".format(ngs_ob.readS0)
        pairs = str("--pe1-s {0} --pe2-s {1} ".format(
            ngs_ob.readF, ngs_ob.readR))
        libs.append(ngs_ob.readS0)
        libs.append(ngs_ob.readF)
        libs.append(ngs_ob.readR)
    reads = str(pairs + singles)

    if prelim:
        cmd = str(
            "{0} --only-assembler --cov-cutoff off --sc --careful {1} " +
            "{2} {3} {4} -o {5}"
        ).format(spades_exe, kmers, reads, alt_contig, addLibs,
                 mapping_ob.assembly_subdir)
    else:
        cmd = "{0} --careful {1} {2} {3} {4} -o {5}".format(
            spades_exe, kmers, reads, alt_contig, addLibs,
            mapping_ob.assembly_subdir)
    if check_libs:
        spades_cmd = make_spades_empty_check(liblist=libs, cmd=cmd,
                                             logger=logger)
    else:
        spades_cmd = cmd
    return python_exe + " " + spades_cmd


def make_spades_empty_check(liblist, cmd, logger):
    """ returns shell/spades cmd as string. All this does is make it a
    conditional shell cmd that depends on the presense of the file
    needed for assembly.  It is needed so we can bin all
    the cmds with multiprocessing.
    """
    logger.debug("constructing shell file check for subprocess cmd")
    prefix = "if "
    for i, lib in enumerate(liblist):
        if i != 0:
            prefix = prefix + "&& "
        check = "[ -s {0} ] ".format(lib)
        prefix = prefix + check
    suffix = str("; then {0} ; else echo 'input lib not found, " +
                 "skipping this SPAdes call' ; fi").format(cmd)
    return str(prefix + suffix)


def exclude_subassembly_based_on_coverage(clu, iteration, logger=None):
    """ if using coverage_exclusion, then return 0 if passing or 2 if not.
    if not, return None
    deal with those excluded from assembly by lack of coverage depth
      ie  (coverage_exclusion=True)
    """
    assert logger is not None, "must use logging"
    if clu.coverage_exclusion is not None:
        assert clu.coverage_exclusion, \
            "this should only be set by the partition_mapping method."
        logger.warning("THIS FEATURE IS NOT WELL TESTED! USE WITH TREPIDATION")
        #  if this is the first iteration, return two
        # Otherwise,  UPDATED: continue with warning
        if iteration > 0:
            clu.mappings[iteration].assembled_contig = \
                clu.mappings[iteration - 1].assembled_contig
            logger.warning(" the coverage is worryingly low, but we " +
                           "will continue.")
            return 0
        else:
            logger.info("the coverage is worryingly low, and as this " +
                        "is the first iteration, we must discard this contig")
            return 2
    else:
        return None

def evaluate_spades_success(clu, mapping_ob, proceed_to_target, target_len,
                            include_short_contigs, min_assembly_len,
                            flank=1000,
                            min_delta=10,
                            keep_best_contig=True,
                            seqname='', logger=None):
    """return success codes:
    0 = include contigs, all good
    DEPRECIATED! 1 = include contigs, but dont keep iterating
    2 = exclude contigs, and keep from iterating
    3 = exclude contigs, error ocurred
    """
    # DANGEROUS_CONTIG_LENGTH_THRESHOLD_FACTOR = 6
    prelog = "{0}-{1}-iter-{2}:".format("SEED_cluster", clu.index,
                                        mapping_ob.iteration)
    assert logger is not None, "Must Use Logging"
    if seqname == '':
        seqname = os.path.splitext(os.path.basename(mapping_ob.ref_fasta))[0]
    cov_exclude_result = exclude_subassembly_based_on_coverage(
        clu=clu, iteration=mapping_ob.iteration, logger=logger)
    if cov_exclude_result is not None:
        return cov_exclude_result
    mapping_ob.assembled_contig = os.path.join(
        mapping_ob.assembly_subdir, "contigs.fasta")
    logger.debug("checking for the following file: \n{0}".format(
        mapping_ob.assembled_contig))
    # check for a spades failure
    if not (os.path.isfile(mapping_ob.assembled_contig) and
            os.path.getsize(mapping_ob.assembled_contig) > 0):
        logger.warning(
            "%s No output from SPAdes this time around! return code 3",
            prelog)
        return 3
    # by default, we keep only the longest, bestest, most fantastic-est contig
    if keep_best_contig:
        logger.debug("reserving first contig")
        try:
            keep_only_first_contig(
                os.path.join(mapping_ob.assembly_subdir, "contigs.fasta"),
                newname=seqname)
        except Exception as f:
            logger.error(f)
            raise f
    # -------------------------- --------------------------- #

    logger.info("%s analyzing  mapping", prelog)
    seed_len = get_fasta_lengths(clu.mappings[0].ref_fasta)[0]
    # seed_len = get_fasta_lengths(mapping_ob.ref_fasta)[0]
    # set proceed_to_target params
    if proceed_to_target:
        if target_len > 0 and 5 > target_len:
            target_seed_len = int(target_len * seed_len)
        elif target_len > 50:
            target_seed_len = int(target_len)
        else:
            logger.error("%s invalid target length provided; must be given " +
                         "as fraction of total length or as an absolute " +
                         "number of base pairs greater than 50", prelog)
            sys.exit(1)
    else:
        pass
    # compare lengths of reference and freshly assembled contig
    contig_len = get_fasta_lengths(mapping_ob.assembled_contig)[0]
    ref_len = get_fasta_lengths(mapping_ob.ref_fasta)[0]
    contig_length_diff = contig_len - ref_len
    logger.info("%s Seed length: %i", prelog, seed_len)
    if proceed_to_target:
        logger.info("Target length: {0}".format(target_seed_len))
    logger.info("%s Length of this iteration's longest contig: %i",
                prelog, contig_len)
    if mapping_ob.iteration != 0:
        logger.info("%s Length of previous longest contig: %i",
                    prelog, ref_len)
        logger.info("%s The new contig differs from the previous " +
                    "iteration by %i bases", prelog, contig_length_diff)
    else:
        logger.info("%s The new contig differs from the reference " +
                    "seed by %i bases", prelog, contig_length_diff)
    if contig_len > (ref_len + (2 * flank)):
        logger.warning(
            "Contig length is exceedingly long!  We set the threshold of " +
            "twice the flanking length as the maximum allowed long-read " +
            "length. This may indicate  bad mapping parameters, so the " +
            "long-read will be discarded.  Return code 2")
        return 2

    # This cuts failing assemblies short
    if min_assembly_len > contig_len:
        logger.warning("The first iteration's assembly's best contig " +
                       "is not greater than length set by " +
                       "--min_assembly_len. Assembly will likely fail if " +
                       "the contig does not meet the length of the seed")
        # if mapping_ob.iteration > 0:
        if include_short_contigs:
            logger.warning("Continuing to , but if this occurs for more " +
                           "than one seed, we reccommend  you abort and " +
                           "retry with longer seeds, a different ref, " +
                           "or re-examine the riboSnag clustering")
            logger.warning("Return code 0")
            return 0
        else:
            logger.warning("Return code 2")
            return 2
    elif proceed_to_target and contig_len >= target_seed_len:
        logger.info("target length threshold! has been reached; " +
                    "skipping future iterations.  return code 0")
        return 0
    # if not first time through, ensure adequate change between iterations to
    # avoid problems with trying to assemble a very small number of reads
    elif min_delta > abs(contig_length_diff) and mapping_ob.iteration != 0:
        logger.warning(str(
            "The length of the assembled contig didn't change more " +
            "more than {0}bp between rounds of iteration. Continuing " +
            "will likely cause error; skipping future iterations. " +
            "return code 0").format(min_delta))
        return 0
    else:
        logger.debug("return code 0")
        return 0


def parse_subassembly_return_code(cluster, logger=None):
    """ given a return code from the above spades success function,
    set object attributes as needed
    20170531 depreciated return code 1: as we have moved to using the
    whole library for mapping, not just the unmapped reads, it is more
    important to keep in even unchanging contigs to allow for
    competition during mapping

    ----------------------
    return success codes:
    0 = include contigs, all good
    1 = include contigs, but dont keep iterating
    2 = exclude contigs, and keep from iterating
    3 = exclude contigs, error ocurred
    """
    assert logger is not None, "must use logging"
    if cluster.assembly_success == 3:
        # TODO other error handling; make a "failed" counter?
        cluster.continue_iterating = False
        cluster.keep_contigs = False
    elif cluster.assembly_success == 2:
        cluster.continue_iterating = False
        cluster.keep_contigs = False
    elif cluster.assembly_success == 1:
        raise ValueError("return code 1 depreciated! a warning can be " +
                         "issued for short asssemblies, but they must " +
                         "remain in the pseudogenome")
    elif cluster.assembly_success == 0:
        cluster.continue_iterating = True
        cluster.keep_contigs = True
    else:
        raise ValueError("Error evaluating spades results return!")


def make_quick_quast_table(pathlist, write=False, writedir=None, logger=None):
    """ given paths to two or more quast reports, this generates dictionary
    where the key is the field in the report and the value is a list of
    the values for each report.   Hand for passing to the logger function.
    This skips any fields not in first report, for better or worse...
    """
    assert logger is not None, "Must Use Logging"
    assert isinstance(pathlist, list) is True,\
        "paths for quast reports must be in a list!"
    filelist = pathlist
    logger.debug("Quast reports to combine: %s", str(filelist))
    mainDict = {}
    counter = 0
    for i in filelist:
        if counter == 0:
            try:
                with open(i, "r") as handle:
                    for dex, line in enumerate(handle):
                        row, val = line.strip().split("\t")
                        if dex in [0]:
                            continue  # skip header
                        else:
                            mainDict[row] = [val]
            except Exception:
                raise ValueError("error parsing %s", i)
        else:
            report_list = []
            try:
                with open(i, "r") as handle:
                    for dex, line in enumerate(handle):
                        row, val = line.strip().split("\t")
                        report_list.append([row, val])
                    # logger.debug("report list: %s", str(report_list))
                    for k, v in mainDict.items():
                        if k in [x[0] for x in report_list]:
                            mainDict[k].append(
                                str([x[1] for x in
                                     report_list if x[0] == k][0]))
                        else:
                            mainDict[k].append("XX")
            except Exception as e:
                logger.warning("error parsing %s", i)
                raise e
        counter = counter + 1
    if write:
        if writedir is None:
            logger.warning("no output dir, cannot write!")
            return mainDict
        try:
            with open(os.path.join(writedir, "combined_quast_report.tsv"),
                      "w") as outfile:
                for k, v in sorted(mainDict.items()):
                    # logger.debug("{0}\t{1}\n".format(k, str("\t".join(v))))
                    outfile.write("{0}\t{1}\n".format(
                        str(k), str("\t".join(v))))
        except Exception as e:
            raise e
    return mainDict


def check_kmer_vs_reads(k, readlen, min_diff=2, logger=None):
    assert logger is not None, "must use logging! "
    # ignore this if we are letting spades
    if k is 'auto':
        logger.debug("no need to check k, we let spades set k")
        return k
    try:
        klist = [int(x) for x in k.split(",")]
    except Exception as e:
        logger.error("error splitting kmers by comma!")
        logger.error(e)
        logger.error(last_exception())
        raise ValueError
    logger.debug(klist)
    new_ks = []
    for i in klist:
        if i > readlen:
            logger.warning("removing %d from list of kmers: exceeds read length",
                        i)
        elif readlen - i <= min_diff:
            logger.warning("removing %d from list of kmers: too close " +
                        "to read length", i)
        elif i % 2 == 0:
            logger.warning("removing %d from list of kmers: must be odd", i)
        else:
            new_ks.append(i)
    return ",".join([str(x) for x in new_ks])


def make_samtools_depth_cmds(exe, bam, chrom, start, end, region=None, prep=False):
    """ this just makes the commands to get the depth from samtools.
    If prep, the sam file gets sorted and indexed first
    """
    prep_cmds = []
    # cmd = "samtools depth ./iter_1_s_mappi.bam -r scannedScaffolds:5000-6000"
    sorted_bam = os.path.join(
        os.path.dirname(bam),
        str(os.path.splitext(os.path.basename(bam))[0] + "_sorted.bam"))
    # sort that bam, just in case
    prep_cmds.append(str("{0} sort {1} > {2}").format(exe, bam, sorted_bam))
    # index that bam!
    prep_cmds.append(str("{0} index {1}").format(exe, sorted_bam))
    if prep:
        bamfile = sorted_bam
    else:
        bamfile = bam
    # extract the depth stats for a region
    if region is None:
        depth_cmd = str("{0} depth -r {2}:{3}-{4} {1}").format(
            exe, bamfile, chrom, start, end)
    else:
        depth_cmd = str("{0} depth -r {2} {1}").format(
            exe, bamfile, region)
    return (prep_cmds, depth_cmd)


def parse_samtools_depth_results(result, logger=None):
    """ parses out the subprocess results from samtools depth
    """
    assert logger is not None, "must use logging"
    try:
        splits = result.stdout.decode("utf-8").split("\n")[0].split("\t")
        if len(splits) != 3:
            logger.warning("unable to split the results from samtools depth")
        else:
            pass
    except Exception as e:
        raise e
    covs = [int(x.split("\t")[2]) for
            x in result.stdout.decode("utf-8").split("\n")[0: -1]]
    if len(covs) == 0:
        if result.returncode != 0:
            logger.warning("Error parsing samtools depth results! " +
                           "Here are the results:")
            logger.warning(result)
            logger.warning("This isn't always fatal, so we will continue. " +
                           "but take a look with IGB or similar so there " +
                           "aren't any suprises down the road.")
        return [[""], 0]

    average = float(sum(covs)) / float(len(covs))
    return [covs, average]


def get_samtools_depths(samtools_exe, bam, chrom, start, end,
                        prep=False, region=None, logger=None):
    """ Use samtools depth and awk to get the average coverage depth of a
    particular region
    """
    prep_cmds, depth_cmd = make_samtools_depth_cmds(
        exe=samtools_exe, bam=bam, chrom=chrom,
        start=start, end=end, region=region, prep=prep)
    logger.debug("running the following commands to get coverage:")
    if prep:
        for i in prep_cmds:  # index and sort
            logger.debug(i)
            subprocess.run(i,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
    else:
        pass
    # get the results from the depth call
    logger.debug(depth_cmd)
    result = subprocess.run(depth_cmd,
                            shell=sys.platform != "win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=False)
    covs, ave = parse_samtools_depth_results(result, logger)
    return [covs, ave]


def prepare_next_mapping(cluster, seedGenome, flank,
                         logger=None):
    """use within partition mapping funtion;
    makes LociMapping, get region coords, write extracted region,
    """
    mapping_subdir = os.path.join(
        seedGenome.output_root, cluster.cluster_dir_name,
        "{0}_cluster_{1}_mapping_iteration_{2}".format(
            cluster.sequence_id, cluster.index, seedGenome.this_iteration))
    assembly_subdir = os.path.join(
        seedGenome.output_root, cluster.cluster_dir_name,
        "{0}_cluster_{1}_assembly_iteration_{2}".format(
            cluster.sequence_id, cluster.index, seedGenome.this_iteration))

    mapping0 = LociMapping(
        name="{0}_cluster_{1}".format(
            cluster.sequence_id, cluster.index),
        iteration=seedGenome.this_iteration,
        assembly_subdir_needed=True,
        mapping_subdir=mapping_subdir,
        assembly_subdir=assembly_subdir)
    # if first time through, get the global start and end coords.
    if cluster.global_start_coord is None or cluster.global_end_coord is None:
        if seedGenome.this_iteration != 0:
            raise ValueError(
                "global start and end should be defined previously! Exiting")
        if sorted([x.start_coord for x in cluster.loci_list]) != \
           [x.start_coord for x in cluster.loci_list]:
            logger.warning("Coords are not in increasing order; " +
                           "you've been warned")
        start_list = sorted([x.start_coord for x in cluster.loci_list])
        logger.debug("Start_list: {0}".format(start_list))

        logger.debug("Finding coords to gather reads from the following loci:")
        for i in cluster.loci_list:
            logger.debug("%s cluster %i -- locus %i -- %s (%i, %i)(%i) %s",
                         i.sequence_id, cluster.index,
                         i.index, i.locus_tag,
                         i.start_coord, i.end_coord, i.strand,
                         i.product)
        #  This works as long as coords are never in reverse order
        cluster.global_start_coord = min([x.start_coord for
                                          x in cluster.loci_list]) - flank
        # if start is negative, just use 1, the beginning of the sequence
        if cluster.global_start_coord < 1:
            logger.warning(
                "Caution! Cannot retrieve full flanking region, as " +
                "the 5' flanking region extends past start of " +
                "sequence. If this is a problem, try using a smaller " +
                "--flanking region, and/or if  appropriate, run with " +
                "--linear.")
            cluster.global_start_coord = 1
        cluster.global_end_coord = max([x.end_coord for
                                        x in cluster.loci_list]) + flank
        # logger.debug("rec len: %i", len(cluster.seq_record.seq))
        if cluster.global_end_coord > len(cluster.seq_record):
            logger.warning(
                "Caution! Cannot retrieve full flanking region, as " +
                "the 5' flanking region extends past start of " +
                "sequence. If this is a problem, try using a smaller " +
                "--flanking region, and/or if  appropriate, run with " +
                "--linear.")
            cluster.global_end_coord = len(cluster.seq_record)
        logger.debug("global start and end: %s %s",
                     cluster.global_start_coord,
                     cluster.global_end_coord)
    #  if not the first time though, fuhgetaboudit.
    #  Ie, the coords have been reassigned by the make_faux_genome function
    #  WE WONT USE A FLANKING REGION BECAUSE NO FLANKING READS ARE AVAILIBLE!
    #  meaning, the overhang is gained from the bits that overhand the end of
    #  the mapping. Because both SMALT and BWA use soft-clipping by defualt, we
    #  recover and use the clipped regions
    else:
        logger.info("using coords from previous iterations 'genome'")
    logger.debug("Coordinates for %s cluster %i:  [%i - %i]",
                cluster.seq_record.id,
                cluster.index,
                cluster.global_start_coord,
                cluster.global_end_coord)

    cluster.extractedSeqRecord = SeqRecord(
        cluster.seq_record.seq[
            cluster.global_start_coord:
            cluster.global_end_coord])

    mapping0.ref_fasta = os.path.join(mapping0.mapping_subdir,
                                      "extracted_seed_sequence.fasta")
    with open(mapping0.ref_fasta, "w") as writepath:
        SeqIO.write(cluster.extractedSeqRecord, writepath, 'fasta')
    cluster.mappings.append(mapping0)


def make_mapped_partition_cmds(cluster, mapping_ob, seedGenome, samtools_exe):
    """ returns cmds and region
    """
    # Prepare for partitioning
    partition_cmds = []
    # sort our source bam
    sort_cmd = str("{0} sort {1} > {2}").format(
        samtools_exe,
        seedGenome.iter_mapping_list[seedGenome.this_iteration].mapped_bam,
        seedGenome.iter_mapping_list[
            seedGenome.this_iteration].sorted_mapped_bam)
    # index it
    index_cmd = str("{0} index {1}").format(
        samtools_exe, seedGenome.iter_mapping_list[
            seedGenome.this_iteration].sorted_mapped_bam)
    partition_cmds.extend([sort_cmd, index_cmd])
    # define the region to extract
    region_to_extract = "{0}:{1}-{2}".format(
        cluster.sequence_id, cluster.global_start_coord,
        cluster.global_end_coord)
    # make a subser from of reads in that region
    view_cmd = str("{0} view -o {1} {2} {3}").format(
        samtools_exe, mapping_ob.mapped_bam,
        seedGenome.iter_mapping_list[
            seedGenome.this_iteration].sorted_mapped_bam,
        region_to_extract)
    partition_cmds.append(view_cmd)
    return (partition_cmds, region_to_extract)


def make_unmapped_partition_cmds(mapped_regions, samtools_exe, seedGenome):
    """ given a list of regions (formatted for samtools view, etc) make a
    list of mapped reads (file path stored under mapped_ids_txt), and
    use the cgrep voodoo to make a sam file from the full library without
    the mapped reads. returns a cmd as a string
    """
    # starting at second iteration, copy previous iterms mapped_ids_txt
    # as a starting point so we can track the reads better.
    unmapped_cmds = []
    if seedGenome.this_iteration > 0:
        shutil.copyfile(
            seedGenome.iter_mapping_list[
                seedGenome.this_iteration - 1].mapped_ids_txt,
            seedGenome.iter_mapping_list[
                seedGenome.this_iteration].mapped_ids_txt)

    # moved to filter_bam_as function
    make_mapped_sam = "{0} view -o {1} -h {2}".format(
        samtools_exe,
        seedGenome.iter_mapping_list[seedGenome.this_iteration].mapped_sam,
        seedGenome.iter_mapping_list[seedGenome.this_iteration].mapped_bam)

    unmapped_cmds.append(make_mapped_sam)

    # for each region, add read names in that region to
    # a list (taken from previous iteration if there has been one)
    for region in mapped_regions:
        unmapped_cmds.append(
            "{0} view {1} {2} | cut -f1 >> {3}".format(
                samtools_exe,
                seedGenome.iter_mapping_list[
                    seedGenome.this_iteration].sorted_mapped_bam,
                region,
                seedGenome.iter_mapping_list[
                    seedGenome.this_iteration].mapped_ids_txt))
    uniquify_list = "sort -u {0} -o {0}".format(
        seedGenome.iter_mapping_list[seedGenome.this_iteration].mapped_ids_txt)
    unmapped_cmds.append(uniquify_list)
    return unmapped_cmds


def pysam_extract_reads(sam, textfile, unmapped_sam, logger=None):
    """ This replaces the "LC_ALL " grep call from the above function. On macs,
    there is no speedup gained.
    """
    qfile = textfile  # contains read names of mapped reads
    sfile = sam  # mapping of all reads to genome
    ofile = unmapped_sam  # destination sam of all reads not in regions
    nunmapped = 0
    total = 0
    # Load query fixed strings as a set
    with open(qfile, 'r') as qfh:
        queries = {q.strip() for q in qfh.readlines() if len(q.strip())}
    # Subset reads
    samfile = pysam.AlignmentFile(sfile, 'r')
    osam = pysam.Samfile(ofile, 'wh', template=samfile)
    for read in samfile.fetch():
        total = total + 1
        if read.qname in queries:
            continue
        else:
            # write out read names not in textfile
            nunmapped = nunmapped + 1
            osam.write(read)
    if logger:
        logger.info("Wrote %i unmapped reads of the %i total from %s to %s",
                    nunmapped, total, sam, ofile)
    osam.close()


def partition_mapping(seedGenome, samtools_exe, flank, min_flank_depth,
                      cluster_list=None, logger=None):
    """ Extract interesting stuff based on coords, not a binary
    mapped/not_mapped condition
    Also, if min_flanking_depth, mark reads with low
    mapping coverage for exclusion
    """
    mapped_regions = []
    logger.info("processing mapping for iteration %i",
                seedGenome.this_iteration)
    for cluster in cluster_list:
        prepare_next_mapping(cluster=cluster, seedGenome=seedGenome,
                             flank=flank, logger=logger)

    mapped_regions = []
    all_depths = []  # each entry is a tuple (idx, start_ave, end_ave)
    filtered_cluster_list = []
    for cluster in cluster_list:
        mapped_partition_cmds, reg_to_extract = make_mapped_partition_cmds(
            cluster=cluster, mapping_ob=cluster.mappings[-1],
            seedGenome=seedGenome, samtools_exe=samtools_exe)
        for cmd in mapped_partition_cmds:
            logger.debug(cmd)
            subprocess.run([cmd],
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        start_depths, start_ave_depth = get_samtools_depths(
            bam=seedGenome.iter_mapping_list[
                seedGenome.this_iteration].sorted_mapped_bam,
            chrom=cluster.sequence_id,
            start=cluster.global_start_coord,
            end=cluster.global_start_coord + flank,
            region=None,
            prep=False,
            samtools_exe=samtools_exe,
            logger=logger)
        end_depths, end_ave_depth = get_samtools_depths(
            bam=seedGenome.iter_mapping_list[
                seedGenome.this_iteration].sorted_mapped_bam,
            chrom=cluster.sequence_id,
            start=cluster.global_end_coord - flank,
            end=cluster.global_end_coord,
            region=None,
            prep=False,
            samtools_exe=samtools_exe,
            logger=logger)
        logger.info("Coverage for cluster " +
                    "%i:\n\t5' %ibp-region: %.2f \n\t3' %ibp-region: %.2f",
                    cluster.index,
                    flank,
                    start_ave_depth,
                    flank,
                    end_ave_depth)
        if start_ave_depth < min_flank_depth:
            logger.warning(str("cluster {0} has insufficient 5' flanking " +
                               "coverage depth for subassembly, and will be " +
                               "removed").format(cluster.index))
            cluster.coverage_exclusion = True
        elif end_ave_depth < min_flank_depth:
            logger.warning(str("cluster {0} has insufficient 3' flanking " +
                               "coverage depth for subassembly, and will be " +
                               "removed").format(cluster.index))
            cluster.coverage_exclusion = True
        else:
            mapped_regions.append(reg_to_extract)
            filtered_cluster_list.append(cluster)
        # regardsless, report stats here
        all_depths.append((cluster.index, start_ave_depth, end_ave_depth))
    logger.info("mapped regions for iteration %i:\n%s",
                seedGenome.this_iteration,
                "\n".join([x for x in mapped_regions]))
    # make commands to extract all the reads NOT mapping to the rDNA regions
    unmapped_partition_cmds = make_unmapped_partition_cmds(
        mapped_regions=mapped_regions, samtools_exe=samtools_exe,
        seedGenome=seedGenome)
    for cmd in unmapped_partition_cmds:
        logger.debug(cmd)
        subprocess.run([cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    logger.info("using pysam to extract a subset of reads ")
    # this may look wierd: for iteration 0, we extract from the mapping.
    # for each one after that, we extract from the previous mapping.
    unmapped_reads_index = 0 if seedGenome.this_iteration == 0 \
        else seedGenome.this_iteration - 1
    pysam_extract_reads(
        sam=seedGenome.iter_mapping_list[unmapped_reads_index].mapped_sam,
        textfile=seedGenome.iter_mapping_list[
            seedGenome.this_iteration].mapped_ids_txt,
        unmapped_sam=seedGenome.iter_mapping_list[
            seedGenome.this_iteration].unmapped_sam, logger=logger)
    # sam_score_list = get_sam_AS
    return (all_depths, filtered_cluster_list)


def add_coords_to_clusters(seedGenome, logger=None):
    """ given a genbank file and some locus tags, add the coordinates, etc,
    to the entry in the seed Genome
    """
    for cluster in seedGenome.loci_clusters:  # for each cluster of loci
        # get seq record that cluster is  from
        try:
            cluster.seq_record = get_rec_from_generator(
                recordID=cluster.sequence_id,
                gen=seedGenome.seq_records,
                method=seedGenome.refresh_seq_rec_generator)
        except Exception as e:
            raise e
        try:  # make coord list
            extract_coords_from_locus(
                cluster=cluster, feature=cluster.feat_of_interest,
                logger=logger)
        except Exception as e:
            raise e
        logger.debug("Here are the detected region,coords, strand, product, " +
                     "locus tag, subfeatures and sequence id of the results:")
        logger.debug(str(cluster.__dict__))


def bool_run_quast(quast_exe, logger):
    if sys.version_info.minor == 6:
        logger.warning("QUAST only supports python3.5 and below. We are " +
                       "skipping QUAST evalutation")
        return False
    if quast_exe is None:
        logger.warning("QUAST executable not provided or available. We are " +
                       "skipping QUAST evalutation")
        return False
    try:
        quast_version = check_version_from_cmd(
            exe=sys.executable + " " + quast_exe,
            cmd='--version', line=1, where='stdout',
            pattern=r".* v(?P<version>[^\n,]+)",
            min_version="4.0", logger=logger,
            coerce_two_digit=False)
        if quast_version =="4.5":
            logger.warning("Due to bugs in QUAST 4.5, we will not run QUAST")
            return False

    except Exception as e:
        logger.error(e)
        logger.warning("Error occured while trying to check QUAST version." +
                       "We are skipping QUAST evalutation")
        return False
    return True


def make_quast_command(exes, output_root, ref, assembly_subdir, name,
                       logger=None):
    assert logger is not None, "must use logging"
    quast_cmd = str("{0} {1} {2} {3} -o {4}").format(
        exes.python,
        exes.quast,
        ref,
        os.path.join(assembly_subdir, "contigs.fasta"),
        os.path.join(output_root, str("quast_" + name)))
    return quast_cmd


def get_final_assemblies_cmds(seedGenome, exes,
                              ref_as_contig,
                              additional_libs,
                              cores,
                              memory,
                              serialize,
                              skip_control=True,
                              kmers="21,33,55,77,99", logger=None):
    """make cmds for runnning of SPAdes and QUAST final assembly and analysis.
    if skip_control, just do the de fere novo assembly.  otherwise, do bother
    returns list of listed cmds
    ([[spades_cmd, quast_cmd], [spades_cmd2, quast_cmd2]])
    """
    logger.info("\n\nStarting Final Assemblies\n\n")
    quast_reports = []
    cmd_list = []
    final_list = ["de_fere_novo"]
    if not skip_control:
        final_list.append("de_novo")
    for j in final_list:
        final_mapping = LociMapping(
            iteration=0,
            name=j,
            mapping_subdir=os.path.join(
                seedGenome.output_root,
                "final_{0}_mapping".format(j)),
            assembly_subdir_needed=True,
            assembly_subdir=os.path.join(
                seedGenome.output_root,
                "final_{0}_assembly".format(j)))
        # logger.info("\n\nRunning %s SPAdes \n" % j)
        if j == "de_novo":
            final_mapping.ref_fasta = ''
            assembly_ref_as_contig = None
        else:
            assert j == "de_fere_novo", \
                "Only valid cases are de novo and de fere novo!"
            final_mapping.ref_fasta = seedGenome.assembled_seeds
            assembly_ref_as_contig = ref_as_contig

        # remove unneeded dir
        os.rmdir(final_mapping.mapping_subdir)

        logger.info("Getting commands for %s SPAdes" % j)
        spades_cmd = generate_spades_cmd(
            single_lib=seedGenome.master_ngs_ob.libtype == "s_1",
            check_libs=False,
            python_exe=exes.python,
            mapping_ob=final_mapping, ngs_ob=seedGenome.master_ngs_ob,
            ref_as_contig=assembly_ref_as_contig, as_paired=True, prelim=False,
            k=kmers, spades_exe=exes.spades, logger=logger)
        modest_spades_cmd = make_modest_spades_cmd(
            cmd=spades_cmd, cores=cores, memory=memory, split=2,
            serialize=serialize, logger=logger)
        ## add additional cmdline args for assembly
        if additional_libs is not None:
            modest_spades_cmd = "{0} {1}".format(
                modest_spades_cmd, additional_libs)
        ref = str("-R %s" % seedGenome.ref_fasta)
        quast_cmd = make_quast_command(
            exes=exes, output_root=seedGenome.output_root, ref=ref,
            assembly_subdir=final_mapping.assembly_subdir, name=j,
            logger=logger)
        if bool_run_quast(exes.quast, logger):
            quast_reports.append(
                os.path.join(seedGenome.output_root,
                             str("quast_" + j), "report.tsv"))
            cmd_list.append([modest_spades_cmd, quast_cmd])
        else:
            quast_reports = None
            cmd_list.append([modest_spades_cmd, None])
    return(cmd_list, quast_reports)


def check_spades_extra_library_input(inp):
    """ check user supplied args against SPAdes args.
     This is pretty brittle, cause its experts only for now.
    will raise an error if any funny business is detected
    It will not however check your syntax, ensure you have correct pairs, etc.
    You are in the deep end now.

    """
    scary_chars = [";", "|"]
    for char in scary_chars:
        if char in inp:
            raise ValueError("'%s' cannot be part of lib input" %char)

    single_valid = "--s"
    simple_valid = ["--sanger", "--pacbio", "--nanopore"]
    comp_valid_prefixes = ["--pe", "--mp", "--hqmp", "--nxmate"]
    comp_valid_suffixes = ["1", "12", "2", "s", "ff", "fr", "rf"]
    splitcmds = inp.split(" ")
    cmds_dict = {}

    MARKER = True
    # attempt to separate out args into key value pairs.
    try:
        for i, cmd in enumerate(splitcmds):
            if MARKER:
                cmds_dict[cmd] = splitcmds[i + 1]
            MARKER = not MARKER
    except IndexError:
        raise IndexError("All args must have exactly 1 value; no flags!")

    for cmd, value in cmds_dict.items():
        VALID = False
        # check if its a single library (coded as --s#)
        if cmd[:-1] == single_valid:
            # does it have a valid numeric library idenifier?
            for i in range(1,10):
                if cmd == single_valid + str(i):
                    VALID = True
        # check if its a simple arg
        elif cmd in simple_valid:
            VALID = True

        # check if its a compound arg
        else:
            # does it have a valid prefix?
            for pref in comp_valid_prefixes:
                if cmd.startswith(pref):
                    # does it have a valid numeric library idenifier?
                    for i in range(1,10):
                        if cmd.startswith(pref + str(i)):
                            # does it have a valid suffix?
                            for suff in comp_valid_suffixes:
                                if cmd == (pref + str(i) + "-" + suff):
                                    VALID = True
                                    break
        if not VALID:
            raise ValueError("Invalid spades arg %s" % cmd)





def make_faux_genome(cluster_list, seedGenome, iteration,
                     output_root, nbuff, logger=None):
    """ stictch together viable assembled contigs.  perhaps more importnatly,
    this also re-write thes coords relative to the new "genome"
    returns path to new faux_genome
    """
    logger.info("preparing extracted region genome for next round of mapping")
    logger.debug("using %i sequences", len(cluster_list))
    nbuffer = "N" * nbuff
    faux_genome = ""
    counter = 0
    new_seq_name = seedGenome.name
    if len(cluster_list) == 0:
        return 1
    for clu in cluster_list:
        if not clu.keep_contigs or not clu.continue_iterating:
            pass
        else:
            clu.global_start_coord = len(faux_genome) + nbuff
            with open(clu.mappings[-1].assembled_contig, 'r') as con:
                contig_rec = list(SeqIO.parse(con, 'fasta'))[0]
            faux_genome = str(faux_genome + nbuffer + contig_rec.seq)
            clu.global_end_coord = len(faux_genome)
            # lastly, set cluster name to new sequence name
            clu.sequence_id = new_seq_name
            counter = counter + 1
    if counter == 0:
        logger.warning("No viable contigs for faux genome construction!")
        return 1
    else:
        logger.info("combined %s records as genome for next round of mapping",
                    counter)
    record = SeqRecord(Seq(str(faux_genome + nbuffer),
                           IUPAC.IUPACAmbiguousDNA()),
                       id=new_seq_name)

    outpath = os.path.join(output_root,
                           "iter_{0}_buffered_genome.fasta".format(iteration))
    with open(outpath, 'w') as outf:
        SeqIO.write(record, outf, 'fasta')
    return (outpath, len(record))


def decide_proceed_to_target(target_len, logger=None):
    assert logger is not None, "Must use logging!"
    if target_len is not None:
        if not target_len > 0 or not isinstance(target_len, float):
            logger.error("--target_len is set to invalid value! Must be a " +
                         "decimal greater than zero, ie where 1.1 would be " +
                         "110% of the original sequence length.")
            raise ValueError
        elif target_len > 5 and 50 > target_len:
            logger.error("We dont reccommend seeding to lengths greater than" +
                         "5x original seed length. Try between 0.5 and 1.5." +
                         "  If you are setting a target number of bases, it " +
                         " must be greater than 50")
            raise ValueError
        else:
            proceed_to_target = True
    else:
        proceed_to_target = False
    return proceed_to_target


def subprocess_run_list(cmdlist, hard=False, logger=None):
    """ This just allows for sequential cmds with multiprocessing.
    It prevents the errors when future commands are looking for and not finding
    a needed file.
    Logger cant be used with multiprocessing
    returns 0 if all is well, otherwise returns 1
    if hard == True, quits instead of returning 1
    """
    for cmd in cmdlist:
        if cmd is not None:
            try:
                subprocess.run([cmd],
                               shell=sys.platform != "win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
            except Exception as e:
                if logger:
                    logger.error(e)
                if hard:
                    sys.exit(1)
                else:
                    return 1
    return 0


def copy_to_handy_dir(outdir, pre, ref_gb, seedGenome,
                      skip_control=False,
                      hard=False, logger=None):
    """ copy the resulting contigs and reference to dir for mauving
    """
    assert logger is not None, "Must use logging"
    os.makedirs(outdir)
    files_to_copy = [
        os.path.join(seedGenome.output_root,
                     "final_de_fere_novo_assembly", "contigs.fasta"),
        os.path.join(seedGenome.output_root,
                     "final_de_novo_assembly", "contigs.fasta"),
        ref_gb]
    new_names = [pre + "_de_fere_novo_contigs.fasta",
                 pre + "_de_novo_contigs.fasta",
                 pre + ".gb"]
    if skip_control:
        new_names = [new_names[i] for i in [0,2]]
        file_to_copy = [file_to_copy[i] for i in [0,2]]
    for idx, f in enumerate(files_to_copy):
        logger.debug("copying %s to %s as %s",
                     f,
                     os.path.join(outdir, os.path.basename(f)),
                     new_names[idx])
        try:
            shutil.copyfile(
                f,
                os.path.join(outdir, new_names[idx]))
        except:
            logger.warning("unable to copy %s to %s.",
                           f, os.path.join(outdir, os.path.basename(f)))
            logger.error(last_exception())
            if hard:
                raise FileNotFoundError


def printPlot(data, line=None, ymax=30, xmax=60, tick=.2,
              title="test", fill=False, pathIfWrite=None, logger=None):
    """ ascii not what your program can do for you...

    Plot a graph of alignment score by read count, using the logger info
    output.  Optionally write this to a file too.

    Args:
        data (list): list of ints representing alignment scores
        line (int): where to draw a cutoff line
        ymax (int): how tall is our graph
        xmax (int): how wide is our graph
        fill (bool): wherther to fill the area under the line
        pathIfWrite (str): path to file to write graph to, if not None
    Raises:
        None
    Returns:
        None
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
    if pathIfWrite is not None:
        with open(pathIfWrite, 'w') as file_handler:
            for item in plotlines:
                file_handler.write("{0}\n".format(item))


# def plotAsScores(score_list, score_min, outdir, logger=None):
#     """ This function is depreciated til I can figure out how to get
#         conda to play nice with matplotlib.  No one (not even me)
#         cares about this graph anyway #singletear
#     """
#     assert logger is not None, "must use logging"
#     basename = os.path.join(outdir, "AS_score_plot")
#     if len(score_list) > 200000:
#         logger.info("Downsampling our pltting data to 20k points")
#         score_list = random.sample(score_list, 200000)
#     xmax = max(score_list)
#     ymax = max(set(score_list), key=score_list.count)
#     logger.info("score_min, xmax and ymax: %s, %s, %s", score_min, xmax, ymax)
#     # plt.figure()  # <- makes a new figure and sets it active (add this)
#     fig = Figure()
#     FigureCanvas(fig)
#     ax1 = fig.add_subplot(121)
#     ax2 = fig.add_subplot(122)

#     # fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False)
#     # ,
#     #                                  gridspec_kw={'height_ratios': [1, 1]})
#     ax1.hist(score_list, bins=50, color='b', alpha=0.9, label='Binned')
#     ax1.hist(score_list, bins=100,
#               cumulative=True, color='r', alpha=0.5,
#               label='Cumulative')
#     ax1.set_title('Read Alignment Score Histogram')
#     ax1.set_xlabel('Alignemnt Score')
#     ax1.set_ylabel('Abundance')
#     ax1.plot([score_min, score_min], [0, len(score_list) * 1.1],
#               color='green', linewidth=5, alpha=0.6)
#     ax1.axis([0, max(score_list), 0, len(score_list) * 1.1])
#     # ax1.set_yscale("log", nonposy='clip')
#     ax1.legend()
#     fig.tight_layout()
#     ax1.grid(True)
#     ax2.grid(True)

#     ax2.scatter(y=sorted(score_list, reverse=True),
#                  x=range(0, len(score_list)))
#     ax2.plot([0, len(score_list) * 1.1], [score_min, score_min],
#               color='green', linewidth=5, alpha=0.6)
#     ax2.axis([0, len(score_list) * 1.1, 0, max(score_list) * 1.1])
#     ax2.set_title('Read Alignment Score, Sorted')
#     ax2.set_ylabel('Alignment Score')
#     ax2.set_xlabel('Index of Sorted Read')
#     fig.set_size_inches(12, 7.5)
#     fig.savefig(str(basename + '.png'), dpi=(200))
#     fig.savefig(str(basename + '.pdf'), dpi=(200))
#     logger.info("Plotting alignment score of mapping:")
#     logger.info("Filled area represents the reads retained after filtering. " +
#                 "If it looks like " +
#                 "this filtering threshold is inapporpriate, consider " +
#                 "adjusting with --score_min")


def set_ref_as_contig(ref_arg, map_percentage, final=False, logger=None):
    """ sets the ref_as_contig arg for spades.
    # DEPRECIATED COMMENT
    # Note this is used to set initial subassembly and final spades assembly.
    # Intermediate runs will always use seeds as "trusted", becasue we do not
    # reuse reads
    this allows the user to submit one of 4 options:
    ignore: do not use regions from reference in initial subassembly
    """
    assert logger is not None, "must use logging"
    if ref_arg == "infer":
        if map_percentage > 80:
            ref_as_contig = "trusted"
        else:
            ref_as_contig = "untrusted"
            logger.info(
                str("unfiltered mapping percentage is %f2 so " +
                    "'ref_as_contigs' is set to %s"),
                map_percentage, ref_as_contig)
    elif ref_arg == "ignore":
        ref_as_contig = None
    else:
        assert ref_arg in [None, "trusted", "untrusted"], \
            "error parsing initial ref_as_contig: %s" % ref_arg
        ref_as_contig = ref_arg
    if final and ref_as_contig is None:
        # if final assembly, --ref-as-contig cannot be none, so we infer
        return set_ref_as_contig(ref_arg="infer",
                                 map_percentage=map_percentage,
                                 final=final, logger=logger)
    return ref_as_contig


def report_region_depths(inp, logger):
    assert logger is not None, "must use logging"
    report_list = []
    nclusters = len(inp[0])
    for i in range(0, nclusters):
        report_list.append("Cluster %i:" % i)
        for itidx, iteration in enumerate(inp):
            for cluster in iteration:
                if cluster[0] == i:
                    report_list.append(
                        "\tIter %i -- 5' coverage: %.2f  3' coverage %.2f" % (
                            itidx, cluster[1], cluster[2]))
    return report_list


def make_modest_spades_cmd(cmd, cores, memory, split=0,
                           serialize=False, logger=None):
    """ adjust spades commands to use set amounts of cores and memory
    returns the command, split on "--careful", cause why would you run
    SPAdes with "--reckless"?
    if you need to split resouces between, say, two commands, use split 2
    """
    assert logger is not None, "Must use logging"
    cmdA, cmdB = cmd.split("--careful")
    if serialize:
        logger.debug("Allocating SPAdes %dgb of memory", memory)
        mem_each = memory
        cores_each = cores
    else:
        if split:
            if split > cores or split > memory:
                logger.error("cannot split cores or memory resources into " +
                             "values less than 1!")
                sys.exit(1)
            mem_each = int(memory / split)
            if mem_each < 1:
                logger.warning(
                    "you must have at least 1gb memory allocated " +
                    "to each spades call! allocating minimum")
                mem_each = 1
            cores_each = int(cores / split)
            logger.info("Running spades with %d cores and %d gb memory",
                        cores_each, mem_each)
        else:
            # make sure spades doesnt hog processors or ram
            mem_each = int(memory / cores)  # should be floor
            if mem_each < 1:
                logger.warning("you must have at least 1gb memory allocated " +
                             "to each spades call!  allocating minimum")
                mem_each = 1
            cores_each = 1
            logger.info(
                "Allocating SPAdes %dgb of memory for each of %d cores",
                mem_each, cores)
    return "{0}-t {1} -m {2} --careful{3}".format(
        cmdA,
        cores_each,
        mem_each,
        cmdB)


def define_score_minimum(args, readlen, iteration, logger):
    """ given the args, define your scoring minimum for mapping
    for smalt, each round of mapping gets more stringent.  For BWA,
    the mapping minimum is deferred to BWA's default settings.

    This will probably be depreciated soon when we finally get rid of the
    smalt options
    """
    assert logger is not None, "must use logging"
    if not args.score_min:
        # This makes it such that score minimum is now more stringent
        # with each mapping.  Too campy? probably.
        if args.mapper == 'smalt':
            scaling_factor = 1.0 - (
                1.0 / (2.0 + float(iteration)))
            score_minimum = int(readlen * scaling_factor)
            logger.info(
                "Mapping with min_score of %f2 (%f2 of read length, %f2)",
                scaling_factor, score_minimum, readlen)
        else:
            assert args.mapper == 'bwa', "must be either smalt or bwa"
            logger.debug("using the default minimum score for BWA")
            score_minimum = None
    else:
        if readlen < args.score_min:
            raise ValueError(
                "--min_score must be smaller than read length {0}".format(
                    readlen))
        score_minimum = args.score_min
        logger.info(
            "Mapping with min_score of %f2 (read length: %f2)",
            score_minimum, readlen)
    return score_minimum


def check_genbank_for_fasta(gb, logger=None):
    """ Ensure user has not provided a fasta by accident
    Yes, its a dumb function, but it was neccsary, given the amount of times
    I would get partway through debugging some error only to realize I
    used a genbank that was actually a fasta
    """
    assert logger is not None, "must use logging"
    with open(gb) as ingb:
        for idx, line in enumerate(ingb):
            if line.startswith(">"):
                logger.error("This genbank file looks like a fasta! Exiting")
                raise ValueError


def get_fasta_consensus_from_BAM(samtools_exe, bcftools_exe, # vcfutils_exe,
                                 outfasta, ref, bam, region=None, logger=None):
    """ run system commands to get consensus fastq from fastq
    """
    cmd_list, consensus_fa = make_get_consensus_cmds(
        samtools_exe=samtools_exe,
        bcftools_exe=bcftools_exe,
        # vcfutils_exe=vcfutils_exe,
        ref=ref, bam=bam,
        region=region,
        outfasta=outfasta,
        old_method=True,
        logger=logger)
    for cmd in cmd_list:  # may have more cmds here in future
        logger.debug(cmd)
        subprocess.run([cmd],
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    # consensus_fa = convert_fastq_to_fasta(
    #     fastq=consensus_fq, outfasta=outfasta,
    #     only_ATCG=True, logger=logger)
    return consensus_fa


def convert_fastq_to_fasta(fastq, outfasta,  # only_first=True,
                           only_ATCG=True, logger=None):
    """ This converts fastq to fasta, as well as converting non ATCG bases
    to N's.  Currently only works on single entry fastqs.
    If you need to convert more, learn how to use awk.
    """
    assert logger is not None, "must use logging"
    from Bio.Alphabet import IUPAC
    counter = 0
    with open(outfasta, "w") as ofile:
        with open(fastq, "r") as ifile:
            for read in SeqIO.parse(ifile, "fastq"):
                if counter > 0:  # and only_first:
                    raise ValueError
                if only_ATCG:
                    new_record = read
                    new_seq = ""
                    for i in range(0, len(read.seq)):
                        if read.seq[i] not in ["A","T", "C", "G"]:
                            new_seq = new_seq + "N"
                        else:
                            new_seq = new_seq + read.seq[i]
                    new_record.seq = Seq(new_seq, IUPAC.ambiguous_dna)
                    SeqIO.write(new_record, ofile, "fasta")
                else:
                    SeqIO.write(read, ofile, "fasta")
                counter = counter + 1
    return outfasta


def make_get_consensus_cmds(samtools_exe, bcftools_exe,# vcfutils_exe,
                            ref, bam, outfasta=None, old_method=False,
                            region=None, logger=None):
    """ use samtoolsa nd vcfutils to get a consesnus fastq from a BAM file
    the old_method referes to how bcftools calls consensuses since 1.2
    """
    sorted_bam = os.path.splitext(bam)[0] + "_sorted.bam"
    if outfasta is None:
        consensus_fa = os.path.splitext(bam)[0] + "_consensus.fa"
    else:
        consensus_fa = outfasta
    faidx_cmd = "{0} faidx {1}".format(samtools_exe, ref)
    if old_method:
        call = "-c"
    else:
        call = "-m"
    if region is None:
        this_region = ""
    else:
        assert isinstance(region, str), \
            "region must be a string in the form of 'chr:start-end'"
        this_region = "-r {0} ".format(region)
    temp_zvcf = os.path.join(os.path.dirname(sorted_bam), "consensus.vcf.gz")

    sort_cmd = "{0} sort {1} > {2}".format(samtools_exe, bam, sorted_bam)
    index_cmd = "{0} index {1}".format(samtools_exe, sorted_bam)
    # -A means ignore strainge pairing (which we need, as we are subsetting a region already)
    # -E means recalculate scores, for imporved accuracy
    # -u for uncocmpressed output
    # -a for all positions, even those with low coverage
    vcf_cmd = str("{0} mpileup -d8000 -EA -uf {1}  {5}{2} | " +
                  "{3} call --ploidy 1 {4} -Oz -o {6} -").format(
                      samtools_exe,  #0
                      ref,  #1
                      sorted_bam,  #2
                      bcftools_exe,  #3
                      call,  #4
                      this_region,  #5
                      temp_zvcf) #6
    tabix_cmd = "{0} {1}".format("tabix", temp_zvcf)
    consensus_cmd = "{0} faidx {1} {2} | {3} consensus {4} > {5}".format(
        samtools_exe,  #0
        ref,  #1
          #2, note this is not "this_region", as we do not need the -r
        region if region is not None else "",
        bcftools_exe,  #3
        temp_zvcf,  #4
        consensus_fa)
    return ([faidx_cmd, sort_cmd, index_cmd, vcf_cmd,
             tabix_cmd, consensus_cmd], consensus_fa)


def main(args, logger=None):
    # allow user to give relative paths
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory already exists; exiting...")
        sys.exit(1)
    t0 = time.time()
    log_path = os.path.join(output_root, "riboSeed.log")
    if logger is None:
        logger = set_up_logging(verbosity=args.verbosity,
                                outfile=log_path,
                                name=__name__)
    # # log version of riboSeed, commandline options, and all settings
    logger.info("riboSeed pipeline package version: %s",
                __version__)

    logger.info("Usage:\n%s\n", " ".join([x for x in sys.argv]))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, str(v))
    logger.debug("current PATH:")
    try:
        logger.debug(os.environ['PATH'])
    except KeyError:
        logger.error("no PATH variable found in system environment.")
        sys.exit(1)
    if args.cores is None:
        args.cores = multiprocessing.cpu_count()
    logger.info("Using %i core(s)", args.cores)
    logger.info("Using %iGB memory", args.memory)
    logger.info("Using %i thread(s)", args.threads)

    logger.info("checking for installations of all required external tools")
    logger.debug("creating an Exes object")
    try:
        sys_exes = Exes(python=sys.executable,
                        samtools=args.samtools_exe,
                        spades=args.spades_exe,
                        bwa=args.bwa_exe,
                        smalt=args.smalt_exe,
                        quast=args.quast_exe,
                        bcftools=args.bcftools_exe,
                        method=args.mapper)
        sys_exes.python = sys_exes.check_spades_python_version(logger=logger)
    except Exception as e:
        logger.error(e)
        sys.exit(1)

    logger.debug("All required system executables found!")
    logger.debug(str(sys_exes.__dict__))
    try:
        samtools_version = check_version_from_cmd(
            exe=sys_exes.samtools, cmd='', line=3, where='stderr',
            pattern=r"\s*Version: (?P<version>[^(]+)",
            min_version=SAMTOOLS_MIN_VERSION, logger=logger)
    except Exception as e:
        logger.info(e)
        logger.info("perhaps conda is giving a build warning?")
        try:
            samtools_version = check_version_from_cmd(
                exe=sys_exes.samtools, cmd='', line=4, where='stderr',
                pattern=r"\s*Version: (?P<version>[^(]+)",
                min_version=SAMTOOLS_MIN_VERSION, logger=logger)
        except Exception as f:
            logger.error(f)
            sys.exit(1)
    if args.additional_libs:
        logger.debug("Checking additional args for final assembly ahead of time")
        check_spades_extra_library_input(args.additional_libs)
    logger.debug("samtools version: %s", samtools_version)
    if samtools_version.startswith("1.5"):
        logger.error("Cannot use samtools 1.5! see samtools github issue #726")
        sys.exit(1)
    # check bambamc is installed proper if using smalt
    if args.mapper == "smalt":
        logger.info("SMALT is the selected mapper")
        test_smalt_cmds = get_smalt_full_install_cmds(smalt_exe=sys_exes.smalt,
                                                      logger=logger)
        test_smalt_bam_install(cmds=test_smalt_cmds, logger=logger)
    else:
        logger.info("BWA is the selected mapper")
    # if --Initial_consensus, ensure our optional programs are in working order
    if args.initial_consensus:
        if any([x is None for x in [sys_exes.bcftools]]):
            logger.error("Must have availible executables for both bcftools " +
                         "and vcfutils if using `--inital_consensus option! " +
                         "Exiting (1)")
            sys.exit(1)

    try:
        check_genbank_for_fasta(gb=args.reference_genbank, logger=logger)
    except ValueError:
        sys.exit(1)
    # if the target_len is set. set needed params
    try:
        proceed_to_target = decide_proceed_to_target(
            target_len=args.target_len, logger=logger)
    except ValueError:
        logger.error("Exiting")
        sys.exit(1)
    # check and warn user about potential RAM issues
    if args.memory < 6 or int(args.memory / args.cores) < 6:
        logger.warning("Danger!  We recommend that you have a minimum of " +
                       "6GB memory (or 6GB memory per core is using " +
                       "multiprocessing) available. If you have less " +
                       "than that per core, SPAdes may run out of memory.")
        logger.warning("You can continue as configured if needed, and if a " +
                       "SPAdes error occurs, you can still use the long " +
                       "reads generated by riboSeed in a standalone assembly")
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

    # make seedGenome object
    logger.debug("constructing the seedGenome object")
    seedGenome = SeedGenome(
        name=os.path.basename(os.path.splitext(args.reference_genbank)[0]),
        # this needs to be zero indexed to access mappings by iter
        this_iteration=0,
        iter_mapping_list=[],
        max_iterations=args.iterations,
        clustered_loci_txt=args.clustered_loci_txt,
        output_root=output_root,
        unmapped_mapping_list=[],
        genbank_path=args.reference_genbank,
        logger=logger)

    seedGenome.iter_mapping_list[0].ref_fasta = seedGenome.ref_fasta
    # add ngslib object for user supplied NGS data
    logger.debug("adding the sequencing libraries to the seedGenome")
    logger.debug(args.fastqS1)
    seedGenome.master_ngs_ob = NgsLib(
        name="master",
        master=True,
        make_dist=args.mapper == "smalt",
        readF=args.fastq1,
        readR=args.fastq2,
        readS0=args.fastqS1,
        logger=logger,
        mapper_exe=sys_exes.mapper,
        ref_fasta=seedGenome.ref_fasta)
    if args.force_kmers:
        logger.info("Skipping kmer check, using kmers provided")
        checked_k = args.kmers
        checked_prek = args.pre_kmers
    else:
        checked_k = check_kmer_vs_reads(
            k=args.kmers,
            readlen=seedGenome.master_ngs_ob.readlen,
            min_diff=2, logger=logger)
        checked_prek = check_kmer_vs_reads(
            k=args.pre_kmers,
            readlen=seedGenome.master_ngs_ob.readlen,
            min_diff=2, logger=logger)
    logger.debug("Using the following kmer values pre and final assemblies:")
    logger.debug(checked_k)
    logger.debug(checked_prek)

    if "pe" in seedGenome.master_ngs_ob.libtype:
        # check equal length fastq.  This doesnt actually check propper pairs
        logger.debug("Checking that the fastq pair have equal number of reads")
        try:
            check_fastqs_len_equal(file1=args.fastq1, file2=args.fastq2)
        except Exception as e:
            # not just value error, whatever file_len throws
            logger.error(last_exception())
            sys.exit(1)
    # read in riboSelect clusters, make a lociCluster ob for each,
    # which get placed in seedGenome.loci_clusters
    try:
        seedGenome.loci_clusters = parse_clustered_loci_file(
            filepath=seedGenome.clustered_loci_txt,
            gb_filepath=seedGenome.genbank_path,
            output_root=output_root,
            circular=args.linear is False,
            logger=logger)
        seedGenome.loci_clusters = add_gb_seqrecords_to_cluster_list(
            cluster_list=seedGenome.loci_clusters,
            gb_filepath=seedGenome.genbank_path)
    except Exception as e:
        logger.error(e)
        logger.error(last_exception())
        sys.exit(1)

    # add coordinates to lociCluster.loci_list
    try:
        add_coords_to_clusters(seedGenome=seedGenome,
                               logger=logger)
    except Exception as e:
        logger.error(e)
        logger.error(last_exception())
        sys.exit(1)
    try:
        logger.info("padding genbank by %i", args.flanking * 3)
        logger.debug("old ref_fasta: %s", seedGenome.ref_fasta)
        seedGenome.pad_genbank(pad=args.flanking * 3,
                               circular=args.linear is False, logger=logger)
        logger.debug("new ref_fasta: %s", seedGenome.ref_fasta)
    except Exception as e:
        logger.error(e)
        logger.error(last_exception())
        sys.exit(1)

    # make first iteration look like future iterations:
    # this should also ensure the mapper uses the padded version
    seedGenome.next_reference_path = seedGenome.ref_fasta
    #
    for cluster in seedGenome.loci_clusters:
        cluster.master_ngs_ob = seedGenome.master_ngs_ob
# ---------------------------------------------------------------------------
    # Performance summary lists
    mapping_percentages = []
    region_depths = []
    # this gets set during the first mapping
    score_minimum = 0
    # now, we need to assemble each mapping object
    # this should exclude any failures
    while seedGenome.this_iteration < args.iterations:
        logger.info("processing iteration %i", seedGenome.this_iteration)
        logger.debug("with new reference: %s", seedGenome.next_reference_path)
        clusters_to_process = [x for x in seedGenome.loci_clusters if
                               x.continue_iterating and
                               x.keep_contigs]
        if len(clusters_to_process) == 0:
            logger.error("No clusters had sufficient mapping! Exiting")
            sys.exit(1)
        if len(clusters_to_process) < len(seedGenome.loci_clusters):
            logger.warning(
                "clusters excluded from this iteration \n%s",
                " ".join([str(x.index) for x in
                          seedGenome.loci_clusters if
                          x.index not in [y.index for
                                          y in clusters_to_process]]))
        # For each (non-inital) iteration
        if seedGenome.this_iteration != 0:
            if seedGenome.this_iteration != 1:
                # clear out old .sam files to save space
                if args.clean_temps:
                    logger.info("removing uneeded files from previous mappings")
                    # delete the read files from the last mapping
                    # dont do this on first iteration cause those be the reads!
                    # and if they aren't backed up you are up a creek and
                    # probably very upset with me.
                    seedGenome.purge_old_files(all_iters=False, logger=logger)

            # seqrecords for the clusters to be gen.next_reference_path
            with open(seedGenome.next_reference_path, 'r') as nextref:
                next_seqrec = list(SeqIO.parse(nextref, 'fasta'))[0]  # next?
            for clu in clusters_to_process:
                clu.seq_record = next_seqrec
            # print qualities of mapped reads
            for clu in clusters_to_process:
                logger.debug("getting mapping scores for cluster %i from %s",
                             clu.index, clu.mappings[-1].mapped_bam)
                mapped_scores = get_bam_AS(
                    inbam=clu.mappings[-1].mapped_bam,
                    logger=logger)
                if len(mapped_scores) > 200000:
                    logger.debug("Downsampling our plotting data to 20k points")
                    mapped_scores = random.sample(mapped_scores, 200000)
                printPlot(data=mapped_scores, line=score_minimum,
                          ymax=18, xmax=60, tick=.2, fill=True,
                          pathIfWrite=None,
                          title=str("Average alignment Scores for cluster " +
                                    "%i\n " % clu.index),
                          logger=logger)
                logger.info(str("-" * 72))

            # make new ngslib from unampped reads
            convert_cmd, unmapped_ngsLib = convert_bam_to_fastqs_cmd(
                mapping_ob=seedGenome.iter_mapping_list[
                    seedGenome.this_iteration - 1],
                samtools_exe=sys_exes.samtools, single=True,
                # ref fasta is used to make index cmd
                ref_fasta=seedGenome.next_reference_path,
                which='unmapped', logger=logger)
            # unless subtract arg is used, use all reads each mapping
            if not args.subtract:
                unmapped_ngsLib = seedGenome.master_ngs_ob

            unmapped_ngsLib.readlen = seedGenome.master_ngs_ob.readlen
            unmapped_ngsLib.smalt_dist_path = \
                seedGenome.master_ngs_ob.smalt_dist_path
            logger.debug("converting unmapped bam into reads:")
            seedGenome.master_ngs_ob.ref_fasta = seedGenome.next_reference_path
            # dont worry about wasting time making these libraries if
            # not subtracting previously mapped reads
            if args.subtract:
                for cmd in [convert_cmd]:  # may have more cmds here in future
                    logger.debug(cmd)
                    subprocess.run([cmd],
                                   shell=sys.platform != "win32",
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   check=True)
        else:
            # start with whole lib if first time through
            unmapped_ngsLib = seedGenome.master_ngs_ob

        score_minimum = define_score_minimum(
            args=args, iteration=seedGenome.this_iteration,
            readlen=unmapped_ngsLib.readlen, logger=logger)

        try:
            nonify_empty_lib_files(unmapped_ngsLib, logger=logger)
        except ValueError:
            logger.error(
                "No reads mapped for this iteration. This could be to an " +
                "error from samtools or elevated mapping stringency.")
            if seedGenome.this_iteration != 0:
                logger.warning(" proceeding to final assemblies")
                break
            else:
                logger.error(" Exiting!")
                sys.exit(1)
        # Run commands to map to the genome
        # the exe argument is Exes.mapper because that is what is checked
        # during object instantiation
        if args.mapper == "smalt":
            # # get rid of bwa mapper default args
            # if args.mapper_args == '-L 0,0 -U 0':
            #     args.mapper_args =
            map_percent = map_to_genome_ref_smalt(
                mapping_ob=seedGenome.iter_mapping_list[
                    seedGenome.this_iteration],
                ngsLib=unmapped_ngsLib,
                cores=(args.cores * args.threads),
                samtools_exe=sys_exes.samtools,
                genome_fasta=seedGenome.next_reference_path,
                smalt_exe=sys_exes.mapper,
                score_minimum=score_minimum,
                step=3, k=5,
                scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
                logger=logger)
        else:
            assert args.mapper == "bwa", "must be either bwa or smalt"
            map_percent, score_list, score_minimum = map_to_genome_ref_bwa(
                mapping_ob=seedGenome.iter_mapping_list[
                    seedGenome.this_iteration],
                ngsLib=unmapped_ngsLib,
                cores=(args.cores * args.threads),
                genome_fasta=seedGenome.next_reference_path,
                samtools_exe=sys_exes.samtools,
                bwa_exe=sys_exes.mapper,
                score_minimum=score_minimum,
                # add_args='-L 0,0 -U 0',
                add_args=args.mapper_args,
                logger=logger)
        mapping_percentages.append([seedGenome.this_iteration, map_percent])
        ##
        # Here we put the logic in for a wee little preliminary assembly to test
        ##
        # if things go really bad on the first mapping, get out while you can
        if len(score_list) == 0:
            logger.error(
                "No reads mapped for this iteration. This could be to an " +
                "error from samtools, bwa mem, a bad reference, " +
                " or elevated mapping stringency. ")
            if seedGenome.this_iteration != 0:
                logger.warning(" proceeding to final assemblies")
                break
            else:
                logger.error("Exiting!")
                sys.exit(1)

        # on first time through, infer ref_as_contig if not
        #   provided via commandline
        if seedGenome.this_iteration == 0:
            # do info for smalt mapping
            if args.mapper == "bwa":
                fig_dir = os.path.join(output_root, "figs")
                os.makedirs(fig_dir)
                # either use defined min or use the same heuristic as mapping
                # if PLOT:
                #     plotAsScores(
                #         score_list,
                #         score_min=score_minimum if score_minimum is not None else
                #         int(round(float(seedGenome.master_ngs_ob.readlen) / 2.0)),
                #         outdir=fig_dir, logger=logger)
                printPlot(data=score_list, line=score_minimum,
                          ymax=30, xmax=60,
                          pathIfWrite=os.path.join(
                              fig_dir, "initial_mapping.txt"),
                          tick=.2, fill=True,
                          title="Average alignment Scores (y) by sorted " +
                          "read index (x)",
                          logger=logger)
                subassembly_ref_as_contig = set_ref_as_contig(
                    ref_arg=args.ref_as_contig,
                    map_percentage=map_percent, logger=logger)
            else:
                subassembly_ref_as_contig = args.ref_as_contig
        else:
            pass

        try:
            # again, this is [(idx, start_depth, end_depth)]
            # clusters_post_partition are the clusters passing the minimum
            # depth on the flanking regions.
            iter_depths, clusters_to_subassemble = partition_mapping(
                seedGenome=seedGenome,
                logger=logger,
                samtools_exe=sys_exes.samtools,
                flank=args.flanking,
                min_flank_depth=args.min_flank_depth,
                cluster_list=clusters_to_process)

        except Exception as e:
            logger.error("Error while partitioning reads from iteration %i",
                         seedGenome.this_iteration)
            logger.error(last_exception())
            logger.error(e)
            sys.exit(1)

        logger.debug(iter_depths)
        region_depths.append(iter_depths)
        extract_convert_assemble_cmds = []
        # generate spades cmds (cannot be multiprocessed becuase of python's
        #  inability to pass objects to multiprocessing)
        # subassembly_ref_as_contig must be 'trusted' here because of the multimapping/
        #  coverage issues
        for cluster in clusters_to_subassemble:
            if args.initial_consensus and seedGenome.this_iteration == 0:
                consensus_fasta = get_fasta_consensus_from_BAM(
                    samtools_exe=sys_exes.samtools,
                    bcftools_exe=sys_exes.bcftools,
                    region="{0}:{1}-{2}".format(
                        cluster.sequence_id,
                        # cluster.global_start_coord - args.flanking,
                        # cluster.global_end_coord + args.flanking),
                        cluster.global_start_coord,
                        cluster.global_end_coord),
                    # region=None,
                    outfasta=os.path.join(
                        cluster.mappings[-1].mapping_subdir,
                        "consensus.fasta"),
                    ref=seedGenome.next_reference_path,
                    bam=seedGenome.iter_mapping_list[0].sorted_mapped_bam,
                    # bam=cluster.mappings[-1].mapped_bam,
                    logger=logger)
                # assign our consensus seqeunce to replace the segment of
                # reference fasta used in initial subassembly
                cluster.mappings[-1].ref_fasta = consensus_fasta
            cmdlist = []
            logger.debug("generating commands to convert bam to fastqs " +
                         "and assemble long reads")
            convert_cmds, new_ngslib = convert_bam_to_fastqs_cmd(
                mapping_ob=cluster.mappings[-1], which='mapped',
                single=True,
                samtools_exe=sys_exes.samtools,
                ref_fasta=cluster.mappings[-1].ref_fasta, logger=logger)
            cmdlist.append(convert_cmds)
            spades_cmd = generate_spades_cmd(
                mapping_ob=cluster.mappings[-1],
                ngs_ob=new_ngslib, single_lib=True,
                ref_as_contig=subassembly_ref_as_contig,
                # only check if files exist if using a subset of reads
                check_libs=args.subtract,
                python_exe=sys_exes.python,
                as_paired=False,
                prelim=True,
                k=checked_prek,
                spades_exe=sys_exes.spades, logger=logger)
            # setting some thread limits here
            modest_spades_cmd = make_modest_spades_cmd(
                cmd=spades_cmd, cores=args.cores, memory=args.memory,
                serialize=args.serialize, logger=logger)
            cmdlist.append(modest_spades_cmd)

            cluster.mappings[-1].mapped_ngslib = new_ngslib
            extract_convert_assemble_cmds.append(cmdlist)

        # run all those commands!
        logger.debug(
            "\n running %i cmds: \n %s",
            len([j for i in extract_convert_assemble_cmds for j in i]),
            "\n".join([j for i in extract_convert_assemble_cmds for j in i]))
        if args.serialize:
            subassembly_return_sum = 0
            logger.info("running without multiprocessing!")
            for cmd in [j for i in extract_convert_assemble_cmds for j in i]:
                logger.debug(cmd)
                result = subprocess.run([cmd],
                                        shell=sys.platform != "win32",
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        check=False)
                logger.debug(result.returncode)
                subassembly_return_sum = \
                    subassembly_return_sum + result.returncode

        else:
            pool = multiprocessing.Pool(processes=args.cores)
            results = [
                pool.apply_async(subprocess_run_list,
                                 (cmds,),
                                 {"logger": None,
                                  "hard": False})
                for cmds in extract_convert_assemble_cmds]
            pool.close()
            pool.join()
            subassembly_return_sum = sum([r.get() for r in results])
        # check return codes
        logger.info("Sum of return codes (should be 0):")
        logger.info(subassembly_return_sum)
        if subassembly_return_sum != 0:
            logger.error(
                "%d error(s) occurred when converting reads and subassembling with SPAdes!",
                subassembly_return_sum)
            if args.damn_the_torpedos or subassembly_return_sum < 2:
                logger.error(
                    "Check the SPAdes and samtools logs to diagnose, especially if " +
                    "this occurs with more than one subassembly. Continuing")
            else:
                sys.exit(1)


        # evaluate mapping (cant be multiprocessed)
        for cluster in clusters_to_process:
            cluster.assembly_success = evaluate_spades_success(
                clu=cluster,
                mapping_ob=cluster.mappings[-1],
                include_short_contigs=args.include_short_contigs,
                keep_best_contig=True,
                min_delta=10,
                flank=args.flanking,
                seqname='', logger=logger,
                min_assembly_len=args.min_assembly_len,
                proceed_to_target=proceed_to_target,
                target_len=args.target_len)
            parse_subassembly_return_code(
                cluster=cluster,
                logger=logger)
        clusters_for_pseudogenome = [
            x for x in seedGenome.loci_clusters if
            x.continue_iterating and x.keep_contigs]
        if len(clusters_for_pseudogenome) != 0:
            faux_genome_path, faux_genome_len = make_faux_genome(
                seedGenome=seedGenome,
                iteration=seedGenome.this_iteration,
                output_root=seedGenome.output_root,
                nbuff=1000,
                cluster_list=[x for x in clusters_for_pseudogenome if
                              x.continue_iterating],
                logger=logger)
            logger.info("Length of buffered 'genome' for mapping: %i",
                        faux_genome_len)
        else:
            faux_genome_path = 1
            seedGenome.this_iteration = args.iterations + 1
        seedGenome.this_iteration = seedGenome.this_iteration + 1
        seedGenome.next_reference_path = faux_genome_path
        if seedGenome.this_iteration >= args.iterations:
            logger.info("moving on to final assemblies!")
        else:
            logger.info("Moving on to iteration: %i",
                        seedGenome.this_iteration)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
    # done with the iterations!  Lets free up some space
    if args.clean_temps:
        seedGenome.purge_old_files(all_iters=True, logger=logger)
    # And add the remaining final contigs to the directory for combination
    contigs_moved_before_list_iter = \
        [x.final_contig_path for x in seedGenome.loci_clusters
         if not x.keep_contigs]
    logger.debug("Contigs moved prior to final iteration:")
    logger.debug(contigs_moved_before_list_iter)
    if len(contigs_moved_before_list_iter) == len(seedGenome.loci_clusters):
        logger.info("all contigs already copied to long_reads dir")
    else:
        # this assumes that all the uable clusters not in
        # clusters_for_pseudogenome have already been copied over
        for clu in [x for x in clusters_for_pseudogenome if x.keep_contigs]:
            try:
                shutil.copyfile(clu.mappings[-1].assembled_contig,
                                os.path.join(
                                    seedGenome.final_long_reads_dir,
                                    "{0}_cluster_{1}_iter_{2}.fasta".format(
                                        clu.sequence_id,
                                        clu.index,
                                        clu.mappings[-1].iteration)))

            except Exception as e:
                logger.error(last_exception())
                sys.exit(1)
    for dirpath, dirnames, files in os.walk(seedGenome.final_long_reads_dir):
        if not files:
            logger.error(
                "No pseudocontigs found in the long reads output " +
                "directory; it appears " +
                "that the subassemblies did not yield pseudocontigs " +
                "of sufficient quality.  Exiting with code 0")
            sys.exit(0)
    logger.info("combining contigs from %s", seedGenome.final_long_reads_dir)
    seedGenome.assembled_seeds = combine_contigs(
        contigs_dir=seedGenome.final_long_reads_dir,
        contigs_name="riboSeedContigs",
        logger=logger)
    logger.info("Combined Seed Contigs: %s", seedGenome.assembled_seeds)
    logger.info("Time taken to run seeding: %.2fm" % ((time.time() - t0) / 60))
    report = report_region_depths(inp=region_depths, logger=logger)
    logger.info("Average depths of mapping for each cluster, by iteration:")
    logger.info("\n" + "\n".join(report))

    logger.info("Average mapping percentage for initial mapping: %d",
                mapping_percentages[0][1])
    logger.info("Average mapping percentage of reads to pseudogenome in " +
                "subsequent mappings:\n" +
                "\n".join(["iter %d: %d%%" % (i[0], i[1]) for \
                           i in mapping_percentages[1:]]))

    if args.just_seed:
        logger.info("Done with riboSeed: %s", time.asctime())
        logger.info("Skipping final assembly")
        logger.info("Combined Contig Seeds (for validation or alternate " +
                    "assembly): %s", seedGenome.assembled_seeds)
        logger.info("Time taken: %.2fm" % ((time.time() - t0) / 60))

    final_ref_as_contig = set_ref_as_contig(
        ref_arg=subassembly_ref_as_contig,
        # initial mapping percentage is used if inferring
        map_percentage=mapping_percentages[0][1], final=True, logger=logger)

    spades_quast_cmds, quast_reports = get_final_assemblies_cmds(
        seedGenome=seedGenome, exes=sys_exes,
        cores=args.cores,
        memory=args.memory,
        additional_libs=args.additional_libs,
        serialize=args.serialize,
        ref_as_contig=final_ref_as_contig,
        skip_control=args.skip_control, kmers=checked_k, logger=logger)

    # run final assembly(-ies)
    if args.serialize:
        logger.info("running without multiprocessing!")
        # unpack nested spades quast list, ignoring Nones
        # hate me yet?
        spades_results = []
        for cmd in [j for i in spades_quast_cmds for j in i if j is not None]:
            logger.debug(cmd)
            spades_result = subprocess.run(
                [cmd],
                shell=sys.platform != "win32",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True)
            spades_results.append(spades_result)
        spades_results_sum = sum([r.returncode for r in spades_results])

    else:
        # split the processors based on how many spades_cmds are on the list
        # dont correct for threads, as Spades defaults to lots of threads
        split_cores = int(args.cores / (len(spades_quast_cmds) / 2))
        if split_cores < 1:
            split_cores = 1
        pool = multiprocessing.Pool(processes=split_cores)
        logger.debug("running the following commands:")
        logger.debug("\n".join([j for i in spades_quast_cmds for
                                j in i if j is not None]))
        spades_results = [
            pool.apply_async(subprocess_run_list,
                             (cmds,),
                             {"logger": None,
                              "hard": False})
            for cmds in spades_quast_cmds]
        pool.close()
        pool.join()
        logger.debug(spades_results)
        logger.info("Sum of return codes (should be 0):")
        spades_results_sum = sum([r.get() for r in spades_results])
        if spades_results_sum == 0:
            logger.info(spades_results_sum)
        else:
            logger.warning(spades_results_sum)

    if spades_results_sum != 0:
        logger.error("%d error(s) occurred when running SPAdes!",
                     spades_results_sum)
        if not args.damn_the_torpedos:
            logger.error("Check the spades logs to diagnose.  Exiting (1)")
            sys.exit(1)
    if not args.skip_control and quast_reports is not None:
        logger.debug("writing combined quast reports")
        logger.info("Comparing de novo and de fere novo assemblies:")
        try:
            quast_comp = make_quick_quast_table(
                quast_reports,
                write=True,
                writedir=seedGenome.output_root,
                logger=logger)
            for k, v in sorted(quast_comp.items()):
                logger.info("%s: %s", k, "  ".join(v))
        except Exception as e:
            logger.error("Error writing out combined quast report")
            logger.error(e)
    # make dir for easy downloading from cluster
    copy_to_handy_dir(outdir=os.path.join(output_root, "mauve"),
                      pre=args.experiment_name,
                      ref_gb=args.reference_genbank,
                      skip_control=args.skip_control,
                      seedGenome=seedGenome,
                      hard=False, logger=logger)
    logger.info("Done with riboSeed: %s", time.asctime())
    logger.info("riboSeed Assembly: %s", seedGenome.output_root)
    logger.info("Combined Contig Seeds (for validation or alternate " +
                "assembly): %s", seedGenome.assembled_seeds)
    logger.info("Time taken: %.2fm" % ((time.time() - t0) / 60))

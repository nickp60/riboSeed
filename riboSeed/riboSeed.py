#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
Minor Version Revisions:
 - no more individual versions; all treated as pipeline from here on
starting at version 0.0.940
Created on Sun Jul 24 19:33:37 2016

See README.md for more info and usage

"""
import argparse
import sys
import time
import re
import logging
import os
import shutil
import multiprocessing
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from pyutilsnrw.utils3_5 import set_up_logging, make_outdir, \
    combine_contigs, run_quast, \
    copy_file, check_installed_tools, get_ave_read_len_from_fastq, \
    get_number_mapped, extract_mapped_and_mappedmates, clean_temp_dir, \
    output_from_subprocess_exists, keep_only_first_contig, get_fasta_lengths, \
    file_len, check_version_from_init, check_version_from_cmd

## GLOBALS
SAMTOOLS_MIN_VERSION = '1.3.1'
#################################### functions ###############################


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="Given regions from riboSnag, assembles the mapped reads",
        add_help=False)  # to allow for custom help
    parser.add_argument("seed_dir", action="store",
                        help="path to roboSnag results directory")

    # taking a hint from http://stackoverflow.com/questions/24180527
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-F", "--fastq1", dest='fastq1', action="store",
                               help="forward fastq reads, can be compressed",
                               type=str, default="", required=True)
    requiredNamed.add_argument("-R", "--fastq2", dest='fastq2', action="store",
                               help="reverse fastq reads, can be compressed",
                               type=str, default="", required=True)
    requiredNamed.add_argument("-r", "--reference_genome",
                               dest='reference_genome',
                               action="store", default='', type=str,
                               help="fasta reference, used to estimate " +
                               "insert sizes, and compare with QUAST",
                               required=True)
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-S", "--fastq_single", dest='fastqS',
                          action="store",
                          help="single fastq reads", type=str, default="")
    optional.add_argument("-n", "--experiment_name", dest='exp_name',
                          action="store",
                          help="prefix for results files; " +
                          "default: %(default)s",
                          default="riboSeed", type=str)
    optional.add_argument("-m", "--method_for_map", dest='method',
                          action="store",
                          help="available mappers: smalt; " +
                          "default: %(default)s",
                          default='smalt', type=str)
    optional.add_argument("-c", "--cores", dest='cores', action="store",
                          default=None, type=int,
                          help="cores for multiprocessing workers" +
                          "; default: %(default)s")
    optional.add_argument("-k", "--kmers", dest='kmers', action="store",
                          default="21,33,55,77,99,127", type=str,
                          help="kmers used for final assembly" +
                          ", separated by commas; default: %(default)s")
    optional.add_argument("-p", "--pre_kmers", dest='pre_kmers',
                          action="store",
                          default="21,33,55", type=str,
                          help="kmers used during seeding assemblies, " +
                          "separated bt commas" +
                          "; default: %(default)s")
    optional.add_argument("-g", "--min_growth", dest='min_growth',
                          action="store",
                          default=0, type=int,
                          help="skip remaining iterations if contig doesnt " +
                          "extend by --min_growth. if 0, ignore" +
                          "; default: %(default)s")
    optional.add_argument("-s", "--min_score_SMALT", dest='min_score_SMALT',
                          action="store",
                          default=None, type=int,
                          help="min score forsmalt mapping; inferred from " +
                          "read length" +
                          "; default: inferred")
    optional.add_argument("--include_shorts", dest='include_shorts',
                          action="store_true",
                          default=False,
                          help="if assembled contig is smaller than  " +
                          "--min_assembly_len, contig will still be included" +
                          " in assembly; default: inferred")
    optional.add_argument("-a", "--min_assembly_len", dest='min_assembly_len',
                          action="store",
                          default=6000, type=int,
                          help="if initial SPAdes assembly largest contig " +
                          "is not at least as long as --min_assembly_len, " +
                          "exit. Set this to the length of the seed " +
                          "sequence; if it is not achieved, seeding across " +
                          "regions will likely fail; default: %(default)s")
    optional.add_argument("--paired_inference", dest='paired_inference',
                          action="store_true", default=False,
                          help="if --paired_inference, mapped read's " +
                          "pairs are included; default: %(default)s")
    optional.add_argument("--subtract", dest='subtract', action="store_true",
                          default=False,
                          help="if --subtract, reads aligned " +
                          "to each reference will not be aligned to future " +
                          "iterations.  Probably you shouldnt do this" +
                          "unless you really happen to want to")
    optional.add_argument("--keep_unmapped", dest='keep_unmapped',
                          action="store_true", default=False,
                          help="if --keep_unmapped, fastqs are generated " +
                          "containing unmapped reads; default: %(default)s")
    optional.add_argument("--ref_as_contig", dest='ref_as_contig',
                          action="store", default="untrusted", type=str,
                          choices=["None", "trusted", "untrusted"],
                          help="if 'trusted', SPAdes will  use the seed " +
                          "sequences as a --trusted-contig; if 'untrusted', " +
                          "SPAdes will treat as --untrusted-contig. if '', " +
                          "seeds will not be used during assembly. " +
                          "See SPAdes docs; default: %(default)s")
    optional.add_argument("--no_temps", dest='no_temps', action="store_true",
                          default=False,
                          help="if --no_temps, mapping files will be " +
                          "removed after all iterations completed; " +
                          " default: %(default)s")
    optional.add_argument("--skip_control", dest='skip_control',
                          action="store_true",
                          default=False,
                          help="if --skip_control, no SPAdes-only de novo " +
                          "assembly will be done; default: %(default)s")
    optional.add_argument("-i", "--iterations", dest='iterations',
                          action="store",
                          default=3, type=int,
                          help="if iterations>1, multiple seedings will " +
                          "occur after assembly of seed regions; " +
                          "if setting --target_len, seedings will continue " +
                          "until --iterations are completed or target_len"
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
    optional.add_argument("--DEBUG", dest='DEBUG', action="store_true",
                          default=False,
                          help="if --DEBUG, test data will be " +
                          "used; default: %(default)s")
    optional.add_argument("--DEBUG_multi", dest='DEBUG_multiprocessing',
                          action="store_true",
                          default=False,
                          help="if --DEBUG_multiprocessing, runs seeding in " +
                          "single loop instead of a multiprocessing pool" +
                          ": %(default)s")
    optional.add_argument("--smalt_scoring", dest='smalt_scoring',
                          action="store",
                          default="match=1,subst=-4,gapopen=-4,gapext=-3",
                          help="submit custom smalt scoring via smalt -S " +
                          "scorespec option; default: %(default)s")
    # had to make this explicitly to call it a faux optional arg
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")

    ##TODO  Make these check a config file
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
    optional.add_argument("--quast_exe", dest="quast_exe",
                          action="store", default="quast.py",
                          help="Path to quast executable; " +
                          "default: %(default)s")
    args = parser.parse_args()
    return args


def check_smalt_full_install(smalt_exe, logger=None):
    smalttestdir = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                "sample_data",
                                "smalt_test", "")
    if logger is None:
        raise ValueError("Must Use Logging")
    logger.debug("looking for smalt test dir: {0}".format(
        smalttestdir))
    if not os.path.exists(smalttestdir):
        raise FileNotFoundError("cannot find smalt_test dir containing " +
                                "files to verify bambamc install!")
    ref = os.path.join(smalttestdir, "ref_to_test_bambamc.fasta")
    index = os.path.join(smalttestdir, "test_index")
    test_bam = os.path.join(smalttestdir, "test_mapping.bam")
    test_reads = os.path.join(smalttestdir, "reads_to_test_bambamc.fastq")
    testindexcmd = str("{0} index {1} {2}".format(smalt_exe, index, ref))
    testmapcmd = str("{0} map -f bam -o {1} {2} {3}".format(smalt_exe,
                                                            test_bam,
                                                            index,
                                                            test_reads))
    logger.debug("testing instalation of smalt and bambamc")
    for i in [testindexcmd, testmapcmd]:
        try:
            logger.debug(i)
            subprocess.run([i],
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        except:
            raise ValueError("Error running test to check bambamc lib is " +
                             "installed! See github.com/gt1/bambamc " +
                             "and the smalt install guide for more details." +
                             "https://sourceforge.net/projects/smalt/files/")
    os.remove(test_bam)
    os.remove(str(index + ".sma"))
    os.remove(str(index + ".smi"))


def estimate_distances_smalt(outfile, smalt_exe, cores, ref_genome,
                             fastq1, fastq2, logger=None):
    """Given fastq pair and a reference, returns path to distance estimations
    used by smalt to help later with mapping. if one already exists,
    return path to it.
    """
    if not os.path.exists(outfile):
        # Index reference for sampling to get PE distances
        if logger:
            logger.info("Estimating insert distances with SMALT")
        # index with default params for genome-sized sequence
        refindex_cmd = str(smalt_exe + " index -k {0} -s {1} {2} " +
                           "{3}").format(20, 10, outfile, ref_genome)
        refsample_cmd = str(smalt_exe + " sample -n {0} -o {1} {2} {3} " +
                            "{4}").format(cores,
                                          outfile,
                                          outfile,
                                          fastq1,
                                          fastq2)
        if logger:
            logger.info("Sampling and indexing {0}".format(
                ref_genome))
        for cmd in [refindex_cmd, refsample_cmd]:
            if logger:
                logger.debug("\t command:\n\t {0}".format(cmd))
            subprocess.run(cmd,
                           shell=sys.platform != "win32",
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           check=True)
    else:
        if logger:
            logger.info("using existing reference file")
        pass
    return outfile

# TODO reimplement bwa, but use BWA-SW instead of MEM
# def map_to_ref_map_mem(ref, fastq_read1, fastq_read2, map_results_prefix,
#                        cores, stdout, stderr, kseed=19):
#     """outputs in bam format to play nice with alternative option, bwa aln
#         Since 0.5:
#         -- removed unpaired penalty (for obvious reasons)|| -U def 9, now 0
#         -- increased mismatch penalty || -B def 4, now 8
#         -- open gap penalty decrease from 6 to 0
#     """
#     print('######  Running BWA MEM...' +
#           str(datetime.time(datetime.now())).split('.')[0])
#     subprocess.call('bwa index %s' % ref, shell=True, stdout=stdout,
#                      stderr=stderr)
#     try:
#         kseed = int(kseed)
#     except ValueError:
#         raise("k must be numeric")
#     subprocess.call('bwa mem -A 1 -d 20  -U 0 -L 100 -B 100 -a -O 6 -t '+
#                     '%s -k %i %s %s %s > %s.sam' %
#                     (cores, kseed, ref, fastq_read1, fastq_read2,
#                      map_results_prefix),
#                     stdout=stdout, stderr=stderr, shell=True)
#     subprocess.call('samtools view -bhS %s.sam > %s.bam' %
#                     (map_results_prefix, map_results_prefix), shell=True,
#                     stdout=stdout, stderr=stderr)


def map_to_ref_smalt(ref, fastq_read1, fastq_read2,
                     distance_results,
                     map_results_prefix, cores, samtools_exe,
                     smalt_exe, score_minimum=None, fastq_readS="",
                     read_len=100, step=3, k=5,
                     scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
                     logger=None):
    """run smalt based on pased args
    requires at least paired end input, but can handle an additional library
    of singleton reads. Will not work on just singletons
    """
    if score_minimum is None:
        score_min = int(read_len * .3)
    else:
        score_min = score_minimum
    logger.debug(str("mapping with smalt using a score min of " +
                     "{0}").format(score_min))
    cmdindex = str("{3} index -k {0} -s {1} {2} {2}").format(
        k, step, ref, smalt_exe)
    cmdmap = str('{7} map -l pe -S {8} ' +
                 '-m {0} -n {1} -g {2} -f bam -o {3}_pe.bam {4} {5} ' +
                 '{6}').format(score_min, cores, distance_results,
                               map_results_prefix, ref, fastq_read1,
                               fastq_read2, smalt_exe, scoring)
    # cmdview = str('{0} view -bhS {1}.sam > {1}_pe.bam').format(
    #     samtools_exe, map_results_prefix)
    # smaltcommands = [cmdindex, cmdmap, cmdview]
    smaltcommands = [cmdindex, cmdmap]
    if fastq_readS != "":
        cmdindexS = str('{0} index -k {1} -s {2} {3} {3}').format(
            smalt_exe, k, step, ref)
        cmdmapS = str("{7} map -S {6} " +
                      "-m {0} -n {1} -g {2} -f bam -o {3}S.bam {4} " +
                      "{5}").format(score_min, cores, distance_results,
                                    map_results_prefix, ref, fastq_readS,
                                    scoring, smalt_exe)
        # cmdviewS = str('{0} view -bhS {1}S.sam > ' +
        #                '{1}S.bam').format(samtools_exe, map_results_prefix)
        cmdmergeS = str('{0} merge -f  {1}.bam {1}_pe.bam ' +
                        '{1}S.bam').format(samtools_exe, map_results_prefix)
        # smaltcommands.extend([cmdindexS, cmdmapS, cmdviewS, cmdmergeS])
        smaltcommands.extend([cmdindexS, cmdmapS, cmdmergeS])
    else:
        cmdmerge = str("{0} view -bh {1}_pe.bam >" +
                       "{1}.bam").format(samtools_exe, map_results_prefix)
        smaltcommands.extend([cmdmerge])
    logger.info("running SMALT:")
    logger.debug("with the following SMALT commands:")
    for i in smaltcommands:
        logger.debug(i)
        subprocess.run(i, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
    if fastq_readS != '':
        logger.info(str("Singleton mapped reads: " +
                        get_number_mapped(str(map_results_prefix + "S.bam"),
                                          samtools_exe=samtools_exe)))
    logger.info(str("PE mapped reads: " +
                    get_number_mapped(str(map_results_prefix + "_pe.bam"),
                                      samtools_exe=samtools_exe)))
    logger.info(str("Combined mapped reads: " +
                    get_number_mapped(str(map_results_prefix + ".bam"),
                                      samtools_exe=samtools_exe)))


def convert_bams_to_fastq(map_results_prefix,
                          fastq_results_prefix,
                          keep_unmapped, samtools_exe, logger=None):
    """ return 6 paths: unmapped F, R, S and mapped F, R. S
    given the prefix for the mapped bam files, write out mapped (and optionally
    unmapped) reads to fasta files
    """
    convert_cmds = []
    bams = ["_unmapped", "_mapped"]
    for i in bams:
        if not os.path.exists("%s%s.bam" % (map_results_prefix, i)):
            if i == '_unmapped' and not keep_unmapped:
                continue
            else:
                if logger:
                    logger.error(str("No {0}{1}.bam file" +
                                     "found").format(map_results_prefix, i))
                raise FileNotFoundError("Cannot find bam mapping; files" +
                                        "must have '_mapped.bam' prefix")
        samfilter = \
            str(samtools_exe + " fastq {0}{1}.bam -1 {2}{1}1.fastq -2 " +
                "{2}{1}2.fastq -s {2}{1}S.fastq").format(map_results_prefix, i,
                                                         fastq_results_prefix)
        convert_cmds.append(samfilter)
    if logger:
        logger.debug("running the following commands to extract reads:")
    for i in convert_cmds:
        if logger:
            logger.debug(i)
        subprocess.run(i, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
    if keep_unmapped:
        fnames = []
        for bam_idx in (0, 1):
            for suffix in ('1', '2', 'S'):
                fnames.append("{0}{1}{2}.fastq".format(fastq_results_prefix,
                                                       bams[bam_idx], suffix))
            return fnames
    else:
        return(None,  # unmapped forward
               None,  # unmapped reverse
               None,  # unmapped Single
               "%s%s1.fastq" % (fastq_results_prefix, bams[1]),  # mapped fwd
               "%s%s2.fastq" % (fastq_results_prefix, bams[1]),  # mapped rev
               "%s%sS.fastq" % (fastq_results_prefix, bams[1]))  # mapped s


def run_spades(output, ref, ref_as_contig, pe1_1='', pe1_2='', pe1_s='',
               as_paired=True, keep_best=True, prelim=False,
               groom_contigs='keep_first',
               k="21,33,55,77,99", seqname='', spades_exe="spades.py",
               logger=None):
    """return path to contigs
    wrapper for common spades setting for long illumina reads
    ref_as_contig should be either blank, 'trusted', or 'untrusted'
    prelim flag is True, only assembly is run, and without coverage correction
    #TODO
    the seqname variable is used only for renaming the resulting contigs
    during iterative assembly.  It would be nice to inheirit from "ref",
    but that is changed with each iteration. This should probably be addressed
    before next major version change
    """
    if logger is None:
        raise ValueError("this must be used with a logger!")
    if groom_contigs not in ['keep_first', 'consensus']:
        raise ValueError("groom_contigs option must be either 'keep_first' " +
                         "or 'consensus'")
    if seqname == '':
        seqname = ref
    kmers = k  # .split[","]
    success = False
    #  prepare reference, if being used
    if not ref_as_contig is None:
        alt_contig = str("--%s-contigs %s" % (ref_as_contig, ref))
    else:
        alt_contig = ''
    # prepare read types, etc
    if as_paired and pe1_s != "":  # for libraries with both
        singles = str("--pe1-s %s " % pe1_s)
        pairs = str("--pe1-1 %s --pe1-2 %s " % (pe1_1, pe1_2))
    elif as_paired and pe1_s == "":  # for libraries with just PE
        singles = ""
        pairs = str("--pe1-1 %s --pe1-2 %s " % (pe1_1, pe1_2))
    elif pe1_s == "":  # for libraries treating paired ends as two single-end libs
        singles = ''
        pairs = str("--pe1-s %s --pe2-s %s " % (pe1_1, pe1_2))
    else:  # for 3 single end libraries
        singles = str("--pe1-s %s " % pe1_s)
        pairs = str("--pe2-s %s --pe3-s %s " % (pe1_1, pe1_2))
    reads = str(pairs + singles)
#    spades_cmds=[]
    if prelim:
        prelim_cmd =\
            str("{4}  --only-assembler --cov-cutoff off --sc --careful -k" +
                " {0} {1} {2} -o {3}").format(kmers, reads, alt_contig,
                                              output, spades_exe)
        logger.debug("Running SPAdes command:\n{0}".format(prelim_cmd))
        subprocess.run(prelim_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        success = output_from_subprocess_exists(os.path.join(output,
                                                             "contigs.fasta"))
        if groom_contigs == "keep_first" and success:
            logger.info("reserving first contig")
            try:
                keep_only_first_contig(str(os.path.join(output,
                                                        "contigs.fasta")),
                                       newname=os.path.splitext(
                                           os.path.basename(seqname))[0])
            except Exception as f:
                logger.error(f)
                raise f
        elif success and groom_contigs == "consensus":
            # get consensus; copy ref for starters to double check later
            contigs_backup = copy_file(current_file=ref, dest_dir=output,
                                       name=str("backedup_contigs.fasta"),
                                       overwrite=True, logger=logger)

            logger.debug("copying {0} to {0} as a backup to self-test " +
                         " consensus".format(ref, contigs_backup))
            # make pileup
            pileupcmd = str("smalt index {0} {0} ; smalt map {0} {1} | " +
                            "samtools sort - | samtools mpileup -f {0} - -o " +
                            "{2}").format(ref,
                                          os.path.join(output,
                                                       'contigs.fasta'),
                                          os.path.join(output,
                                                       'contigs_pileup.txt'))
            subprocess.run(pileupcmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
            # test pileup
            try:
                pileup = check_samtools_pileup(
                    os.path.join(output, 'contigs_pileup.txt'))
            except Exception as e:
                logger.error(e)
                sys.exit(1)
            # parse pileup
            try:
                consensus = reconstruct_seq(refpath=ref, pileup=pileup,
                                            verbose=False, veryverb=False,
                                            logger=logger)
            except Exception as e:
                logger.error(e)
                sys.exit(1)

            with open(os.path.join(output, 'contigs.fasta'), 'w') as new_seqs:
                seqrec = Seq(consensus, IUPAC.IUPACAmbiguousDNA())
                SeqIO.write(SeqRecord(seqrec,
                                      id="contigs_consensus_riboSeed",
                                      description=""), new_seqs, 'fasta')
        else:
            logger.warning("No output from SPAdes this time around")
    else:
        spades_cmd = str(spades_exe + " --careful -k {0} {1} {2} -o " +
                         "{3}").format(kmers, reads, alt_contig, output)
        logger.info("Running the following command:\n{0}".format(spades_cmd))
        subprocess.run(spades_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
        # not check=True; dont know spades return codes
        success = output_from_subprocess_exists(os.path.join(output,
                                                             "contigs.fasta"))
    return("{0}contigs.fasta".format(os.path.join(output, "")), success)


def check_samtools_pileup(pileup):
    """checks that the file is a valid samtools pileup. outputs a list
    """
    res = []
    try:
        with open(pileup, "r") as file:
            for f in file:
                f = re.split(r'\t+', f.strip())
                if len(f) != 6:
                    print("6 columns not found")
                if not isinstance(int(f[1]), int):
                    print("second column isnt int")
                res.append(f)
    except:
        raise ValueError("Error with reading pileup file")
    return res


def reconstruct_seq(refpath, pileup, verbose=True, veryverb=False,
                    logger=None):
    """ This is a bit of a mess, to say the least.  Given a list from
    check samtools pileup, and a reference fasta, this reconstructs ambiguous
    regions
    """
    if logger is None:
        raise ValueError("Logger needed for the 'reconstruct_seqs' function!")
    logger.warning("This function is sketchy at best. Here be dragons!")
    if verbose:
        logger.debug(str("reconstucting consensus sequence " +
                         "from {0} and pileup").format(refpath))
    seqfile = SeqIO.parse(open(refpath, "r"), "fasta")
    for i in seqfile:
        ref = str(i.seq)
    ref = "${0}".format(ref)  # make 0 based
    new = ""
    skip = 0
    indels = 0
    N_deletions, N_insertions = 0, 0
    insert = re.compile('[,\\.]{0,1}\\+[0-9]+[ACGTNacgtn]+')
    delete = re.compile('[\\.,]{0,1}-[0-9]+[ACGTNacgtn]+')
    j = 0  # counter for pileup
    for i in range(0, len(ref)):  # counter for ref index
        if veryverb and verbose:
            try:
                logger.debug("ref. index: %i\n pile index: %i" % (i, j))
                logger.debug(ref[i])
                logger.debug(pileup[j])
            except:
                pass
        # This is how we handle deletions; decrement skip, and next iteration
        if skip > 0:
            skip = skip - 1
            if verbose:
                logger.debug("skipping {0}".format(i))
        # if j is greater then length of pileup, go with ref.
        # This should avoid out of range issues
        elif j > len(pileup) - 1:
            new = "".join([new, ref[i]])
        #  if index isnt in second col of pileup, skip, filling with ref
        # note because the reference is now zero base, no correction needed
        elif i != int((pileup[j][1])):
            if verbose:
                logger.debug("no entry in pileup for %i" % i)
            new = "".join([new, ref[i]])
            # this should keep pileup counter the same when
            j = j - 1
        # if (N in ref, or pilup differs from ref), and
        # (pilup has single value, or all values same) and
        # pileup char isnt $*, go with pileup
        # *(start char is ^W, which is two chars, breaks 3rd line of conditions
        # NOTE: lowercase letters converted to upper, because orientation
        #       is already handled by samtools.
        elif (ref[i] == "N" or pileup[j][4][0] != ref[i]) and \
             (len(pileup[j][4]) == 1 or
              all(x == pileup[j][4][0] for x in list(pileup[j][4]))) and \
            pileup[j][4][0] not in [",", ".", "^", "$"]:
            new = "".join([new, pileup[j][4][0].upper()])  # append  upper
        # This is tp handle insetions;  could use a lamda?
        elif re.match(insert, pileup[j][4]) is not None and \
            all([hits == re.findall(insert, pileup[j][4])[0] for hits in
                 re.findall(insert, pileup[j][4])]):
            if verbose:
                logger.debug("found insert!")
            insert_seq = re.search('[ACGTNacgtn]+', pileup[j][4]).group(0)
            insert_N = int(re.search('[0-9]+', pileup[j][4]).group(0))
            if not len(insert_seq) == insert_N:
                raise ValueError("error parsing insert")
            new = "".join([new, insert_seq])
            indels = indels + insert_N
            N_insertions = N_insertions + insert_N
        # deletions
        elif re.match(delete, pileup[j][4]) is not None and \
            all([hits == re.findall(insert, pileup[j][4])[0] for hits in
                 re.findall(insert, pileup[j][4])]):
            if verbose:
                logger.debug("found deletion! {0}".format(pileup[j][4]))
            delete_N = int(re.search('[0-9]+', pileup[j][4]).group(0))
            skip = delete_N
            indels = indels + delete_N
            N_deletions = N_deletions + delete_N
        # Most cases fall in this category
        elif pileup[j][4][0] in [",", ".", "^", "$"] or \
            not all(x == pileup[j][4][0] for x in list(pileup[j][4])):
            if verbose:
                logger.debug("using ref")
            new = "".join([new, ref[i]])
        else:
            raise ValueError("Error parsing pileup: Case Not covered!")
        j = j + 1  # increment the pileup counter
    if verbose:
        logger.info(str("total indels: {0}\n\tdeletions {1}\n\tinsetions: " +
                        "{2}").format(indels, N_deletions, N_insertions))
    else:
        print(str("total indels: {0}\n\tdeletions {1}\n\tinsetions: " +
                  "{2}").format(indels, N_deletions, N_insertions))
    return new[1:]  # [1:] gets rid of starting dollar character


def make_quick_quast_table(pathlist, write=False, writedir=None, logger=None):
    """This skips any fields not in first report, for better or worse...
    """
    if logger is None:
        raise ValueError("Logging must be enabled for make_quick_quast_table")
    if not isinstance(pathlist, list):
        logger.warning("paths for quast reports must be in a list!")
        return None
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
            except Exception as e:
                raise e("error parsing %s", i)
        else:
            report_list = []
            try:
                with open(i, "r") as handle:
                    for dex, line in enumerate(handle):
                        row, val = line.strip().split("\t")
                        report_list.append([row, val])
                    logger.debug("report list: %s", str(report_list))
                    for k, v in mainDict.items():
                        if k in [x[0] for x in report_list]:
                            mainDict[k].append(
                                str([x[1] for x in
                                     report_list if x[0] == k][0]))
                        else:
                            mainDict[k].append("XX")
            except Exception as e:
                raise e("error parsing %s", i)
        counter = counter + 1
    logger.info(str(mainDict))
    if write:
        if writedir is None:
            logger.warning("no output dir, cannot write!")
            return mainDict
        try:
            with open(os.path.join(
                    writedir, "combined_quast_report.tsv"), "w") as outfile:
                for k, v in sorted(mainDict.items()):
                    logger.debug("{0}\t{1}\n".format(k, str("\t".join(v))))
                    outfile.write("{0}\t{1}\n".format(
                        str(k), str("\t".join(v))))
        except Exception as e:
            raise e("Error wrting out combined quast report")

    return mainDict


def main(fasta, num, results_dir, exp_name, mauve_path, map_output_dir, method,
         reference_genome, fastq1, fastq2, fastqS, ave_read_length, cores,
         subtract_reads, ref_as_contig, fetch_mates, keep_unmapped_reads,
         paired_inference, smalt_scoring, min_growth, max_iterations, kmers,
         no_temps, distance_estimation, proceed_to_target, target_len,
         score_minimum, min_contig_len, include_short_contigs):
    """
    process each fasta seed to parallelize time comsuming parts
    """
    prelog = "{0}-{1}:".format("SEED", num)
    logger.info("%s processing %s", prelog, fasta)
    logger.info("%s item %i of %i", prelog, fastas.index(fasta) + 1, nfastas)
    spades_dir = str(results_dir + "SPAdes_" +
                     os.path.split(fasta)[1].split(".fasta")[0])
    mapping_dir = str(map_output_dir + "mapping_" +
                      os.path.split(fasta)[1].split(".fasta")[0])
    logger.debug("%s output dirs: \n%s\n%s", prelog, spades_dir, mapping_dir)
    #  make appropriate output directories
    for i in [spades_dir, mapping_dir]:
        if not os.path.isdir(i):
            os.makedirs(i)
    map_results_prefix = os.path.join(
        mapping_dir, str(exp_name + "_" +
                         os.path.split(fasta)[1].split(".fasta")[0]))
    fastq_results_prefix = os.path.join(
        results_dir, str(exp_name + "_" +
                         os.path.split(fasta)[1].split(".fasta")[0]))
    logger.debug("%s copying seed file to mapping directory", prelog)
    new_reference = copy_file(current_file=fasta, dest_dir=mapping_dir,
                              name='', overwrite=False, logger=logger)

    seed_len = get_fasta_lengths(fasta)[0]
    # set proceed_to_target params
    if proceed_to_target:
        if target_len > 0 and 5 > target_len:
            target_seed_len = int(target_len * seed_len)
        elif target_len > 50:
            target_seed_len = int(target_len)
        else:
            logger.error("%s invalid taget length provided; must be given" +
                         " as fraction of total length or as an absolute " +
                         "number of base pairs greater than 50", prelog)
            sys.exit(1)
    else:
        pass

    this_iteration = 1  # counter for iterations.
    # perform iterative mapping
    proceed = True  # this is set to false in while loop if spades fails
    keep_contig = True  # this becomes false if assembly fails,
    while this_iteration <= max_iterations and proceed:
        # if not the first round, replace ref with extended contigs
        if this_iteration != 1:
            logger.info("%s copying contigs for next iteration of assembly",
                        prelog)
            new_reference = copy_file(current_file=new_reference,
                                      dest_dir=mapping_dir,
                                      name=str(os.path.basename(fasta) +
                                               "_iter_" + str(this_iteration) +
                                               ".fasta"),
                                      overwrite=True, logger=logger)
        logger.info("%s Iteration %i of %i, item %i out of %i",
                    prelog, this_iteration, max_iterations,
                    fastas.index(fasta) + 1, len(fastas))
        map_to_ref_smalt(ref=new_reference,  # ref_genome=reference_genome,
                         fastq_read1=fastq1, fastq_read2=fastq2,
                         fastq_readS=fastqS, read_len=average_read_length,
                         map_results_prefix=map_results_prefix, cores=cores,
                         step=3, k=5,
                         distance_results=distance_estimation,
                         score_minimum=score_minimum,
                         scoring=smalt_scoring, smalt_exe=args.smalt_exe,
                         samtools_exe=args.samtools_exe, logger=logger)
        extract_mapped_and_mappedmates(map_results_prefix,
                                       fetch_mates=fetch_mates,
                                       samtools_exe=args.samtools_exe,
                                       keep_unmapped=keep_unmapped_reads,
                                       logger=logger)

        logger.info("%s Converting mapped results to fastqs", prelog)
        try:
            new_fastq1, new_fastq2, new_fastqS, \
                mapped_fastq1, mapped_fastq2, mapped_fastqS = \
                    convert_bams_to_fastq(map_results_prefix,
                                          fastq_results_prefix,
                                          keep_unmapped=keep_unmapped_reads,
                                          logger=logger,
                                          samtools_exe=args.samtools_exe)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        if subtract_reads and keep_unmapped_reads:
            logger.warning("%s using reduced reads for next iteration", prelog)
            fastq1, fastq2, fastqS = new_fastq1, new_fastq2, new_fastqS
        logger.info("%s Running SPAdes", prelog)
        try:
            contigs_path, proceed = \
                run_spades(pe1_1=mapped_fastq1, pe1_2=mapped_fastq2,
                           pe1_s=mapped_fastqS, prelim=True,
                           as_paired=paired_inference,
                           groom_contigs="keep_first",
                           output=spades_dir, keep_best=True,
                           ref=new_reference, ref_as_contig=ref_as_contig,
                           k=kmers, logger=logger,
                           seqname=fasta, spades_exe=args.spades_exe)
        except Exception as e:
            logger.error("%s SPAdes error:", prelog)
            logger.error(e)
            sys.exit(1)
        if not proceed:
            logger.warning("%s Assembly failed: no spades output for %s",
                           prelog, os.path.basename(fasta))

        # compare lengths of reference and freshly assembled contig
        contig_len = get_fasta_lengths(contigs_path)[0]
        ref_len = get_fasta_lengths(new_reference)[0]
        contig_length_diff = contig_len - ref_len
        logger.info("%s Seed length: %i", prelog, seed_len)
        if proceed_to_target:
            logger.info("Target length: {0}".format(target_seed_len))
        logger.info("%s Length of this iteration's longest contig: %i",
                    prelog, contig_len)
        if this_iteration != 1:
            logger.info("%s Length of previous longest contig: %i",
                        prelog, ref_len)
            logger.info("%s The new contig differs from the previous " +
                        "iteration by %i bases", prelog, contig_length_diff)
        else:
            logger.info("%s The new contig differs from the reference " +
                        "seed by %i bases", prelog, contig_length_diff)

        # This cuts failing assemblies short
        if this_iteration == 1 and min_contig_len > contig_len:
            logger.warning("The first iteration's assembly's best contig" +
                           " is not greater than length set by " +
                           "--min_assembly_len. Assembly will likely fail if" +
                           " the contig does not meet the length of the seed")
            if include_short_contigs:
                logger.warning("Continuing, but if this occurs for more " +
                               "than one seed, we reccommend  you abort and " +
                               "retry with longer seeds, a different ref, " +
                               "or re-examine the riboSnag clustering")
            else:
                keep_contig = False  # flags contig for exclusion
            this_iteration = max_iterations + 1  # skip remaining iterations
        else:
            pass
        # This is a feature that is supposed to help skip unneccesary
        # iterations. If the difference is negative (new contig is shorter)
        # continue, as this may happen (especially in first mapping if
        # reference is not closely related to Sample), continue to map.
        # If the contig length increases, but not as much as min_growth,
        # skip future iterations
        if contig_length_diff > 0 and contig_length_diff < min_growth and \
           min_growth > 0:  # ie, ignore by default
            logger.info("the length of the new contig was only 0bp changed " +
                        "from previous iteration; skipping future iterations")
            this_iteration = max_iterations + 1  # skip remaining iterations
        # if continuing til reaching the target lenth of the seed
        elif proceed_to_target and contig_len >= target_seed_len:
            logger.info("target length threshold! has been reached; " +
                        "skipping future iterations")
            this_iteration = max_iterations + 1  # skip remaining iterations
        else:
            this_iteration = this_iteration + 1

        # use contigs_path as new reference
        new_reference = contigs_path

    #  Now, after iterative seeding
    if not keep_contig:
        logger.warning("Excluding contig seeded by %s!", fasta)
        return(1)
    else:
        try:
            contigs_new_path = copy_file(current_file=contigs_path,
                                         dest_dir=mauve_path,
                                         name=str(os.path.basename(fasta) +
                                                  "_final_iter_" +
                                                  str(this_iteration) + ".fasta"),
                                         logger=logger)
        except:
            logger.warning(str("no contigs moved for {0}! Check  SPAdes log " +
                               "in results dir if worried").format(fasta))
    logger.debug("moved {0} to {1}".format(contigs_path, contigs_new_path))
    if no_temps:
        logger.info("removing temporary files from {0}".format(mapping_dir))
        clean_temp_dir(mapping_dir)
    if proceed:
        return 0
    else:
        return 1


#%%
if __name__ == "__main__":
    args = get_args()
    # smalt checks this to estimate PE distance; see estimate_smalt_distances
    mapped_genome_sam = "genome_distance_est.sam"
    # allow user to give relative paths
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root)
    except OSError:
        raise OSError(str("Output directory already exists"))
    map_output_dir = os.path.join(output_root, 'map', "")
    results_dir = os.path.join(output_root, 'results', "")
    mauve_dir = os.path.join(output_root, 'results', "mauve", "")
    t0 = time.time()
    log_path = os.path.join(output_root,
                            str("{0}_riboSeed_log.txt".format(
                                time.strftime("%Y%m%d%H%M"))))
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=log_path,
                            name=__name__)
    package_init = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "__init__.py")

    # log version of riboSeed, commandline options, and all settings
    logger.info("riboSeed pipeine package version {0}".format(
        check_version_from_init(init_file=package_init, min_version="0.0.0")))
    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    if args.cores is None:
        args.cores = multiprocessing.cpu_count()
        logger.info("Using %i cores", multiprocessing.cpu_count())
    logger.debug(str("\noutput root {0}\nmap_output_dir: {1}\nresults_dir: " +
                     "{2}\n").format(output_root, map_output_dir, results_dir))

    # Cannot set Nonetype objects via commandline directly and I dont want None
    # to be default dehaviour, so here we convert 'None' to None.
    # I have not moral compass
    if args.ref_as_contig == 'None':
        args.ref_as_contig = None

    # TODO Look into resupporting bwa, as it plays better with BAM files,
    # though smalt beats it on the overhangs. The main issue is that BWA mem
    # doesnt work with overhangs well, and bwasw is slower than smalt.
    # if args.method is not in  ["smalt", "bwa"]:
    #     logger.error("'smalt' and  'bwa' only methods currently supported")
    #     sys.exit(1)
    if args.method not in ["smalt"]:
        logger.error("'smalt' only method currently supported")
        sys.exit(1)
    logger.debug("checking for installations of all required external tools")
    executables = [args.samtools_exe, args.spades_exe, args.quast_exe]
    if args.method == "smalt":
        executables.append(args.smalt_exe)
    # elif args.method == "bwa":
    #     executables.append(args.bwa_exe)
    else:
        logger.error("Mapping method not found!")
        sys.exit(1)
    logger.debug(str(executables))
    test_ex = [check_installed_tools(x, logger=logger) for x in executables]
    if all(test_ex):
        logger.debug("All needed system executables found!")
        logger.debug(str([shutil.which(i) for i in executables]))

    # check samtools verison
    try:
        samtools_verison = check_version_from_cmd(
            exe='samtools',
            cmd='',
            line=3,
            pattern=r"\s*Version: (?P<version>[^(]+)",
            where='stderr',
            min_version=SAMTOOLS_MIN_VERSION)
    except Exception as e:
        logger.error(e)
        sys.exit(1)
    logger.debug("samtools version: %s", samtools_verison)
    # check bambamc is installed proper if using smalt
    if args.method == "smalt":
        try:
            check_smalt_full_install(smalt_exe=args.smalt_exe, logger=logger)
        except Exception as e:
            logger.error(e)
            sys.exit(1)

    # check equal length fastq.  This doesnt actually check propper pairs
    if file_len(args.fastq1) != file_len(args.fastq2):
        logger.error("Input Fastq's are of unequal length! Try " +
                     "fixing with this script: " +
                     "github.com/enormandeau/Scripts/fastqCombinePairedEnd.py")
        sys.exit(1)

    # if the target_len is set. set needed params
    if args.target_len is not None:
        if not args.target_len > 0 or not isinstance(args.target_len, float):
            logger.error("--target_len is set to invalid value! Must be a " +
                         "decimal greater than zero, ie where 1.1 would be " +
                         "110% of the original sequence length.")
            sys.exit(1)
        elif args.target_len > 5 and 50 > args.target_len:
            logger.error("We dont reccommend seeding to lengths greater than" +
                         "5x original seed length. Try between 0.5 and 1.5." +
                         "  If you are setting a target number of bases, it " +
                         " must be greater than 50")
            sys.exit(1)
        else:
            proceed_to_target = True
    else:
        proceed_to_target = False

    for i in [map_output_dir, results_dir, mauve_dir]:
        make_outdir(i)
    average_read_length = get_ave_read_len_from_fastq(args.fastq1,
                                                      N=36, logger=logger)
    fastq_results_prefix = os.path.join(results_dir, args.exp_name)

    fastas = [os.path.join(args.seed_dir, x) for
              x in os.listdir(os.path.join(args.seed_dir, "")) if
              x.endswith('.fasta')]
    if len(fastas) == 0:
        logger.error(str("no files found in {0} ending with " +
                         "'.fasta'").format(args.seed_dir))

    nfastas = len(fastas)
    logger.debug(fastas)

    ### if using smalt (which you are), check for mapped reference
    if args.method == 'smalt':
        path_to_distance_file = os.path.join(results_dir, mapped_genome_sam)
        dist_est = estimate_distances_smalt(outfile=path_to_distance_file,
                                            smalt_exe=args.smalt_exe,
                                            ref_genome=args.reference_genome,
                                            fastq1=args.fastq1,
                                            fastq2=args.fastq2,
                                            cores=args.cores, logger=logger)
    else:
        logger.error("As of v 0.88, only supported mapper is 'smalt'")
        sys.exit(1)

    ### Main function call
    if args.DEBUG_multiprocessing:
        logger.warning("running without multiprocessing!")
        for num, fasta in enumerate(fastas):
            main(fasta=fasta,
                 num=num,
                 results_dir=results_dir,
                 exp_name=args.exp_name,
                 mauve_path=mauve_dir,
                 map_output_dir=map_output_dir,
                 method=args.method,
                 reference_genome=args.reference_genome,
                 fastq1=args.fastq1,
                 fastq2=args.fastq2,
                 fastqS=args.fastqS,
                 ave_read_length=average_read_length,
                 cores=args.cores,
                 subtract_reads=args.subtract,
                 ref_as_contig=args.ref_as_contig,
                 fetch_mates=args.paired_inference,
                 keep_unmapped_reads=args.keep_unmapped,
                 paired_inference=args.paired_inference,
                 smalt_scoring=args.smalt_scoring,
                 min_growth=args.min_growth,
                 max_iterations=args.iterations,
                 kmers=args.pre_kmers,
                 no_temps=args.no_temps,
                 distance_estimation=dist_est,
                 proceed_to_target=proceed_to_target,
                 target_len=args.target_len,
                 score_minimum=args.min_score_SMALT,
                 min_contig_len=args.min_assembly_len,
                 include_short_contigs=args.include_shorts)

    else:
        #  default is now to get cores available for worker and
        pool = multiprocessing.Pool(processes=args.cores)
        results = [pool.apply_async(main, (fasta,),
                                    {"num": num,
                                     "results_dir": results_dir,
                                     "map_output_dir": map_output_dir,
                                     "exp_name": args.exp_name,
                                     "method": args.method,
                                     "reference_genome": args.reference_genome,
                                     "fastq1": args.fastq1,
                                     "fastq2": args.fastq2,
                                     "fastqS": args.fastqS,
                                     "cores": 1,  # args.cores,
                                     "mauve_path": mauve_dir,
                                     "ave_read_length": average_read_length,
                                     "fetch_mates": args.paired_inference,
                                     "keep_unmapped_reads": args.keep_unmapped,
                                     "paired_inference": args.paired_inference,
                                     "subtract_reads": args.subtract,
                                     "ref_as_contig": args.ref_as_contig,
                                     "smalt_scoring": args.smalt_scoring,
                                     "min_growth": args.min_growth,
                                     "max_iterations": args.iterations,
                                     "kmers": args.pre_kmers,
                                     "no_temps": args.no_temps,
                                     "distance_estimation": dist_est,
                                     "proceed_to_target": proceed_to_target,
                                     "target_len": args.target_len,
                                     "score_minimum": args.min_score_SMALT,
                                     "min_contig_len": args.min_assembly_len,
                                     "include_short_contigs": args.include_shorts})
                   for num, fasta in enumerate(fastas)]
        pool.close()
        pool.join()
        logger.info(results)
        logger.info(sum([r.get() for r in results]))

    logging.info("combinging contigs from %s" % mauve_dir)
    new_contig_file = combine_contigs(contigs_dir=mauve_dir,
                                      contigs_name="riboSeedContigs",
                                      logger=logger)
    logger.info("Combined Seed Contigs: {0}".format(new_contig_file))
    logger.info("Time taken to run seeding: %.2fm" % ((time.time() - t0) / 60))
    # logger.info("Time taken to run seeding: %.2fm" % (time.time() - t0) / 60)
    logger.info("\n\n Starting Final Assemblies\n\n")

    quast_reports = []
    if not args.skip_control:
        final_list = ["de_novo", "de_fere_novo"]
    else:
        final_list = ["de_fere_novo"]
    for j in final_list:
        logging.info("\n\nRunning %s SPAdes \n" % j)
        if j == "de_novo":
            assembly_ref = ''
            assembly_ref_as_contig = None
        elif j == "de_fere_novo":
            assembly_ref = new_contig_file
            assembly_ref_as_contig = 'trusted'
        else:
            logger.error("Only valid cases are de novo and de fere novo!")
            sys.exit(1)
        logger.info("Running %s SPAdes" % j)
        output_contigs, \
            final_success = run_spades(pe1_1=args.fastq1, pe1_2=args.fastq2,
                                       output=os.path.join(results_dir, j),
                                       ref=assembly_ref,
                                       ref_as_contig=assembly_ref_as_contig,
                                       prelim=False, keep_best=False,
                                       k=args.kmers, logger=logger)
        if final_success:
            logger.info("Running %s QUAST" % j)
            run_quast(contigs=output_contigs,
                      output=os.path.join(results_dir, str("quast_" + j)),
                      quast_exe=args.quast_exe,
                      threads=args.cores,
                      ref=args.reference_genome,
                      logger=logger)
        quast_reports.append(os.path.join(results_dir, str("quast_" + j),
                                          "report.tsv"))

    if not args.skip_control:
        logger.debug("writing combined quast reports")
        try:
            quast_comp = make_quick_quast_table(quast_reports,
                                                write=True,
                                                writedir=results_dir,
                                                logger=logger)
            for k, v in sorted(quast_comp.items()):
                logger.info("{0}: {1}".format(k, "  ".join(v)))
        except Exception as e:
            logger.warning(e)
        logger.info("Comparing de novo and de fere novo assemblies:")
    # Report that we've finished
    logger.info("Done: %s." % time.asctime())
    logger.info("Time taken: %.2fm" % ((time.time() - t0) / 60))

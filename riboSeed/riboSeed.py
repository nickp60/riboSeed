#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
version 0.9.3

Minor Version Revisions:
 - added min improvement option

Created on Sun Jul 24 19:33:37 2016
The goal of this script will be to use a small fasta sequence to build indices
and use as references for BWA. Mapped sequences will be outputted from BWA
 including the paired-end overhanging regions. The reads will be alligned
 de-novo, and the resulting consensus will be outputted.  This make be done
 iteratively, where the consensus is then indexed and used with BWA. In the
 ideal world, the outputted sequences will be the size of the original sequence
 plus twice the length of the reads, as it should have tails on either side.
@author: nicholas, but I stole some stuff from github/x/stx_subtyping

"""
import gzip
import sys
import time
import re
import errno
import logging
import traceback
import os
import shutil
import argparse
import multiprocessing
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import glob

from pyutilsnrw import utils3_5
from pyutilsnrw.utils3_5 import set_up_logging, make_outdir, \
    make_output_prefix, combine_contigs, run_quast, \
    copy_file, check_installed_tools, get_ave_read_len_from_fastq, \
    get_number_mapped, extract_mapped_and_mappedmates, clean_temp_dir, \
    output_from_subprocess_exists, keep_only_first_contig, get_fasta_lengths

#################################### functions ###############################


def get_args(DEBUG=False):
    parser = argparse.ArgumentParser(
        description="Given regions from riboSnag, assembles the mapped reads")
    parser.add_argument("seed_dir", action="store",
                        help="path to roboSnag results directory")
    parser.add_argument("-F", "--fastq1", dest='fastq1', action="store",
                        help="forward fastq reads, can be compressed",
                        type=str, default="")
    parser.add_argument("-R", "--fastq2", dest='fastq2', action="store",
                        help="reverse fastq reads, can be compressed",
                        type=str, default="")
    parser.add_argument("-S", "--fastq_single", dest='fastqS', action="store",
                        help="single fastq reads", type=str, default="")
    parser.add_argument("-n", "--experiment_name", dest='exp_name',
                        action="store",
                        help="prefix for results files; default: %(default)s",
                        default="riboSeed", type=str)
    parser.add_argument("-m", "--method_for_map", dest='method',
                        action="store",
                        help="availible mappers: smalt; default: %(default)s",
                        default='smalt', type=str)
    parser.add_argument("-c", "--cores", dest='cores', action="store",
                        default=1, type=int,
                        help="cores for multiprocessing workers" +
                        "; default: %(default)s")
    parser.add_argument("-g", "--min_growth", dest='min_growth',
                        action="store",
                        default=0, type=int,
                        help="skip remaining iterations if contig doesnt " +
                        "extend by --min_growth. if 0, ignore" +
                        "; default: %(default)s")
    parser.add_argument("-r", "--reference_genome", dest='reference_genome',
                        action="store", default='', type=str,
                        help="fasta reference genome, used for estimating " +
                        "insert sizes, QUAST, and SPAdes")
    parser.add_argument("-o", "--output", dest='output', action="store",
                        help="output directory; " +
                        "default: %(default)s", default=os.getcwd(), type=str)
    parser.add_argument("--paired_inference", dest='paired_inference',
                        action="store_true", default=False,
                        help="if --paired_inference, mapped read's " +
                        "pairs are included; default: %(default)s")
    parser.add_argument("--subtract", dest='subtract', action="store_true",
                        default=False,
                        help="if --subtract, reads aligned " +
                        "to each reference will not be aligned to future " +
                        "iterations.  Probably you shouldnt do this" +
                        "unless you really happen to want to")
    parser.add_argument("--keep_unmapped", dest='keep_unmapped',
                        action="store_true", default=False,
                        help="if --keep_unmapped fastqs are generated " +
                        "containing the unmapped reads; default: %(default)s")
    parser.add_argument("--ref_as_contig", dest='ref_as_contig',
                        action="store", default="", type=str,
                        help="if 'trusted', SPAdes will  use the seed " +
                        "sequences as a --trusted-contig; if 'untrusted', " +
                        "SPAdes will treat as --untrusted-contig. if '', " +
                        "seeds will not be used during assembly. " +
                        "See SPAdes docs; default: %(default)s")
    parser.add_argument("--temps", dest='temps', action="store_true",
                        default=False,
                        help="if --temps, intermediate files will be " +
                        "kept; default: %(default)s")
    parser.add_argument("-i", "--iterations", dest='iterations',
                        action="store",
                        default=2, type=int,
                        help="if iterations>1, multiple seedings will " +
                        "occur after assembly of seed regions; " +
                        "default: %(default)s")
    parser.add_argument("-v", "--verbosity", dest='verbosity', action="store",
                        default=2, type=int,
                        help="1 = debug(), 2 = info(), 3 = warning(), " +
                        "4 = error() and 5 = critical(); default: %(default)s")
    parser.add_argument("--DEBUG", dest='DEBUG', action="store_true",
                        default=False,
                        help="if --DEBUG, test data will be " +
                        "used; default: %(default)s")
    parser.add_argument("--DEBUG_multi", dest='DEBUG_multiprocessing',
                        action="store_true",
                        default=False,
                        help="if --DEBUG_multiprocessing, runs processes in " +
                        "single loop instead of a multiprocessing pool" +
                        ": %(default)s")
    parser.add_argument("--smalt_scoring", dest='smalt_scoring',
                        action="store",
                        default="match=1,subst=-4,gapopen=-4,gapext=-3",
                        help="submit custom smalt scoring via the smalt -S " +
                        "scorespec option; default: %(default)s")
    ##TODO  Make these check a config file
    parser.add_argument("--spades_exe", dest="spades_exe",
                        action="store", default="spades.py",
                        help="Path to spades executable; default: %(default)s")
    parser.add_argument("--samtools_exe", dest="samtools_exe",
                        action="store", default="samtools",
                        help="Path to bwa executable; default: %(default)s")
    # parser.add_argument("--bcftools_exe", dest="bcftools_exe",
    #                     action="store", default="bcftools",
    #                     help="Path to bwa executable; default: %(default)s")
    parser.add_argument("--smalt_exe", dest="smalt_exe",
                        action="store", default="smalt",
                        help="Path to smalt executable; default: %(default)s")
    parser.add_argument("--quast_exe", dest="quast_exe",
                        action="store", default="quast.py",
                        help="Path to quast executable; default: %(default)s")
    args = parser.parse_args()
    return(args)

# def is_non_zero_file(fpath):
#     """http://stackoverflow.com/questions/2507808/
#        python-how-to-check-file-empty-or-not
#     """
#     return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def check_smalt_full_install(smalt_exe, logger=None):
    smalttestdir = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                "sample_data",
                                "smalt_test","")
    if logger:
        logger.debug("looking for smalt test dir: {0}".format(
            smalttestdir))
    if not os.path.exists(smalttestdir):
        if logger:
            logger.error("cannot find smalt_test dir containing " +\
                         "files to verify bambamc install!")
            sys.exit(1)
    ref = os.path.join(smalttestdir, "ref_to_test_bambamc.fasta")
    index = os.path.join(smalttestdir, "test_index")
    test_bam = os.path.join(smalttestdir, "test_mapping.bam")
    test_reads = os.path.join(smalttestdir, "reads_to_test_bambamc.fastq")
    testindexcmd = str("{0} index {1} {2}".format(smalt_exe, index, ref))
    testmapcmd = str("{0} map -f bam -o {1} {2} {3}".format(smalt_exe,
                                                            test_bam,
                                                            index,
                                                            test_reads))
    if logger:
        logger.debug("testing instalation of smalt and bambamc")
    for i in [testindexcmd, testmapcmd]:
        try:
            if logger:
                logger.debug(i)
            subprocess.run([i],
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        except:
            if logger:
                logger.error("Error running test to check bambamc library is " +\
                             "installed! See https://github.com/gt1/bambamc " +\
                             "and the smalt install guide for more details." +\
                             "https://sourceforge.net/projects/smalt/files/")
            sys.exit(1)
    os.remove(test_bam)
    os.remove(str(index + ".sma"))
    os.remove(str(index + ".smi"))


def map_to_ref_smalt(ref, ref_genome, fastq_read1, fastq_read2,
                     distance_results,
                     map_results_prefix, cores, samtools_exe,
                     smalt_exe, fastq_readS="",
                     read_len=100, step=3, k=5,
                     scoring="match=1,subst=-4,gapopen=-4,gapext=-3"):
    """run smalt based on pased args
    requires at least paired end input, but can handle an additional library
    of singleton reads. Will not work on just singletons
    """
    score_min = int(read_len * .3)
    logger.debug(str("mapping with smalt using a score min of " +
                     "{0}").format(score_min))
    cmdindex = str("{3} index -k {0} -s {1} {2} {2}").format(
        k, step, ref, args.smalt_exe)
    cmdmap = str('{7} map -l pe -S {8} ' +
                 '-m {0} -n {1} -g {2} -f bam -o {3}_pe.bam {4} {5} ' +
                 '{6}').format(score_min, cores, distance_results,
                               map_results_prefix, ref, fastq_read1,
                               fastq_read2, args.smalt_exe, scoring)
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
        #  is there a better. safer way than using mv -f?
        # cmdmerge = str("cp -f {0}_pe.bam " +
        #                "{0}.bam").format(map_results_prefix)
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
                                      samtools_exe=args.samtools_exe)))
    logger.info(str("PE mapped reads: " +
                    get_number_mapped(str(map_results_prefix + "_pe.bam"),
                                      samtools_exe=args.samtools_exe)))
    logger.info(str("Combined mapped reads: " +
                    get_number_mapped(str(map_results_prefix + ".bam"),
                                      samtools_exe=args.samtools_exe)))


def convert_bams_to_fastq(map_results_prefix, fastq_results_prefix,
                          keep_unmapped):
    """ new version with samtools.  Why didnt I do this before?
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
                logger.error(str("No {0}{1}.bam file" +
                                 "found").format(map_results_prefix, i))
                sys.exit(1)
        samfilter = \
            str(args.samtools_exe + " fastq {0}{1}.bam -1 {2}{1}1.fastq -2 " +
                "{2}{1}2.fastq -s {2}{1}S.fastq").format(map_results_prefix, i,
                                                         fastq_results_prefix)
        convert_cmds.append(samfilter)
    logger.debug("running the following commands to extract reads:")
    for i in convert_cmds:
        logger.debug(i)
        subprocess.run(i, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
    if keep_unmapped:
        return("%s%s1.fastq" % (fastq_results_prefix, bams[0]),  # unmapped fwd
               "%s%s2.fastq" % (fastq_results_prefix, bams[0]),  # unmapped rev
               "%s%sS.fastq" % (fastq_results_prefix, bams[0]),  # unmapped s
               "%s%s1.fastq" % (fastq_results_prefix, bams[1]),  # mapped fwd
               "%s%s2.fastq" % (fastq_results_prefix, bams[1]),  # mapped rev
               "%s%sS.fastq" % (fastq_results_prefix, bams[1]))  # mapped s
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
               k="21,33,55,77,99", seqname='', spades_exe="spades.py"):
    """wrapper for common spades setting for long illumina reads
    ref_as_contig should be either blank, 'trusted', or 'untrusted'
    prelim flag is True, only assembly is run, and without coverage correction
    #TODO
    the seqname variable is used only for renaming the resulting contigs
    during iterative assembly.  It would be nice to inheirit from "ref",
    but that is changed with each iteration. This should probably be addressed
    before next major version change
    """
    if groom_contigs not in ['keep_first', 'consensus']:
        logger.error("groom_contigs option must be either keep first or " +
                     "consensus")
        sys.exit(1)
    if seqname == '':
        seqname = ref
    kmers = k  # .split[","]
    success = False
    #  prepare reference, if being used
    if not ref_as_contig == "":
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
        logger.info("Running the following command:\n{0}".format(prelim_cmd))
        subprocess.run(prelim_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        success = output_from_subprocess_exists(os.path.join(output,
                                                             "contigs.fasta"))
        if groom_contigs == "keep_first" and success:
            logger.info("reserving first contig")
            keep_only_first_contig(str(os.path.join(output, "contigs.fasta")),
                                   newname=os.path.splitext(
                                       os.path.basename(seqname))[0])
        elif success and groom_contigs == "consensus":
            # get consensus; copy ref for starters to double check later
            contigs_backup = copy_file(current_file=ref, dest_dir=output,
                                       name=str("backedup_contigs.fasta"),
                                       overwrite=True, logger=logger)

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
            pileup = check_samtools_pileup(os.path.join(output,
                                                        'contigs_pileup.txt'))
            # parse pileup
            consensus = reconstruct_seq(refpath=ref, pileup=pileup,
                                        verbose=False, veryverb=False,
                                        logger=logger)
            with open(os.path.join(output, 'contigs.fasta'), 'w') as new_seqs:
                seqrec = Seq(consensus, IUPAC.IUPACAmbiguousDNA())
                SeqIO.write(SeqRecord(seqrec,
                                      id="contigs_consensus_riboSeed",
                                      description=""), new_seqs, 'fasta')
        else:
            logger.warning("No output from SPAdes this time around")
    else:
        spades_cmd = str(args.spades_exe + " --careful -k {0} {1} {2} -o " +
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
    return(res)


def reconstruct_seq(refpath, pileup, verbose=True, veryverb=False,
                    logger=None):
    """ This is a bit of a mess, to say the least.  Given a list from
    check samtools pileup, and a reference fasta, this reconstructs ambiguous
    regions
    """
    if verbose and logger:
        log_status = logger.info
    elif verbose or veryverb:
        log_status = print
    else:
        pass
    if verbose:
        log_status(str("reconstucting consensus sequence " +
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
                log_status("ref. index: %i\n pile index: %i" % (i, j))
                log_status(ref[i])
                log_status(pileup[j])
            except:
                pass
        # This is how we handle deletions; decrement skip, and next iteration
        if skip > 0:
            skip = skip - 1
            if verbose:
                log_status("skipping {0}".format(i))
        # if j is greater then length of pileup, go with ref.
        # This should avoid out of range issues
        elif j > len(pileup) - 1:
            new = "".join([new, ref[i]])
        #  if index isnt in second col of pileup, skip, filling with ref
        # note because the reference is now zero base, no correction needed
        elif i != int((pileup[j][1])):
            if verbose:
                log_status("no entry in pileup for %i" % i)
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
             (len(pileup[j][4]) == 1 or \
              all(x == pileup[j][4][0] for x in list(pileup[j][4]))) and \
            pileup[j][4][0] not in [",", ".", "^", "$"]:
            new = "".join([new, pileup[j][4][0].upper()])  # append  upper
        # This is tp handle insetions;  could use a lamda?
        elif re.match(insert, pileup[j][4]) is not None and \
            all([hits == re.findall(insert, pileup[j][4])[0] for hits in \
                 re.findall(insert, pileup[j][4])]):
            if verbose:
                log_status("found insert!")
            insert_seq = re.search('[ACGTNacgtn]+', pileup[j][4]).group(0)
            insert_N = int(re.search('[0-9]+', pileup[j][4]).group(0))
            if not len(insert_seq) == insert_N:
                raise ValueError("error parsing insert")
            new = "".join([new, insert_seq])
            indels = indels + insert_N
            N_insertions = N_insertions + insert_N
        # deletions
        elif re.match(delete, pileup[j][4]) is not None and \
            all([hits == re.findall(insert, pileup[j][4])[0] for hits in \
                 re.findall(insert, pileup[j][4])]):
            if verbose:
                log_status("found deletion! {0}".format(pileup[j][4]))
            delete_N = int(re.search('[0-9]+', pileup[j][4]).group(0))
            skip = delete_N
            indels = indels + delete_N
            N_deletions = N_deletions + delete_N
        # Most cases fall in this category
        elif pileup[j][4][0] in [",", ".", "^", "$"] or \
            not all(x == pileup[j][4][0] for x in list(pileup[j][4])):
            if verbose:
                log_status("using ref")
            new = "".join([new, ref[i]])
        else:
            if verbose:
                log_status("Case Not covered!")
            sys.exit(1)
        j = j + 1  # increment the pileup counter
    if verbose:
        log_status(str("total indels: {0}\n\tdeletions {1}\n\tinsetions: " +
                       "{2}").format(indels, N_deletions, N_insertions))
    else:
        print(str("total indels: {0}\n\tdeletions {1}\n\tinsetions: " +
                  "{2}").format(indels, N_deletions, N_insertions))
    return(new[1:])  # [1:] gets rid of starting dollar character


def main(fasta, results_dir, exp_name, mauve_path, map_output_dir, method,
         reference_genome, fastq1, fastq2, fastqS, ave_read_length, cores,
         subtract_reads, ref_as_contig, fetch_mates, keep_unmapped_reads,
         paired_inference, smalt_scoring, min_growth, max_iterations):
    """
    essentially a 'main' function,  to parallelize time comsuming parts
    """
    logger.info("processing {0}".format(fasta))
    logger.info("\nITEM %s of %s\n" % (str(fastas.index(fasta) + 1), nfastas))
    spades_dir = str(results_dir + "SPAdes_" +
                     os.path.split(fasta)[1].split(".fasta")[0])
    quast_dir = str(results_dir + "QUAST_" +
                    os.path.split(fasta)[1].split(".fasta")[0])
    mapping_dir = str(map_output_dir + "mapping_" +
                      os.path.split(fasta)[1].split(".fasta")[0])
    logger.debug(str("this fasta's output dirs: " +
                     "\n{0}\n{1}\n{2}").format(spades_dir, quast_dir,
                                               mapping_dir))
    for i in [spades_dir, quast_dir, mapping_dir]:
        if not os.path.isdir(i):
            os.makedirs(i)
    # if not os.path.isdir(quast_dir):
    #     os.makedirs(quast_dir)
    # if not os.path.isdir(mapping_dir):
    #     os.makedirs(mapping_dir)
    map_results_prefix = os.path.join(
        mapping_dir, str(exp_name + "_" +
                         os.path.split(fasta)[1].split(".fasta")[0]))
    fastq_results_prefix = os.path.join(
        results_dir, str(exp_name + "_" +
                         os.path.split(fasta)[1].split(".fasta")[0]))
    new_reference = copy_file(current_file=fasta, dest_dir=mapping_dir,
                              name='', overwrite=False, logger=logger)
    this_iteration = 1  # counter for iterations.
    # perform iterative mapping
    proceed = True  # this is set to false in while loop if spades fails
    while this_iteration <= max_iterations and proceed:
        # if not the first round, replace ref with extended contigs
        if this_iteration != 1:
            logger.info("copying contigs file for next iteration of assembly")
            new_reference = copy_file(current_file=new_reference,
                                      dest_dir=mapping_dir,
                                      name=str(os.path.basename(fasta) +
                                               "_iter_" + str(this_iteration) +
                                               ".fasta"),
                                      overwrite=True, logger=logger)
        logger.info("Iteration {0} for {1}, item {2} out of {3}".format(
            this_iteration,
            os.path.basename(fasta), fastas.index(fasta) + 1, len(fastas)))
        map_to_ref_smalt(ref=new_reference, ref_genome=reference_genome,
                         fastq_read1=fastq1, fastq_read2=fastq2,
                         fastq_readS=fastqS, read_len=average_read_length,
                         map_results_prefix=map_results_prefix, cores=cores,
                         step=3, k=5,
                         distance_results=os.path.join(results_dir,
                                                       mapped_genome_sam),
                         scoring=smalt_scoring, smalt_exe=args.smalt_exe,
                         samtools_exe=args.samtools_exe)
        extract_mapped_and_mappedmates(map_results_prefix,
                                       fetch_mates=fetch_mates,
                                       samtools_exe=args.samtools_exe,
                                       keep_unmapped=keep_unmapped_reads,
                                       logger=logger)
        logger.info("Converting mapped results to fastqs")
        new_fastq1, new_fastq2, new_fastqS, \
            mapped_fastq1, mapped_fastq2, mapped_fastqS = \
            convert_bams_to_fastq(map_results_prefix, fastq_results_prefix,
                                  keep_unmapped=keep_unmapped_reads,
                                  )
        if subtract_reads and keep_unmapped_reads:
            logger.warning("using reduced reads with next iteration")
            fastq1, fastq2, fastqS = new_fastq1, new_fastq2, new_fastqS
        logger.info("Running SPAdes")
        last_time_through = False
        if this_iteration == args.iterations:
            last_time_through = True
        contigs_path, proceed = \
            run_spades(pe1_1=mapped_fastq1, pe1_2=mapped_fastq2,
                       pe1_s=mapped_fastqS, prelim=True,
                       as_paired=paired_inference,
                       groom_contigs="keep_first",
                       output=spades_dir, keep_best=last_time_through,
                       ref=new_reference, ref_as_contig=ref_as_contig,
                       k="21,33,55",
                       seqname=fasta, spades_exe=args.spades_exe)
        if not proceed:
            logger.warning("Assembly failed: no spades output for {0}".format(
                           os.path.basename(fasta)))
        # compare lengths of reference and freshly assembled contig
        contig_length_diff = get_fasta_lengths(contigs_path)[0] - \
                             get_fasta_lengths(new_reference)[0]
        logger.info("The new contig differs from the reference (or previous " +
                    "iteration) by {0} bases".format(contig_length_diff))
        #  This is a feature that is supposed to help skip unneccesary
        # iterations.  If the difference is negative (new contig is shorter)
        # continue, as this may happen (especially in first mapping if
        # reference is not closely related to Sample), continue to map.
        # If the contig length increases, but not as much as min_growth,
        # skip future iterations
        if contig_length_diff > 0 and contig_length_diff < min_growth and \
           min_growth != 0:  # ie, ignore by default
            this_iteration = max_iterations
        else:
            this_iteration = this_iteration + 1

        # use contigs_path as new reference
        new_reference = contigs_path

    #  Now, after iterative seeding
    try:
        contigs_new_path = copy_file(current_file=contigs_path,
                                     dest_dir=mauve_path,
                                     name=str(os.path.basename(fasta) +
                                              "_final_iter_" +
                                              str(this_iteration) + ".fasta"),
                                     logger=logger)
    except:
        logger.warning("no contigs moved! {0}".format(fasta))
    if not args.temps:
        logger.info("removing temporary files from {0}".format(mapping_dir))
        clean_temp_dir(mapping_dir)
        # clean_temp_dir(mapping_dir, logger=logger)
    if proceed:
        return(0)
    else:
        return(1)


#%%
if __name__ == "__main__":
    args = get_args(DEBUG=False)
    mapped_genome_sam = "genome_distance_est.sam"
    # allow user to give rel path
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root)
    except OSError:
        raise OSError(str("Output directory already exists"))
    map_output_dir = os.path.join(output_root, 'map', "")
    results_dir = os.path.join(output_root, 'results', "")
    mauve_dir = os.path.join(output_root, 'results', "mauve", "")
    t0 = time.time()

    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=str("%s_riboSeed_log.txt" %
                                        os.path.join(output_root,
                                                     time.strftime("%Y%m%d%H%M"))),
                            name=__name__)
    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))
    logger.debug(str("\noutput root {0}\nmap_output_dir: {1}\nresults_dir: " +
                     "{2}\n").format(output_root, map_output_dir, results_dir))
    # check cases of switch-typ args
    if args.ref_as_contig not in ["", 'trusted', 'untrusted']:
        logger.error(str("--ref_as_contig can only be 'trusted', " +
                         "'untrusted', or unused, defaulting to ignoring " +
                         " the reference during assembly. See Spades docs."))
        sys.exit(1)
    if args.method != "smalt":
        logger.error("'smalt' only method currently supported")
        sys.exit(1)
    logger.debug("checking for installations of all required external tools")
    check_installed_tools([args.smalt_exe, args.samtools_exe,
                           args.spades_exe, args.quast_exe], logger=logger)
    # check bambamc is installed proper
    check_smalt_full_install(smalt_exe=args.smalt_exe, logger=logger)
    for i in [map_output_dir, results_dir, mauve_dir]:
        make_outdir(i)
    average_read_length = get_ave_read_len_from_fastq(args.fastq1,
                                                      N=36, logger=logger)
#    map_results_prefix = os.path.join(map_output_dir, args.exp_name)
    fastq_results_prefix = os.path.join(results_dir, args.exp_name)
    #TODO make this with listdirs
    fastas = subprocess.run("ls %s*.fasta" % os.path.join(args.seed_dir, ""),
                            shell=sys.platform != "win32",
                            stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            check=True).stdout.decode().strip().split("\n")
    nfastas = len(fastas)
    logger.debug(fastas)
    ### if using smalt (which you are), check for mapped reference
    if args.method == 'smalt':
        if not os.path.exists(os.path.join(results_dir, mapped_genome_sam)):
            # Index reference for sampling to get PE distances
            logger.info("Estimating insert distances with SMALT")
            # index with default params for genome-sized sequence
            refindex_cmd = str(args.smalt_exe + " index -k {0} -s {1} {2} " +
                               "{2}").format(20, 10, args.reference_genome)
            refsample_cmd = \
                str(args.smalt_exe + " sample -n {0} -o {1} {2} {3} " +
                    "{4}").format(args.cores,
                                  os.path.join(results_dir, mapped_genome_sam),
                                  args.reference_genome,
                                  args.fastq1, args.fastq2)
            logger.info("Sampling and indexing {0}".format(
                args.reference_genome))
            for cmd in [refindex_cmd, refsample_cmd]:
                logger.debug("\t command:\n\t {0}".format(cmd))
                subprocess.run(cmd,
                               shell=sys.platform != "win32",
                               stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               check=True)
        else:
            logger.info("using existing reference file")
    else:
        logger.error("As of v 0.88, only supported mapper is 'smalt'")
        sys.exit(1)
    ### Main function call
    if args.DEBUG_multiprocessing:
        logger.warning("running without multiprocessing!")
        for i in fastas:
            main(fasta=i, results_dir=results_dir, exp_name=args.exp_name,
                 mauve_path=mauve_dir, map_output_dir=map_output_dir,
                 method=args.method, reference_genome=args.reference_genome,
                 fastq1=args.fastq1, fastq2=args.fastq2, fastqS=args.fastqS,
                 ave_read_length=average_read_length, cores=args.cores,
                 subtract_reads=args.subtract,
                 ref_as_contig=args.ref_as_contig,
                 fetch_mates=args.paired_inference,
                 keep_unmapped_reads=args.keep_unmapped,
                 paired_inference=args.paired_inference,
                 smalt_scoring=args.smalt_scoring,
                 min_growth=args.min_growth,
                 max_iterations=args.iterations)
    else:
        pool = multiprocessing.Pool(processes=args.cores)
        # cores_per_process =
        results = [pool.apply_async(main, (fasta,),
                                    {"results_dir": results_dir,
                                     "map_output_dir": map_output_dir,
                                     "exp_name": args.exp_name,
                                     "method": args.method,
                                     "reference_genome": args.reference_genome,
                                     "fastq1": args.fastq1,
                                     "fastq2": args.fastq2,
                                     "fastqS": args.fastqS,
                                     "cores": 4,  # cores": args.cores,
                                     "mauve_path": mauve_dir,
                                     "ave_read_length": average_read_length,
                                     "fetch_mates": args.paired_inference,
                                     "keep_unmapped_reads": args.keep_unmapped,
                                     "paired_inference": args.paired_inference,
                                     "subtract_reads": args.subtract,
                                     "ref_as_contig": args.ref_as_contig,
                                     "smalt_scoring": args.smalt_scoring,
                                     "min_growth": args.min_growth,
                                     "max_iterations": args.iterations})
                   for fasta in fastas]
        pool.close()
        pool.join()
        logger.info(results)
        logger.info(sum([r.get() for r in results]))
    logging.info("combinging contigs from %s" % mauve_dir)
    new_contig_file = combine_contigs(contigs_dir=mauve_dir,
                                      contigs_name="riboSeedContigs",
                                      logger=logger)
    logger.info("Combined Seed Contigs: {0}".format(new_contig_file))
    logger.info("Time taken to run seeding: %.2fs" % (time.time() - t0))
    logger.info("\n\n Starting Final Assemblies\n\n")
    for j in "de_novo", "de_fere_novo":
        logging.info("\n\nRunning %s SPAdes \n" % j)
        if j == "de_novo":
            assembly_ref = ''
            assembly_ref_as_contig = ''
        elif j == "de_fere_novo":
            assembly_ref = new_contig_file
            assembly_ref_as_contig = 'trusted'
        else:
            raise ValueError("onl valid cases are de novo and de fere novo")
        output_contigs, \
            final_success = run_spades(pe1_1=args.fastq1, pe1_2=args.fastq2,
                                       output=os.path.join(results_dir, j),
                                       ref=assembly_ref,
                                       ref_as_contig=assembly_ref_as_contig,
                                       prelim=False, keep_best=False,
                                       k="21,33,55,77,99,127")
        if final_success:
            logger.info("\n\nRunning %s QUAST" % j )
            run_quast(contigs=output_contigs,
                      output=os.path.join(results_dir, str("quast_" + j)),
                      quast_exe=args.quast_exe,
                      ref=args.reference_genome,
                      logger=logger)
    # Report that we've finished
    logger.info("Done: %s." % time.asctime())
    logger.info("Time taken: %.2fs" % (time.time() - t0))

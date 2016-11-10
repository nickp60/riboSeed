#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
Minor Version Revisions:
 - no more individual versions; all treated as pipeline from here on
starting at version 0.0.940
Created on Sun Jul 24 19:33:37 2016

See README.md for more info and usage


### iooutline
read in (multi)gb, and a riboSelect output
constroct a seedGenome
- append seq record
- append all needed loci_clusters
- append all needed Locus

map all reads to genome
extract and link files to loci

extract



###



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


from riboSeed.riboSnag import parse_clustered_loci_file, \
    extract_coords_from_locus, \
    stitch_together_target_regions, get_genbank_rec_from_multigb, \
    pad_genbank_sequence, prepare_prank_cmd, prepare_mafft_cmd, \
    calc_Shannon_entropy, calc_entropy_msa, \
    annotate_msa_conensus, plot_scatter_with_anno, \
    profile_kmer_occurances, plot_pairwise_least_squares, make_msa, \
    LociCluster, Locus

## GLOBALS
SAMTOOLS_MIN_VERSION = '1.3.1'
#################################### classes ###############################


class SeedGenome(object):
    """ organizes the clustering process instead of dealing with individual
    seeds alone
    This holds all he data pertaining to te clustering of a scaffold
    """
    def __init__(self, genbank_path, name=None, fasta_path=None,
                 loci_clusters=None, output_root=None,
                 riboSelect_path=None, seq_records=None,
                 initial_map_prefix=None, ngs_ob=None, logger=None):
        self.name = name  # get from commsanline in case running multiple
        self.output_root = output_root
        self.genbank_path = genbank_path
        self.fasta_path = fasta_path  # this is made dynamically
        self.riboSelect_path = riboSelect_path
        self.loci_clusters = loci_clusters  # this holds LociCluster objects
        self.seq_records = seq_records  # this is set
        self.initial_map_prefix = initial_map_prefix  # set this dynamically
        self.ngs_ob = ngs_ob  # for ngslib object
        self.logger = logger
        self.write_fasta_genome()
        self.attach_genome_seqRecord()
        self.check_records()
        self.make_map_prefix()

    def make_map_prefix(self):
        initial_map_dir = os.path.join(self.output_root,
                                       str(self.name + "_initial_mapping"))
        if not os.path.isdir(initial_map_dir):
            os.makedirs(initial_map_dir)
        self.initial_map_prefix = os.path.join(
            self.output_root,
            str(self.name + "_initial_mapping"),
            "initial_mapping")

    def write_fasta_genome(self):
        """
        """
        self.name = os.path.splitext(
            os.path.basename(self.genbank_path))[0]
        self.fasta_path = os.path.join(self.output_root,
                                       str(self.name + ".fasta"))
        with open(self.genbank_path, 'r') as fh:
            with open(self.fasta_path, 'w') as outfh:
                sequences = SeqIO.parse(fh, "genbank")
                count = SeqIO.write(sequences, outfh, "fasta")
                # print("re-wrote %i sequences as fasta" % count)

    def attach_genome_seqRecord(self):
        """
        """
        with open(self.genbank_path, 'r') as fh:
            self.seq_records = list(SeqIO.parse(fh, "genbank"))

    def check_records(self):
        assert len(list(SeqIO.parse(self.fasta_path, "fasta"))) == \
            len(self.seq_records), "Error parsing genbank file!"


class ngsLib(object):
    """paired end data object
    """
    def __init__(self, name, readF=None, readR=None,
                 readS0=None, readS1=None,
                 # readS1=None, readS2=None, readS3=None,
                 smalt_dist_path=None, readlen=None,
                 libtype=None, logger=None, smalt_exe=None,
                 ref_fasta=None):
        self.name = name
        self.libtype = libtype  # set this dynamically
        self.readF = readF
        self.readR = readR
        self.readS0 = readS0
        self.readS1 = readS1
        self.readlen = readlen  # set this dynamically
        self.logger = logger
        self.smalt_exe = smalt_exe
        self.ref_fasta = ref_fasta
        self.smalt_dist_path = smalt_dist_path  # set this dynamically
        self.set_libtype()
        self.get_readlen()
        self.smalt_insert_file()

    def set_libtype(self):
        """sets to either s_1, pe, pe_s
        TODO: add mate support and support for multiple single libraries
        """
        if self.readF is None:
            if self.readR is None:
                if self.readS0 is not None:
                    self.libtype = "s_1"
                else:
                    raise ValueError("cannot set library type")
            else:
                raise ValueError("cannot set library type from one PE read")
        else:
            if self.readS0 is not None:
                self.libtype = "pe_s"
            else:
                self.libtype = "pe"

    def get_readlen(self):
        if self.libtype in ['pe', 'pe_s']:
            self.readlen = get_ave_read_len_from_fastq(
                self.readF, N=36, logger=self.logger)
        else:
            self.readlen = get_ave_read_len_from_fastq(
                self.readS0, N=36, logger=self.logger)

    def smalt_insert_file(self):
        if self.libtype in ['pe', 'pe_s']:
            self.smalt_dist_path = estimate_distances_smalt(
                outfile=os.path.join(os.path.dirname(self.ref_fasta),
                                     "smalt_distance_est.sam"),
                smalt_exe=self.smalt_exe,
                ref_genome=self.ref_fasta,
                fastq1=self.readF, fastq2=self.readR,
                logger=self.logger)
        else:
            return None


# class LociCluster(object):
#     """ organizes the clustering process instead of dealing with nested lists
#     This holds the whole cluster of one to several individual loci
#     """
#     def __init__(self, index, sequence_id, loci_list, padding=None,
#                  global_start_coord=None, global_end_coord=None,
#                  seq_record=None, feat_of_interest=None, mappings=None,
#                  extractedSeqRecord=None, cluster_dir_name=None,
#                  circular=False):
#         self.index = index
#         self.sequence_id = sequence_id
#         self.loci_list = loci_list  # this holds the Locus objects
#         self.global_start_coord = global_start_coord
#         self.global_end_coord = global_end_coord
#         self.padding = padding
#         self.circular = circular
#         self.cluster_dir_name = cluster_dir_name  # named dynamically
#         self.mappings = mappings
#         self.feat_of_interest = feat_of_interest
#         self.extractedSeqRecord = extractedSeqRecord
#         self.name_cluster_dir()

#     def name_mapping_dir(self):
#         self.cluster_dir_name = str("{0}_cluster_{1}").format(
#             self.sequence_id, self.index)


class LociMapping(object):
    """ order of operations: map to reference, extract and convert,
    assemble, save results here
    """
    def __init__(self, iteration, mapping_subdir,
                 mapping_success=False, assembly_success=False,
                 ref_path=None, pe_map_bam=None, s_map_bam=None,
                 sorted_map_bam=None, merge_map_sam=None, mapped_bam=None,
                 merge_map_bam=None, mapped_sam=None, spades_subdir=None,
                 unmapped_sam=None, mappedF=None, mappedR=None,
                 mapped_ids_txt=None, unmapped_ids_txt=None, unmapped_bam=None,
                 mappedS=None, assembled_contig=None, assembly_subdir=None,
                 unmappedF=None, unmappedR=None, unmappedS=None):
        self.iteration = iteration
        self.mapping_success = mapping_success
        self.assembly_success = assembly_success
        self.mapping_subdir = mapping_subdir
        self.assembly_subdir = assembly_subdir
        self.ref_path = ref_path
        self.pe_map_bam = pe_map_bam  # all reads from pe mapping
        self.s_map_bam = s_map_bam  # all reads from singltons mapping
        self.merge_map_bam = merge_map_bam  # combined pe and s as bam
        self.merge_map_sam = merge_map_sam  # combined pe and s as sam
        self.mapped_sam = mapped_sam  # mapped reads only, sam
        self.mapped_bam = mapped_bam  # mapped reads only, bam
        self.mapped_ids_txt = mapped_ids_txt
        self.unmapped_ids_txt = unmapped_ids_txt
        self.unmapped_sam = unmapped_sam
        self.unmapped_bam = unmapped_bam
        #  The fastqs
        self.mappedF = mappedF
        self.mappedR = mappedR
        self.mappedS = mappedS
        self.unmappedF = unmappedF
        self.unmappedR = unmappedR
        self.unmappedS = unmappedS
        self.sorted_map_bam = sorted_map_bam   # used with intial mapping
        # self.sorted_map_subset_bam = sorted_map_subset_bam
        self.assembled_contig = assembled_contig
        ###
        self.make_mapping_subdir()
        self.make_assembly_subdir()
        # self.name_id_file()

    def make_mapping_subdir(self):
        if not os.path.isdir(self.mapping_subdir):
            os.makedirs(self.mapping_subdir)

    def make_assembly_subdir(self):
        if self.assembly_subdir is not None:
            if not os.path.isdir(self.assembly_subdir):
                os.makedirs(self.assembly_subdir)


# class Locus(object):
#     """ this holds the info for each individual Locus"
#     """
#     def __init__(self, index, sequence_id, locus_tag, strand=None,
#                  start_coord=None, end_coord=None, rel_start_coord=None,
#                  rel_end_coord=None, product=None):
#         # self.parent ??
#         self.index = index
#         self.sequence_id = sequence_id  # is this needed? I dont think so as long
#         # as a locus is never decoupled from the LociCluster
#         self.locus_tag = locus_tag
#         self.strand = strand  # 1 is +, -1 is -
#         self.start_coord = start_coord
#         self.end_coord = end_coord
#         # self.rel_start_coord = rel_start_coord  # start relative to length of region
#         # self.rel_end_coord = rel_end_coord  # end relative to length of region
#         self.product = product

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
    optional.add_argument("-l", "--flanking_length",
                          help="length of flanking regions, can be colon-" +
                          "separated to give separate upstream and " +
                          "downstream flanking regions; default: %(default)s",
                          default='1000', type=str, dest="flanking")
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

########################  funtions adapted from elsewhere ###################


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


def estimate_distances_smalt(outfile, smalt_exe, ref_genome,
                             fastq1, fastq2, cores=None, logger=None):
    """Given fastq pair and a reference, returns path to distance estimations
    used by smalt to help later with mapping. if one already exists,
    return path to it.
    """
    if cores is None:
        cores = multiprocessing.cpu_count()
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


def map_to_genome_smalt(
        seed_genome, ngsLib, map_results_prefix, cores,
        samtools_exe, smalt_exe, score_minimum=None, step=3,
        k=5, scoring="match=1,subst=-4,gapopen=-4,gapext=-3", logger=None):
    """run smalt based on pased args
    requires at least paired end input, but can handle an additional library
    of singleton reads. Will not work on just singletons
    """
    # check min score
    if score_minimum is None:
        score_min = int(ngsLib.readlen * .3)
    else:
        score_min = score_minimum
    logger.debug(str("mapping with smalt using a score min of " +
                     "{0}").format(score_min))
    cmdindex = str("{3} index -k {0} -s {1} {2} {2}").format(
        k, step, seed_genome.fasta_path, smalt_exe)
    cmdmap = str('{7} map -l pe -S {8} ' +
                 '-m {0} -n {1} -g {2} -f bam -o {3}_pe.bam {4} {5} ' +
                 '{6}').format(score_min, cores, ngsLib.smalt_dist_path,
                               map_results_prefix, seed_genome.fasta_path,
                               ngsLib.readF,
                               ngsLib.readR, smalt_exe, scoring)
    smaltcommands = [cmdindex, cmdmap]

    if ngsLib.readS0 is not None:
        cmdindexS = str('{0} index -k {1} -s {2} {3} {3}').format(
            smalt_exe, k, step, seed_genome.fasta_path)
        cmdmapS = str("{7} map -S {6} " +
                      "-m {0} -n {1} -g {2} -f bam -o {3}S.bam {4} " +
                      "{5}").format(score_min, cores, ngsLib.smalt_dist_path,
                                    map_results_prefix, seed_genome.fasta_path,
                                    ngsLib.readS0,
                                    scoring, smalt_exe)
        cmdmergeS = str('{0} merge -f  {1}.bam {1}_pe.bam ' +
                        '{1}S.bam').format(samtools_exe, map_results_prefix)
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
    if ngsLib.readS0 is not None:
        logger.info(str("Singleton mapped reads: " +
                        get_number_mapped(str(map_results_prefix + "S.bam"),
                                          samtools_exe=samtools_exe)))
    logger.info(str("PE mapped reads: " +
                    get_number_mapped(str(map_results_prefix + "_pe.bam"),
                                      samtools_exe=samtools_exe)))
    logger.info(str("Combined mapped reads: " +
                    get_number_mapped(str(map_results_prefix + ".bam"),
                                      samtools_exe=samtools_exe)))


def convert_bams_to_fastq(mapping_ob,
                          keep_unmapped, samtools_exe, logger=None):
    """ return 6 paths: unmapped F, R, S and mapped F, R. S
    given the prefix for the mapped bam files, write out mapped (and optionally
    unmapped) reads to fasta files
    """
    output_paths = {'mapped_bam': ['mappedF', 'mappedR', 'mappedS'],
                    'unmapped_bam': ['unmappedF', 'unmappedR', 'unmappedS']}
    for key, value in output_paths.items():
        for fastq in value:
            setattr(mapping_ob, fastq,
                    str(os.path.splitext(mapping_ob.merge_map_bam)[0] +
                        "_" + fastq + '.fastq'))
    convert_cmds = []
    for key, val in output_paths.items():
        if not os.path.exists(getattr(mapping_ob, key)):
            if key == 'unmapped_bam' and not keep_unmapped:
                continue
            else:
                if logger:
                    logger.error(str("No {0} file found").format(
                        getattr(mapping_ob, key)))
                raise FileNotFoundError("No {0} file found".format(
                    getattr(mapping_ob, key)))
        samfilter = "{0} fastq {1} -1 {2} -2 {3} -s {4}".format(
            samtools_exe, getattr(mapping_ob, key),
            getattr(mapping_ob, val[0]),
            getattr(mapping_ob, val[1]),
            getattr(mapping_ob, val[2]))
        convert_cmds.append(samfilter)
    if logger:
        logger.debug("running the following commands to extract reads:")
    for i in convert_cmds:
        if logger:
            logger.debug(i)
        subprocess.run(i, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)



def run_spades(
        mapping_ob, ref_as_contig, as_paired=True, keep_best=True,
        prelim=False, groom_contigs='keep_first', k="21,33,55,77,99",
        seqname='', spades_exe="spades.py", logger=None):
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
        seqname = mapping_ob.ref_path
    kmers = k  # .split[","]
    #  prepare reference, if being used
    if not ref_as_contig is None:
        alt_contig = "--{0}-contigs {1}".format(
            ref_as_contig, mapping_ob.ref_path)
    else:
        alt_contig = ''
    # prepare read types, etc
    if as_paired and mapping_ob.mappedS is not None:  # for lib with both
        singles = "--pe1-s {0} ".format(mapping_ob.mappedS)
        pairs = "--pe1-1 {0} --pe1-2 {1} ".format(
            mapping_ob.mappedF, mapping_ob.mappedR)
    elif as_paired and mapping_ob.mappedS is None:  # for lib with just PE
        singles = ""
        pairs = "--pe1-1 {0} --pe1-2 {1} ".format(
            mapping_ob.mappedF, mapping_ob.mappedR)
    # for libraries treating paired ends as two single-end libs
    elif not as_paired and mapping_ob.mappedS is None:
        singles = ''
        pairs = "--pe1-s {0} --pe2-s {1} ".format(
            mapping_ob.mappedF, mapping_ob.mappedR)
    else:  # for 3 single end libraries
        singles = "--pe1-s {0} ".format(mapping_ob.mappedS)
        pairs = str("--pe2-s {0} --pe3-s {1} ".format(
            mapping_ob.mappedF, mapping_ob.mappedR))
    reads = str(pairs + singles)
#    spades_cmds=[]
    if prelim:
        prelim_cmd = str(
            "{0} --only-assembler --cov-cutoff off --sc --careful -k {1} " +
            "{2} {3} -o {4}").format(spades_exe, kmers, reads, alt_contig,
                                     mapping_ob.assembly_subdir)
        logger.debug("Running SPAdes command:\n%s", prelim_cmd)
        subprocess.run(prelim_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        mapping_ob.assembly_success = output_from_subprocess_exists(
            os.path.join(mapping_ob.assembly_subdir, "contigs.fasta"))
        if groom_contigs == "keep_first" and mapping_ob.assembly_success:
            logger.info("reserving first contig")
            try:
                keep_only_first_contig(
                    os.path.join(mapping_ob.assembly_subdir, "contigs.fasta"),
                    newname=os.path.splitext(os.path.basename(seqname))[0])
            except Exception as f:
                logger.error(f)
                raise f
        else:
            logger.warning("No output from SPAdes this time around")
    else:
        spades_cmd = "{0} --careful -k {1} {2} {3} -o {4}".format(
            spades_exe, kmers, reads, alt_contig, mapping_ob.assembly_subdir)
        logger.warning("Running the following command:\n{0}".format(spades_cmd))
        subprocess.run(spades_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
        # not check=True; dont know spades return codes
        mapping_ob.assembly_success = output_from_subprocess_exists(
            os.path.join(mapping_ob.assembly_subdir, "contigs.fasta"))


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
            raise e

    return mainDict


def partition_mapped_reads(seedGenome, samtools_exe,
                           flank=[0, 0], logger=None):
    """ This handles the first round of mapped reads
    """
    # make a list of regions mapping
    mapped_regions = []
    for cluster in seedGenome.loci_clusters:
        """ make iteration 0 LociMapping object
        """
        if logger is None:
            raise ValueError("Must have logger for this function")

        ## make maping objec
        mapping_subdir = os.path.join(
            seedGenome.output_root,
            cluster.cluster_dir_name,
            "{0}_cluster_{1}_mapping_iteration_0".format(
                cluster.sequence_id, cluster.index))
        assembly_subdir = os.path.join(
            seedGenome.output_root,
            cluster.cluster_dir_name,
            "{0}_cluster_{1}_assembly_iteration_0".format(
                cluster.sequence_id, cluster.index))
        merge_map_bam = str("{0}{1}{2}_{3}_{4}").format( #merged_map bam so it integrates
            mapping_subdir, os.path.sep, "cluster", cluster.index,
            "sorted_subset.bam")
        sorted_map_bam = str(seedGenome.initial_map_prefix +
                             "_sorted.bam")

        mapping0 = LociMapping(iteration=0,
                               mapping_subdir=mapping_subdir,
                               assembly_subdir=assembly_subdir,
                               # sorted_map_subset_bam=sorted_map_subset_bam,
                               sorted_map_bam=sorted_map_bam,
                               ref_path=None,
                               pe_map_bam=None,
                               s_map_bam=None,
                               merge_map_bam=merge_map_bam,
                               mapped_sam=None,
                               unmapped_sam=None, mappedF=None,
                               mappedR=None,
                               mappedS=None)
        if sorted([x.start_coord for x in cluster.loci_list]) != \
           [x.start_coord for x in cluster.loci_list]:
            logger.warning("Coords are not in increasing order; " +
                           "you've been warned")
        start_list = sorted([x.start_coord for x in cluster.loci_list])
        logger.debug("Start_list: {0}".format(start_list))

        logger.debug("Find coordinates to gather reads from the following coords:")
        for i in cluster.loci_list:
            logger.debug(str(i.__dict__))
        #  This works as long as coords are never in reverse order
        cluster.global_start_coord = min([x.start_coord for
                                          x in cluster.loci_list]) - flank[0]
        # if start is negative, just use 1, the beginning of the sequence
        if cluster.global_start_coord < 1:
            logger.warning(
                "Caution! Cannot retrieve full flanking region, as " +
                "the 5' flanking region extends past start of " +
                "sequence. If this is a problem, try using a smaller " +
                "--flanking region, and/or if  appropriate, run with " +
                "--circular.")
            cluster.global_start_coord = 1
        cluster.global_end_coord = max([x.end_coord for
                                        x in cluster.loci_list]) + flank[1]
        if cluster.global_end_coord > len(cluster.seq_record):
            logger.warning(
                "Caution! Cannot retrieve full flanking region, as " +
                "the 5' flanking region extends past start of " +
                "sequence. If this is a problem, try using a smaller " +
                "--flanking region, and/or if  appropriate, run with " +
                "--circular.")
            cluster.global_end_coord = len(cluster.seq_record)
        logger.debug("global start and end: %s %s", cluster.global_start_coord,
                     cluster.global_end_coord)

        logger.warning("Extracting the sequence: %s %s",
                       cluster.global_start_coord,
                       cluster.global_end_coord)
        cluster.extractedSeqRecord = SeqRecord(
            cluster.seq_record.seq[
                cluster.global_start_coord:
                cluster.global_end_coord])
        mapping0.ref_path = os.path.join(mapping0.mapping_subdir,
                                         "extracted_seed_sequence.fasta")
        with open(mapping0.ref_path, "w") as writepath:
            SeqIO.write(cluster.extractedSeqRecord, writepath, 'fasta')

        # Prepare for partitioning
        partition_cmds = []
        if not os.path.exists(mapping0.sorted_map_bam):
            sort_cmd = str("{0} sort {1} > {2}").format(
                samtools_exe, str(seedGenome.initial_map_prefix + ".bam"),
                mapping0.sorted_map_bam)
            index_cmd = str("{0} index {1}").format(
                samtools_exe, mapping0.sorted_map_bam)
            partition_cmds.extend([sort_cmd, index_cmd])
        region_to_extract = "{0}:{1}-{2}".format(
            cluster.sequence_id, cluster.global_start_coord,
            cluster.global_end_coord)
        view_cmd = str("{0} view -o {1} {2} {3}").format(
            samtools_exe, mapping0.merge_map_bam,
            mapping0.sorted_map_bam,
            region_to_extract)
        partition_cmds.append(view_cmd)
        mapped_regions.append(region_to_extract)
        ### run cmds
        cluster.mappings.append(mapping0)
        for cmd in partition_cmds:
            subprocess.run([cmd],
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
    # get unmapped poo
    init_unmapped = LociMapping(iteration=0,
                                mapping_subdir=os.path.join(
                                    seedGenome.output_root, str(
                                        seedGenome.name +
                                        "_unmapped_iteration_0")),
                                ref_path=None,
                                pe_map_bam=None,
                                s_map_bam=None,
                                sorted_map_bam=str(
                                    seedGenome.initial_map_prefix +
                                    "_sorted.bam"),
                                merge_map_bam=None, mapped_sam=None,
                                unmapped_sam=None, mappedF=None,
                                mappedR=None,
                                mappedS=None)
    init_unmapped.merge_map_bam = str(
        init_unmapped.mapping_subdir +
        os.path.sep +
        "unmapped_subset.bam")

    unmapped_view_cmd = str("{0} view -o {1} {2} -U {3}").format(
        samtools_exe, init_unmapped.merge_map_bam,
        init_unmapped.sorted_map_bam,
        ' '.join([x for x in mapped_regions]))
    subprocess.run([unmapped_view_cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def add_coords_to_clusters(seedGenome, logger=None):
    for cluster in seedGenome.loci_clusters:  # for each cluster of loci
        # get seq record that cluster is  from
        try:
            cluster.seq_record = \
                get_genbank_rec_from_multigb(
                    recordID=cluster.sequence_id,
                    genbank_records=seedGenome.seq_records)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        # make coord list
        try:
            extract_coords_from_locus(
                cluster=cluster, feature=cluster.feat_of_interest,
                logger=logger)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        logger.info(str(cluster.__dict__))


def main(fasta, fastas, num, results_dir, exp_name, mauve_path, map_output_dir, method,
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

    fastq_results_prefix = os.path.join(results_dir, args.exp_name)

    ### #new logic starts here
    # make seedGenome object
    seedGenome = SeedGenome(
        riboSelect_path=args.riboSelect_path,
        output_root=output_root,
        genbank_path=args.genbank_genome,
        logger=logger)

    ### add ngsobject
    seedGenome.ngs_ob = ngsLib(name="",
                               readF=args.fastq1,
                               readR=args.fastq2,
                               readS0=args.fastqS,
                               logger=logger,
                               smalt_exe=args.smalt_exe,
                               ref_fasta=seedGenome.fasta_path)

    # read in riboSelect clusters
    seedGenome.loci_clusters = parse_clustered_loci_file(
        filepath=seedGenome.riboSelect_path,
        gb_filepath=seedGenome.genbank_path,
        padding=args.padding,
        circular=args.circular,
        logger=logger)

    # add coordinates
    add_coords_to_clusters(seedGenome=seedGenome,
                           logger=logger)

    # Run commands to map to the genome
    map_to_genome_smalt(seed_genome=seedGenome,
                        ngsLib=seedGenome.ngs_ob,
                        map_results_prefix=seedGenome.initial_map_prefix,
                        cores=args.cores,
                        samtools_exe=args.samtools_exe,
                        smalt_exe=args.smalt_exe,
                        score_minimum=None,
                        step=3, k=5,
                        scoring="match=1,subst=-4,gapopen=-4,gapext=-3",
                        logger=logger)
    partition_mapped_reads(
        seedGenome=seedGenome,
        samtools_exe=args.samtools_exe,
        flank=[0, 0],
        logger=logger)

    # now, we need to assemble each mapping object
    if args.DEBUG_multiprocessing:
        logger.warning("running without multiprocessing!")
        for cluster in seedGenome.loci_clusters:
            assemble_initial_mapping(
                cluster,
                nseqs=len(seedGenome.loci_clusters),
                fetch_mates=False,
                target_len=6000,
                samtools_exe=args.samtools_exe,
                keep_unmapped_reads=False,
                logger=logger)
    else:
        pool = multiprocessing.Pool(processes=args.cores)
        results = [pool.apply_async(assemble_initial_mapping,
                                    (cluster,),
                                    {"nseqs": len(seedGenome.loci_clusters),
                                     "fetch_mates": False,
                                     "target_len": 6000,
                                     "keep_unmapped_reads": False,
                                     "samtools_exe": args.samtools_exe,
                                     "logger": logger})
                   for cluster in seedGenome.loci_clusters]


        # logger.info("%s Converting mapped results to fastqs", prelog)
        # try:
        #     new_fastq1, new_fastq2, new_fastqS, \
        #         mapped_fastq1, mapped_fastq2, mapped_fastqS = \
        #             convert_bams_to_fastq(map_results_prefix,
        #                                   fastq_results_prefix,
        #                                   keep_unmapped=keep_unmapped_reads,
        #                                   logger=logger,
        #                                   samtools_exe=args.samtools_exe)
        # except Exception as e:
        #     logger.error(e)
        #     sys.exit(1)

def assemble_initial_mapping(clu, nseqs, target_len, fetch_mates, logger,
                             samtools_exe, keep_unmapped_reads=False,
                             prelim=True):
    prelog = "{0}-{1}:".format("SEED_cluster", clu.index)
    logger.info("%s processing initial mapping", prelog)
    logger.info("%s item %i of %i", prelog, clu.index + 1, nseqs)
    logger.debug("%s output dirs: \n%s\n%s", prelog,
                 clu.mappings[0].assembly_subdir,
                 clu.mappings[0].mapping_subdir)
    try:
        extract_mapped_reads(mapping_ob=clu.mappings[0],
                             fetch_mates=fetch_mates,
                             samtools_exe=samtools_exe,
                             keep_unmapped=keep_unmapped_reads,
                             logger=logger)
    except Exception as e:
        logger.error(e)
        sys.exit(1)
    try:
        convert_bams_to_fastq(mapping_ob=clu.mappings[0],
                              keep_unmapped=False,
                              samtools_exe=samtools_exe,
                              logger=logger)

    except Exception as e:
        logger.error(e)
        sys.exit(1)
    logger.info("%s Running SPAdes", prelog)
    try:
        run_spades(
            mapping_ob=clu.mappings[0], ref_as_contig='trusted',
            as_paired=False, keep_best=True, prelim=prelim,
            groom_contigs='keep_first', k="21,33,55,77,99",
            seqname='', spades_exe="spades.py", logger=logger)

    except Exception as e:
        logger.error("SPAdes error:")
        logger.error(e)
        sys.exit(1)
    # if not clu.mapping.assembly_success:
    #     logger.warning("%s Assembly failed: no spades output for %s",
    #                    prelog, os.path.basename(fasta))
    return 0


def extract_mapped_reads(mapping_ob, fetch_mates,
                         keep_unmapped, samtools_exe, logger=None):
    """
    Take a prefix for a dir containing your mapped bam file.
    IF fetch_mates is true, mapped reads are extracted,
    and mates are feteched with the LC_ALL line.If not, that part is
    skipped, and just the mapped reads are extracted.
    Setting keep_unmapped to true will output a bam file with
    all the remaining reads.  This could be used if you are really confident
    there are no duplicate mapping_obs you are interested in.
     -F 4 option selects mapped reads
    Note that the umapped output includes reads whose pairs were mapped.
    This is to try to catch the stragglers.
    LC_ALL=C  call from pierre lindenbaum. No idea how it does, but its magic
    """
    all_files = {'sam': ['merge_map_sam', 'unmapped_sam', 'mapped_sam'],
                 'txt': ['unmapped_ids_txt', 'mapped_ids_txt'],
                 'bam': ['unmapped_bam', 'mapped_bam', 'unmapped_bam']}
    extract_cmds = []
    for key, values in all_files.items():
        for value in values:
            setattr(mapping_ob, value,
                    str(os.path.splitext(mapping_ob.merge_map_bam)[0] +
                        "_" + value + "." + key))

    # Either get nates or ignore mates
    if fetch_mates:
        makesam = "{0} view -o {1}".format(samtools_exe, mapping_ob.merge_map_sam)
        samview = str("{0} view -h -F 4 {1} | cut -f1 > {2}").format(
            samtools_exe, mapping_ob.merge_map_bam, mapping_ob.mapped_ids_txt)
        lc_cmd = str("LC_ALL=C grep -w -F -f {0}  < {1} > {2}").format(
            mapping_ob.mapped_ids_txt, mapping_ob.merge_map_sam, mapping_ob.mapped_sam)
        extract_cmds.extend([makesam, samview, lc_cmd])
    else:
        samview = str("{0} view -hS -F 4 {1} > {2}").format(
            samtools_exe, mapping_ob.merge_map_bam, mapping_ob.mapped_sam)
        extract_cmds.extend([samview])
    samsort = str("{0} view -bhS {1} | samtools sort - > {2}").format(
        samtools_exe, mapping_ob.mapped_sam, mapping_ob.mapped_bam)
    samindex = " {0} index {1}".format(samtools_exe, mapping_ob.mapped_bam)
    extract_cmds.extend([samsort, samindex])
    if keep_unmapped:
        samviewU = str("{0} view -f 4 {1} | cut -f1 > {2}").format(
            samtools_exe, mapping_ob.mapped_bam, mapping_ob.unmapped_ids_txt)
        lc_cmdU = str("LC_ALL=C grep -w -F -f {0} < {1} > {2}").format(
            mapping_ob.unmapped_ids_txt, mapping_ob.merge_map_sam,
            mapping_ob.unmapped_sam)
        samindexU = str("{0} view -bhS {1} " +
                        "| {0} sort - -o {2} && {0} index {2}").format(
            samtools_exe, mapping_ob.unmapped_sam, mapping_ob.unmapped_bam,)
        extract_cmds.extend([samviewU, lc_cmdU, samindexU])
    if logger:
        logger.debug("running the following commands to extract reads:")
    for i in extract_cmds:
        if logger:
            logger.debug(i)
        subprocess.run(i, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
    return 0



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
            final_success = run_spades(
                pe1_1=args.fastq1, pe1_2=args.fastq2,
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
        else:
            logger.error("Some error occured during final assemblies; " +
                         "SPAdes logs")
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
            logger.error("Error writing out combined quast report")
            logger.error(e)
        logger.info("Comparing de novo and de fere novo assemblies:")
    # Report that we've finished
    logger.info("Done: %s." % time.asctime())
    logger.info("Time taken: %.2fm" % ((time.time() - t0) / 60))

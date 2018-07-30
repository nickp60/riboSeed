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
import random
import os
import re
import shutil
import multiprocessing
import subprocess
import traceback
import math
import itertools

from bisect import bisect
from itertools import chain
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from distutils.version import StrictVersion


# GLOBALS
SAMTOOLS_MIN_VERSION = '1.3.1'
# --------------------------- classes --------------------------- #


class SeedGenome(object):
    """ This object is the "master" object which holds slots for each
    cluster's mappings, the sequencing library, and keeps track of the
    current iteration as execution progresses.

    When instantiated, self.check_mands() checks that all required attributes
    are present, self.attach_genome_seqRecords() parses and loads the gb data,
    self.write_fasta_genome() writes a .fasta version of the genome for later
    use, and self.make_map_paths_and_dir() sets up the required directories for
    each iteration
    """
    def __init__(self, genbank_path, final_long_reads_dir=None,
                 this_iteration=0, ref_fasta=None, next_reference_path=None,
                 loci_clusters=None, output_root=None, initial_map_bam=None,
                 unmapped_ngsLib=None, name=None, iter_mapping_list=None,
                 reads_mapped_txt=None, unmapped_mapping_list=None,
                 max_iterations=None, initial_map_sam=None, unmapped_sam=None,
                 clustered_loci_txt=None, seq_records=None, master_ngs_ob=None,
                 initial_map_sorted_bam=None, initial_map_prefix=None,
                 assembled_seeds=None, seq_records_count=None, logger=None):
        self.name = name  # get from commsanline in case running multiple
        self.this_iteration = this_iteration  # this should always start at 0
        self.max_iterations = max_iterations
        self.iter_mapping_list = iter_mapping_list  # holds each mapping object
        # The main output for resulting files, all other are relative
        self.output_root = output_root
        # from command line
        self.genbank_path = genbank_path
        # This is created from genbank
        self.ref_fasta = ref_fasta  # this is made dynamically
        # output from riboSelect
        self.clustered_loci_txt = clustered_loci_txt
        # holds a list of LociCluster objects
        self.loci_clusters = loci_clusters
        # this set below with attach_genome_seqRecor
        self.seq_records = seq_records  # this is set
        # this will hold prefix for initial mapping
        self.initial_map_prefix = initial_map_prefix  # set this dynamically
        # extracting reads by position requires an indexed, sorted bam
        self.initial_map_sorted_bam = initial_map_sorted_bam  # set dynamically
        # inial mapping result (combined s and pe)
        self.initial_map_bam = initial_map_bam  # set this dynamically
        # see above
        self.initial_map_sam = initial_map_sam  # set this dynamically
        # holds user-provided sequencing data. Keep intact for final assembly
        self.master_ngs_ob = master_ngs_ob  # for ngslib object
        # each round of seeding results in a ml list for remaining reads
        self.unmapped_mapping_list = unmapped_mapping_list
        # holds dynamically updated list of remaining unmapped
        self.unmapped_sam = unmapped_sam
        # after partitioning, the last mapping list is extracted into this ngs ob
        self.unmapped_ngsLib = unmapped_ngsLib
        # path to file mapped read names are appended to
        self.reads_mapped_txt = reads_mapped_txt
        # destination for seeded contigs prior to final assemblies
        self.final_long_reads_dir = final_long_reads_dir
        # after faux genome construction, store path here
        self.next_reference_path = next_reference_path
        # where to put the combined contigs at the end:
        self.assembled_seeds = assembled_seeds
        # this is set dynamically as well
        self.seq_records_count = seq_records_count
        # a logger
        self.logger = logger
        self.check_mands()
        # self.attach_genome_seqrecords()  # this method comes first,
        self.refresh_seq_rec_generator()  # this method comes first,
        # sets self.seq_records_count
        self.count_seq_records()
        self.write_fasta_genome()  # because this method relies on it
        self.make_map_paths_and_dir()

    def check_mands(self):
        """ check that all mandatory arguments are not none
        """
        mandatory = [self.genbank_path, self.max_iterations,
                     self.output_root, self.clustered_loci_txt]
        if None in mandatory:
            raise ValueError("SeedGenome must be instantiated with at least "
                             "genbank_path, max_iterations, cluster file, " +
                             "and output_root")

    def make_map_paths_and_dir(self):
        """ Given a output root, prepare all the needed subdirs and paths
        """
        self.iter_mapping_list = []
        for i in range(0, self.max_iterations):
            self.iter_mapping_list.append(LociMapping(
                name="{0}_mapping_iteration_{1}".format(self.name, i),
                iteration=i,
                mapping_subdir=os.path.join(
                    self.output_root,
                    "{0}_mapping_for_iteration_{1}".format(self.name, i)),
                assembly_subdir_needed=False))
        if self.final_long_reads_dir is None:
            self.final_long_reads_dir = os.path.join(self.output_root,
                                                     "final_long_reads")
        if not os.path.isdir(self.final_long_reads_dir):
            os.makedirs(self.final_long_reads_dir)

    def write_fasta_genome(self):
        """Given a genbank file, write out as (multi)fasta
        """
        self.name = os.path.splitext(
            os.path.basename(self.genbank_path))[0]
        self.ref_fasta = os.path.join(self.output_root,
                                      str(self.name + ".fasta"))
        with open(self.genbank_path, 'r') as fh:
            with open(self.ref_fasta, 'w') as outfh:
                sequences = SeqIO.parse(fh, "genbank")
                count = SeqIO.write(sequences, outfh, "fasta")
        assert count == self.seq_records_count, "Error parsing genbank file!"
        self.refresh_seq_rec_generator()

    def pad_genbank(self, pad=5000, circular=False, logger=None):
        """ if the genome is circular (which it is, by default) adjust the
        cluster coordinates and rewrite the reference fasta as padded.
        """
        from .shared_methods import pad_genbank_sequence
        if circular:
            for clu in self.loci_clusters:
                clu.padding = pad
                clu = pad_genbank_sequence(cluster=clu, logger=logger)
            logger.info("rewriting fasta-formated genome padded " +
                        "by %i bases on each end", pad)
            new_fasta_ref = os.path.join(
                os.path.basename(self.ref_fasta),
                str(os.path.splitext(self.ref_fasta)[0] +
                    "_padded.fasta"))
            with open(self.genbank_path, 'r') as fh:
                with open(new_fasta_ref, 'w') as outfh:
                    recs = SeqIO.parse(fh, "genbank")
                    count = 0
                    for rec in recs:
                        new_rec = SeqRecord(
                            id=rec.id,
                            seq=Seq(str(rec.seq[-pad:] + rec.seq +
                                    rec.seq[0: pad]),
                                    IUPAC.IUPACAmbiguousDNA()))
                        SeqIO.write(new_rec, outfh, "fasta")
                        count = count + 1
                    assert count == self.seq_records_count, \
                        "Error parsing genbank file!"
            self.ref_fasta = new_fasta_ref
        else:
            pass

    def count_seq_records(self):
        """ Nothing fancy; just a quick way to get the number of records
        """
        self.seq_records_count = sum(1 for x in self.seq_records)
        self.refresh_seq_rec_generator()

    def refresh_seq_rec_generator(self):
        """instead of using stored genbank records, this method restarts
        the generator so that each time self.seq_records is accessed
        this method can get things going again
        """
        self.seq_records = SeqIO.parse(self.genbank_path, "genbank")

    def purge_old_files(self, all_iters=False, logger=None):
        """ remove bulky files from two iterations ago, or if
        all_iters, remove all big mapping files
        """
        assert logger is not None, "Must use logging"
        if all_iters:
            target_iters = range(0, self.max_iterations)
        else:
            target_iters = [self.this_iteration - 2]
            # ensure that you are deleting files from not the previous, but the
            # TWICE previous iterations's files, (which cant be less than 0)
            assert target_iters[0] < self.this_iteration - 1 and \
                target_iters[0] >= 0, \
                "previous mapping is required, can only purge 2nd previous"
        logger.debug("Looking for temp files to remove. " +
                     "To keep all temp files, use the --keep_temps flag")
        for i in target_iters:
            for f in [self.iter_mapping_list[i].pe_map_bam,
                      self.iter_mapping_list[i].s_map_bam,
                      self.iter_mapping_list[i].mapped_sam,
                      # self.iter_mapping_list[i].mapped_bam,  # for riboStack
                      self.iter_mapping_list[i].unmapped_sam,
                      self.iter_mapping_list[i].unmapped_bam,
                      self.iter_mapping_list[i].sorted_mapped_bam]:
                if f is not None:
                    if os.path.isfile(f):
                        os.unlink(f)
                        logger.debug("deleting %s", f)


class NgsLib(object):
    """ NgsLib objects are used to hold the sequencing data supplied by the
    user (master) and the seq data extracted from each iteration. Currently the
    software requires either a paired-end or single-end library, but this
    should handle more diverse library types in the future.

    If ngsLib is master, read lengths are determined and if smalt is used for
    mapping, a distance estimation file is generated.

    when instatiated,self.check_mands() ensure required attributes have values,
    self.set_libtype() sets the libtype attribute based on libraries present,
    self.get_readlen() determines the read length if master, and  if master,
    self.smalt_insert_file() creates a distance estimation file for smalt

    """
    def __init__(self, name, master=False, readF=None, readR=None,
                 readS0=None, readS1=None, mapping_success=False,
                 smalt_dist_path=None, readlen=None, make_dist=False,
                 libtype=None, logger=None, mapper_exe=None, liblist=None,
                 ref_fasta=None):
        self.name = name
        # Bool: whether this is a master record
        self.master = master
        # holds libtype detected from non-None libraries
        self.libtype = libtype  # set this dynamically
        # forward Fastq path
        self.readF = readF
        # reverse Fastq path
        self.readR = readR
        # singleton Fastq path
        self.readS0 = readS0
        # other singleton fastq path (thining ahead, but not really...)
        self.readS1 = readS1
        # detected from Forward fastq with gget_readlen below
        self.readlen = readlen  # set this dynamically
        # bool: did mapping have errors?
        self.mapping_success = mapping_success
        # needed to generate distance file if master
        self.mapper_exe = mapper_exe
        # whether to make a distance file for smalt
        self.make_dist = make_dist
        # also needed to generate distance file if master
        self.ref_fasta = ref_fasta
        # results of distance mapping
        self.smalt_dist_path = smalt_dist_path  # set this dynamically
        self.liblist = liblist  # make dynamically
        self.logger = logger
        self.check_mands()
        self.set_libtype()
        self.get_readlen()
        self.smalt_insert_file()
        self.listLibs()
        self.check_dups()

    def check_mands(self):
        """ checks that all mandatory arguments are not none
        """
        mandatory_other = [self.name, self.ref_fasta]
        if None in mandatory_other:
            raise ValueError("SeedGenome must be instantiated with name, "
                             "and ref fasta")

    def set_libtype(self):
        """sets to either s_1, pe, pe_s based on Non-none libraries
        TODO: add mate support and support for multiple single libraries
        """
        if self.readF is None:
            if self.readR is None:
                if self.readS0 is None:
                    raise ValueError(
                        "Must have at least either a paired library or " +
                        "single library")
                else:
                    self.libtype = "s_1"  # single library
            else:
                raise ValueError("cannot set library type from one PE read")
        else:
            if self.readS0 is not None:
                self.libtype = "pe_s"  # paired end and single
            else:
                self.libtype = "pe"  # just paired end

    def get_readlen(self):
        """ If NgsLib is master, estimate the read lengh using
        get_ave_read_len_from_fastq.
        """
        if self.master is not True:
            return None

        if self.libtype in ['pe', 'pe_s']:
            self.readlen = get_ave_read_len_from_fastq(
                self.readF, N=36, logger=self.logger)
        else:
            self.readlen = get_ave_read_len_from_fastq(
                self.readS0, N=36, logger=self.logger)

    def smalt_insert_file(self):
        """ Smalt mapper uses a subset of mapped reads to estimate distribution
        of insert sizes.  This file is used along with mapping, and is created
        if make_dist, master, and lib_type indicates paired data.
        """
        if self.master is not True:
            return None
        if not self.make_dist:
            return None
        if self.libtype in ['pe', 'pe_s']:
            self.smalt_dist_path = estimate_distances_smalt(
                outfile=os.path.join(os.path.dirname(self.ref_fasta),
                                     "smalt_distance_est.sam"),
                smalt_exe=self.mapper_exe,
                ref_genome=self.ref_fasta,
                fastq1=self.readF, fastq2=self.readR,
                logger=self.logger)
        else:
            print("cannot create distance estimate for lib type %s" %
                  self.libtype)
            return None

    def purge_old_files(self, master, logger=None):
        """ before reasigning unmapped lib, delete
        useless files that were used in the previous iteration
        return codes:
         0) all is well
         1) no files deleted due to conflict
        """
        assert master is not None, \
            str("Must submit your master ngs lib as a control" +
                " to prevent undesired deletions!")
        assert logger is not None, "must use logging"
        if self.master:
            logger.warning("cannot remove master NgsLib!  skipping purge")
            return 1

        for f in [self.readF,
                  self.readR,
                  self.readS0]:
            if f is not None:
                if os.path.isfile(f):
                    if f not in master.liblist:
                        os.unlink(f)
                        logger.debug("Deleted %s", f)
                    else:
                        logger.error(str("%s is in the master ngs object; " +
                                         "cannot delete!"), f)
                        return 1
        return 0

    def listLibs(self):
        self.liblist = [x for x in
                        [self.readF, self.readR, self.readS0, self.readS1]
                        if x is not None]

    def check_dups(self):
        if len(set(self.liblist)) < len(self.liblist):

            raise ValueError(
                "duplicate libraries detected! SPAdes will complain. " +
                "See library list below \n{0}\n. Exiting".format(
                    "\n".join(self.liblist)))


class LociMapping(object):
    """
    instantiate with iteration, mapping subdir, and ref
    map_ref_genome_will use it and an ngslib for mapping
    extract_
    order of operations: map to reference, extract and convert,
    assemble, save results here
    """
    def __init__(self, name, iteration, mapping_subdir, mapping_prefix=None,
                 assembly_success=False,
                 ref_fasta=None, pe_map_bam=None, s_map_bam=None,
                 sorted_mapped_bam=None,
                 mapped_bam=None, mapped_sam=None,
                 mapped_bam_unfiltered=None,
                 unmapped_sam=None,
                 mapped_ids_txt=None, unmapped_bam=None,
                 mappedS=None, assembled_contig=None, assembly_subdir=None,
                 mapped_ngsLib=None, unmapped_ngsLib=None,
                 assembly_subdir_needed=True):
        # int: current iteration (0 is initial)
        self.iteration = iteration
        self.name = name
        # bool: did the eassembly run without errors
        self.assembly_success = assembly_success
        # results for the
        self.mapping_subdir = mapping_subdir
        self.assembly_subdir = assembly_subdir
        self.assembly_subdir_needed = assembly_subdir_needed
        self.ref_fasta = ref_fasta
        # ----  do not ever name these directly ---- #
        self.pe_map_bam = pe_map_bam  # all reads from pe mapping
        self.s_map_bam = s_map_bam  # all reads from singltons mapping
        self.mapping_prefix = mapping_prefix  # set dynamically
        # added to makes filtering with pysam easier
        self.mapped_bam_unfiltered = mapped_bam_unfiltered
        self.mapped_sam = mapped_sam  # mapped reads only, sam
        self.mapped_bam = mapped_bam  # mapped reads only, bam
        self.mapped_ids_txt = mapped_ids_txt
        self.unmapped_sam = unmapped_sam
        self.unmapped_bam = unmapped_bam
        self.sorted_mapped_bam = sorted_mapped_bam  # used with intial mapping
        self.mapped_ngsLib = mapped_ngsLib
        self.unmapped_ngsLib = unmapped_ngsLib
        self.assembled_contig = assembled_contig
        #
        self.check_mands()
        self.make_mapping_subdir()
        self.make_assembly_subdir()
        self.name_bams_and_sams()

    def check_mands(self):
        """ checks that all mandatory arguments are not none
        """
        mandatory = [self.name, self.iteration, self.mapping_subdir]
        if None in mandatory:
            raise ValueError("mapping ob must be instantiated with name, "
                             "iteration, mapping_subdir name")

    def make_mapping_subdir(self):
        """make a subdirectory for mapping results """
        if not os.path.isdir(self.mapping_subdir):
            os.makedirs(self.mapping_subdir)
        else:
            pass

    def name_bams_and_sams(self):
        """ make a prefix and use it to name the future output files """
        mapping_prefix = os.path.join(
            self.mapping_subdir,
            self.name)
        self.pe_map_bam = str(mapping_prefix + "_pe.bam")
        self.s_map_bam = str(mapping_prefix + "_s.bam")
        self.mapped_bam_unfiltered = str(mapping_prefix + "_unfiltered.bam")
        self.mapped_bam = str(mapping_prefix + ".bam")
        # self.unampped_bam = str(mapping_prefix + "unmapped.bam")
        self.sorted_mapped_bam = str(mapping_prefix + "_sorted.bam")
        self.mapped_sam = str(mapping_prefix + ".sam")
        self.unmapped_sam = str(mapping_prefix + "_unmapped.sam")
        self.mapped_ids_txt = str(mapping_prefix + "_mapped.txt")

    def make_assembly_subdir(self):
        """ make a subdirectory for assembly if it is needed """
        if self.assembly_subdir_needed:
            if self.assembly_subdir is not None:
                if not os.path.isdir(self.assembly_subdir):
                    os.makedirs(self.assembly_subdir)
        else:
            pass


class Exes(object):
    """
    given the amount of system tools that riboSeed requires, this object
    holds the paths to the executables after expanding the user-supplied
    path and verifying with shutil.which that the executable is availible
    to the program.


    """
    def __init__(self, python, samtools, method, spades, quast,
                 smalt, bwa, bcftools=None, vcfutils=None, check=True, mapper=None):
        self.python = python
        self.samtools = samtools
        self.method = method
        self.mapper = mapper
        self.spades = spades
        self.quast = quast
        self.smalt = smalt
        self.bwa = bwa
        self.bcftools = bcftools
        self.check = check
        self.check_mands()
        self.set_mapper()
        self.check_expand_mand_exes()
        self.check_expand_opt_exes()

    def check_mands(self):
        """ checks that all mandatory arguments are not none
        """
        mandatory = [self.python, self.spades, self.method,
                     self.samtools, self.bwa]
        assert None not in mandatory, \
            "must instantiate with python, samtools, spades, method!"

    def set_mapper(self):
        """Exes.mapper attribute is set here to avoid further
        "if method =='smalt' clauses later.
        """
        if self.method == "smalt":
            self.mapper = self.smalt
        elif self.method == "bwa":
            self.mapper = self.bwa
        else:
            raise ValueError("Mapping method not found!")

    def check_expand_mand_exes(self):
        """ for each executable, expand wildcards and use shutil.which
        to get full path to executable.  If not found, throw an error
        """
        if self.check:
            for exe in ["mapper", "samtools", "spades",
                        "mapper"]:
                exe_groomed = os.path.expanduser(getattr(self, exe))
                exe_groomed = shutil.which(exe_groomed)
                if exe_groomed is None:
                    raise ValueError("%s not found in PATH!" % exe)
                setattr(self, exe, exe_groomed)
        else:
            pass

    def check_expand_opt_exes(self):
        """ for each executable, expand wildcards and use shutil.which
        to get full path to executable.
        """
        if self.check:
            for exe in ["bcftools", "quast"]:
                if getattr(self, exe) is not None:
                    print(exe)
                    exe_groomed = os.path.expanduser(getattr(self, exe))
                    exe_groomed = shutil.which(exe_groomed)
                    setattr(self, exe, exe_groomed)
        else:
            pass

    def check_spades_python_version(self, logger):
        """   ensure that we have a availible python to
        ruin spades with. See details in method docstrings
        Note this version will be used for quast as well
        """
        from .riboSeed import fiddle_with_spades_exe
        spades_python = fiddle_with_spades_exe(
            spades_exe=self.spades, logger=logger)
        return spades_python


class LociCluster(object):
    """ organizes the clustering process instead of dealing with nested lists
    This holds the whole cluster of one to several individual loci
    """
    newid = itertools.count()

    def __init__(self, sequence_id, loci_list, padding=None,
                 global_start_coord=None, global_end_coord=None,
                 seq_record=None, feat_of_interest=None, mappings=None,
                 extractedSeqRecord=None, cluster_dir_name=None,
                 coverage_exclusion=None,
                 circular=False, output_root=None, final_contigs_path=None,
                 continue_iterating=True, keep_contigs=True):
        # int: unique identifier for cluster
        self.index = next(LociCluster.newid)
        # str: sequence name, usually looks like 'NC_17777373.1' or similar
        self.sequence_id = sequence_id
        # list: hold locus objects for each item in cluster
        self.loci_list = loci_list  # this holds the Locus objects
        # int: bounds ____[___.....rRNA....rRNA..rRNA....__]_________
        self.global_start_coord = global_start_coord
        self.global_end_coord = global_end_coord
        # int: how much to pad sequences y if treating as circular
        self.padding = padding
        # str: feature for filtering: rRNA, cDNA, exon, etc
        self.feat_of_interest = feat_of_interest
        # Bool: treat seqs as circular by padding the ends
        self.circular = circular
        # path: where your cluster-specific output goes
        self.cluster_dir_name = cluster_dir_name  # named dynamically
        # path: where the overall output goes
        self.output_root = output_root
        # list: lociMapping objects that hold mappinging paths
        self.mappings = mappings
        # SeqRecord: holds SeqIO Seqrecord for sequence_id
        self.seq_record = seq_record
        # SeqRecord: holds SeqIO Seqrecord for seq extracted from global coords
        self.extractedSeqRecord = extractedSeqRecord
        # path: for best contig after riboseed2 iterations
        self.keep_contigs = keep_contigs  # by default, include all
        self.continue_iterating = continue_iterating  # by default, keep going
        self.coverage_exclusion = coverage_exclusion
        self.final_contig_path = final_contigs_path
        self.name_mapping_dir()

    def name_mapping_dir(self):
        self.cluster_dir_name = str("{0}_cluster_{1}").format(
            self.sequence_id, self.index)


class Locus(object):
    """ this holds the info for each individual Locus"
    """
    def __init__(self, index, sequence_id, locus_tag, strand=None,
                 start_coord=None, end_coord=None,
                 product=None,
                 feature_type=None):
        # int: unique identifier for cluster
        self.index = index
        # str: sequence name, usually looks like 'NC_17777373.1' or similar
        self.sequence_id = sequence_id  # is this needed? I dont think so as long
        # str: unique identifier from \locus_tag= of gb file
        self.locus_tag = locus_tag
        # int: 1 is +strand, -1 is -strand
        self.strand = strand
        # int:
        self.start_coord = start_coord
        self.end_coord = end_coord
        # str: from \product= of gb file
        self.product = product
        # str "feature" genbank identifier
        self.feature_type = feature_type



# --------------------------- methods --------------------------- #

def get_ave_read_len_from_fastq(fastq1, N=50, logger=None):
    """from LP; return average read length in fastq1 file from first N reads
    """
    count, tot = 0, 0
    if os.path.splitext(fastq1)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    with open_fun(fastq1, "rt") as file_handle:
        data = SeqIO.parse(file_handle, "fastq")
        for read in data:
            count += 1
            tot += len(read)
            if count >= N:
                break
    if logger:
        logger.info(str("From the first {0} reads in {1}, " +
                        "mean length is {2}").format(N,
                                                     os.path.basename(fastq1),
                                                     float(tot / count)))
    file_handle.close()
    return float(tot / count)


def estimate_distances_smalt(outfile, smalt_exe, ref_genome, fastq1, fastq2,
                             cores=None, logger=None):  # pragma: no cover
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

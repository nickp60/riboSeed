#!/usr/bin/env python
"""
Minor version changes:
 - cleaned up
Input:
- genbank file
- dictionary
- specific features : 16S, 5S
- upstream, downstream widths

Output:
-dir containing DNA fastas in their

"""
# import re
import os
# import csv
import subprocess
import datetime
import time
import argparse
import sys
import math
import re

from collections import defaultdict  # for calculating kmer frequency
from itertools import product  # for getting all possible kmers

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# import seaborn as sns
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from heatmapcluster import heatmapcluster

#from pyutilsnrw import utils3_5
from pyutilsnrw.utils3_5 import set_up_logging, check_installed_tools,\
    combine_contigs


class loci_cluster(object):
    """ organizes the clustering process instead of dealing with nested lists
    This holds the whole cluster of one to several individual loci
    """
    def __init__(self, index, sequence, loci_list, padding=None,
                 global_start_coord=None, global_end_coord=None,
                 SeqRecord=None, extractedSeqRecord=None, replace=False,
                 circular=False):
        self.index = index
        self.sequence = sequence
        self.loci_list = loci_list  # this holds the locus objects
        self.global_start_coord = global_start_coord
        self.global_end_coord = global_end_coord
        self.padding = padding
        self.replace = replace
        self.circular = circular
        self.SeqRecord = SeqRecord
        self.extractedSeqRecord = extractedSeqRecord


class locus(object):
    """ this holds the info for each individual locus"
    """
    def __init__(self, index, sequence, locus_tag, strand=None,
                 start_coord=None, end_coord=None, rel_start_coord=None,
                 rel_end_coord=None, product=None):
        # self.parent ??
        self.index = index
        self.sequence = sequence  # is this needed? I dont think so as long
        # as a locus is never decoupled from the loci_cluster
        self.locus_tag = locus_tag
        self.strand = strand  # 1 is +, -1 is -
        self.start_coord = start_coord
        self.end_coord = end_coord
        self.rel_start_coord = rel_start_coord  # start relative to length of region
        self.rel_end_coord = rel_end_coord  # end relative to length of region
        self.product = product


def get_args():  # pragma: no cover
    """get the arguments as a main parser with subparsers
    for named required arguments and optional arguments
    """
    parser = argparse.ArgumentParser(description="Use to extract regions " +
                                     "of interest based on supplied locus " +
                                     " tags.", add_help=False)
    parser.add_argument("genbank_genome", help="Genbank file (WITH SEQUENCE)")
    parser.add_argument("clustered_loci", help="output from riboSelect")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output",
                               help="output directory; default: %(default)s",
                               default=os.getcwd(),
                               type=str, dest="output")

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
    parser.add_argument("-f", "--feature", help="Feature, such as CDS,tRNA, " +
                        "rRNA; default: %(default)s",
                        default='rRNA', dest="feature",
                        action="store", type=str)
    optional.add_argument("-w", "--within_feature_length",
                          help="bp's to include within the region; " +
                          "default: %(default)s",
                          default=0, dest="within", action="store", type=int)
    optional.add_argument("-m", "--minimum_feature_length",
                          help="if --replace, and sequence is shorter than " +
                          "2x --within_feature_length, --within will be " +
                          "modified so that only -m bp of sequnece are" +
                          "turned to N's " +
                          "default: %(default)s",
                          default=100, dest="minimum",
                          action="store", type=int)
    optional.add_argument("-n", "--name",
                          help="rename the contigs with this prefix" +
                          # "default: %(default)s",
                          "default: date (YYYYMMDD)",
                          default=None, dest="name",
                          action="store", type=str)
    optional.add_argument("-l", "--flanking_length",
                          help="length of flanking regions, can be colon-" +
                          "separated to give separate upstream and " +
                          "downstream flanking regions; default: %(default)s",
                          default='1000', type=str, dest="flanking")
    optional.add_argument("-r", "--replace",
                          help="replace sequence with N's; " +
                          "default: %(default)s",
                          default=False, action="store_true", dest="replace")
    optional.add_argument("--msa_kmers",
                          help="calculate kmer similarity based on aligned " +
                          "sequences instead of raw sequences;" +
                          "default: %(default)s",
                          default=False, action="store_true", dest="msa_kmers")
    optional.add_argument("-c", "--circular",
                          help="if the genome is known to be circular, and " +
                          "an region of interest (including flanking bits) " +
                          "extends past chromosome end, this extends the " +
                          "seqence past chromosome origin forward by 5kb; " +
                          "default: %(default)s",
                          default=False, dest="circular", action="store_true")
    optional.add_argument("-p", "--padding", dest='padding', action="store",
                          default=5000, type=int,
                          help="if treating as circular, this controls the " +
                          "length of sequence added to the 5' and 3' ends " +
                          "to allow for selecting regions that cross the " +
                          "chromosom's origin; default: %(default)s")
    optional.add_argument("-v", "--verbosity",
                          dest='verbosity', action="store",
                          default=2, type=int,
                          help="1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("--clobber",
                          help="overwrite previous output files" +
                          "default: %(default)s", action='store_true',
                          default=False, dest="clobber")
    optional.add_argument("--no_revcomp",
                          help="default returns reverse complimented seq " +
                          "if majority of regions on reverse strand. if  " +
                          "--no_revcomp, this is overwridden" +
                          "default: %(default)s",
                          action='store_true',
                          default=False, dest="no_revcomp")
    optional.add_argument("--skip_check",
                          help="Dont bother calculating Shannon Entropy; " +
                          "default: %(default)s",
                          action='store_true',
                          default=False, dest="skip_check")
    optional.add_argument("--msa_tool", dest="msa_tool",
                          choices=["mafft", "prank"],
                          action="store", default="mafft",
                          help="Path to PRANK executable; " +
                          "default: %(default)s")
    optional.add_argument("--prank_exe", dest="prank_exe",
                          action="store", default="prank",
                          help="Path to PRANK executable; " +
                          "default: %(default)s")
    optional.add_argument("--mafft_exe", dest="mafft_exe",
                          action="store", default="mafft",
                          help="Path to MAFFT executable; " +
                          "default: %(default)s")
    optional.add_argument("--barrnap_exe", dest="barrnap_exe",
                          action="store", default="barrnap",
                          help="Path to barrnap executable; " +
                          "default: %(default)s")
    optional.add_argument("--kingdom", dest="kingdom",
                          action="store", default="bac",
                          choices=["mito", "euk", "arc", "bac"],
                          help="kingdom for barrnap; " +
                          "default: %(default)s")
    # had to make this explicitly to call it a faux optional arg
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    args = parser.parse_args()
    return args


def get_genbank_rec_from_multigb(recordID, genbank_records):
    """ given a record ID and and list of genbank records, return sequence of
    genbank record that has all the loci.
    If on different sequences, return error
    """
    for record in genbank_records:
        if recordID == record.id:
            return record
        else:
            pass
    # if none found, raise error
    raise ValueError("no record found matching record id!")


def parse_clustered_loci_file(filepath, gb_filepath,
                              padding, circular, logger=None):
    """Given a file from riboSelect or manually created (see specs in README)
    this parses the clusters and returns a list where [0] is sequence name
    and [1] is a list of loci in that cluster
    As of 20161028, this returns a list of loci_cluster objects!
    """
    if logger is None:
        raise ValueError("logging must be used!")
    if not (os.path.isfile(filepath) and os.path.getsize(filepath) > 0):
        raise ValueError("Cluster File not found!")
    clusters = []
    cluster_index = 0
    # this covers common case where user submits genbank and cluster file
    # in the wrong order.
    if filepath.endswith(("gb", "genbank", "gbk")):
        logger.error("Hmm, this cluster file looks like genbank; " +
                     "it ends in {0}".format(os.path.splitext(filepath)[1]))
        raise FileNotFoundError
    try:
        with open(filepath, "r") as f:
            file_contents = list(f)
    except Exception as e:
        logger.error("Cluster file could not be parsed!")
        raise e
    for line in file_contents:
        try:
            if line.startswith("#") or line.strip() == '':
                continue
            seqname = line.strip("\n").split(" ")[0]
            lt_list = [x for x in
                       line.strip("\n").split(" ")[1].split(":")]
        except Exception as e:
            logger.error("error parsing line: %s" % line)
            raise e
        # make and append the locus objects
        loci_list = []
        for i, loc in enumerate(lt_list):
            loci_list.append(locus(index=i,
                                   locus_tag=loc,
                                   sequence=seqname))
        # make and append loci_cluster objects
        clusters.append(loci_cluster(index=cluster_index,
                                     sequence=seqname,
                                     loci_list=loci_list,
                                     padding=padding,
                                     circular=circular))
        cluster_index = cluster_index + 1
    if len(clusters) == 0:
        raise ValueError("No Clusters Found!!")
    # match up seqrecords
    with open(gb_filepath) as fh:
        gb_records = list(SeqIO.parse(fh, 'genbank'))
    for clu in clusters:
        clu.SeqRecord = get_genbank_rec_from_multigb(
            recordID=clu.sequence,
            genbank_records=gb_records)
    return clusters


def extract_coords_from_locus(cluster,
                              feature="rRNA", logger=None):
    """given a list of locus_tags, return a list of
    [loc_number,[start_coord, end_coord], strand, product,
    locus_tag, record.id]
    20161028 returns a loci_cluster
    """
    if logger is None:
        raise ValueError("logging must be used!")
    loc_number = 0  # index for hits
    locus_tags = [x.locus_tag for x in cluster.loci_list]
    for feat in cluster.SeqRecord.features:
        if not feat.type in feature:
            continue
        logger.debug("found {0} in the following feature : \n{1}".format(
            feature, feat))

        try:
            locus_tag = feat.qualifiers.get("locus_tag")[0]
        except:
            logger.error(str("found a feature ({0}), but there is no" +
                             "locus tag associated with it! Try formatting " +
                             "it by running scanScaffolds.sh").format(feat))
            raise ValueError

        # quick way of checking without using whole object
        if locus_tag in locus_tags:
            # make this_locus point to locus we are adding info to
            this_locus = next((x for x in cluster.loci_list if
                               x.locus_tag == locus_tag), None)
            #  SeqIO makes coords 0-based; the +1 below undoes that
            this_locus.start_coord = feat.location.start.position + 1
            this_locus.end_coord = feat.location.end.position
            this_locus.strand = feat.strand
            this_locus.product = feat.qualifiers.get("product")
            # assert cluster.SeqRecord.id == this_locus.sequence # sanity check, probably not needed
            # loc_list.append([loc_number, coords, strand,
            #                  product, locus_tag, record.id])
            logger.debug("Added attributes for %s", this_locus.locus_tag)
            logger.debug(str(this_locus.__dict__))
            loc_number = loc_number + 1
        else:
            pass
    if not loc_number > 0:
        logger.error("no hits found in any record with feature %s! Double " +
                     "check your genbank file", feature)
        raise ValueError
    logger.debug("Here are the detected region,coords, strand, product, " +
                 "locus tag, subfeatures and sequence id of the results:")
    logger.debug(str(cluster.__dict__))
    return cluster


def pad_genbank_sequence(cluster, logger=None):
    """coords in coords list should be the 1st list item, with the index
    being 0th. Given a genbank record and a coord_list. this returns a seq
    padded on both ends by --padding bp, and returns a coord_list with coords
    adjusted accordingly.  Used to capture regions across origin.
    # as of 20161028, cluster object, not coord_list, is used
    """
    ### take care of the coordinates
    if logger:
        logger.info(str("adjusting coordinates by {0} to account for " +
                    "padding").format(cluster.padding))
    for loc in cluster.loci_list:
        logger.debug("pre-padded")
        logger.debug(str(loc.__dict__))
        start, end = loc.start_coord, loc.end_coord
        loc.start_coord, loc.end_coord = [start + cluster.padding,
                                          end + cluster.padding]
        logger.debug("post-padded")
        logger.debug(str(loc.__dict__))
    ### take care of the sequence
    old_seq = cluster.SeqRecord.seq
    if cluster.padding > len(old_seq):
        raise ValueError("padding cannot be greater than length of sequence")
    new_seq = str(old_seq[-cluster.padding:]
                  + old_seq
                  + old_seq[0: cluster.padding])
    assert len(new_seq) == (len(old_seq) + (2 * cluster.padding)), \
        "Error within function! new seq should be len of " + \
        "seq plus 2x padding"
    cluster.SeqRecord = SeqRecord(Seq(new_seq))
    return cluster


def stitch_together_target_regions(cluster,
                                   flanking="500:500",
                                   within=50, minimum=50, replace=True,
                                   logger=None, verbose=True, circular=False,
                                   revcomp=False):
    """
    given a list from get_coords, usually of length 3 (16,5,and 23 rRNAs),
    return a string with the sequence of the region, replacing coding
    sequences with N's (or not, if replace=False), and including the flanking
    regions upstream and down.

    revamped 20161004
    """
    if logger is None:
        raise ValueError("Must have logger for this function")
    if replace is True:
        # raise ValueError("--replace no longer supported")
        logger.warning("--replace no longer supported")
    try:
        flank = [int(x) for x in flanking.split(":")]
        if len(flank) == 1:  # if only one value use for both up and downstream
            flank.append(flank[0])
        assert len(flank) == 2
    except:
        raise ValueError("Error parsing flanking value; must either be " +
                         "integer or two colon-seapred integers")

    #TODO : make this safer. coord list is constructed sequentially but this
    # is a backup. Throws sort of a cryptic error. but as I said, its a backup
    if sorted([x.start_coord for x in cluster.loci_list]) != \
       [x.start_coord for x in cluster.loci_list]:
        logger.warning("Coords are not in increasing order; " +
                       "you've been warned")
    start_list = sorted([x.start_coord for x in cluster.loci_list])
    logger.debug("Start_list: {0}".format(start_list))
    # if not strictly_increasing([x for x in start_list]):
    #     raise ValueError("coords are not increasing!")
    smallest_feature = min([x.end_coord - x.start_coord for
                            x in cluster.loci_list])
    if smallest_feature < (minimum) and replace:
        raise ValueError("invalid minimum of {0}! cannot exceed half of " +
                         "smallest feature, which is {1} in this case".format(
                             minimum, smallest_feature))

    logger.debug("stitching together the following coords:")
    for i in cluster.loci_list:
        logger.debug(str(i.__dict__))
    #  This works as long as coords are never in reverse order
    cluster.global_start_coord = min([x.start_coord for
                                      x in cluster.loci_list]) - flank[0]
    #
    # if start is negative, just use 1, the beginning of the sequence
    if cluster.global_start_coord < 1:
        logger.warning("Caution! Cannot retrieve full flanking region, as " +
                       "the 5' flanking region extends past start of " +
                       "sequence. If this is a problem, try using a smaller " +
                       "--flanking region, and/or if  appropriate, run with " +
                       "--circular.")
        cluster.global_start_coord = 1
    cluster.global_end_coord = max([x.end_coord for
                                    x in cluster.loci_list]) + flank[1]
    if cluster.global_end_coord > len(cluster.SeqRecord):
        logger.warning("Caution! Cannot retrieve full flanking region, as " +
                       "the 5' flanking region extends past start of " +
                       "sequence. If this is a problem, try using a smaller " +
                       "--flanking region, and/or if  appropriate, run with " +
                       "--circular.")
        cluster.global_end_coord = len(cluster.SeqRecord)

    logger.debug("global start and end: %s %s", cluster.global_start_coord,
                 cluster.global_end_coord)
    #  the minus one makes things go from 1 based to zero based
    seq_with_ns = str(cluster.SeqRecord.seq[cluster.global_start_coord - 1:
                                            cluster.global_end_coord])
    seq_len = len(seq_with_ns[:])
    #
    # loop to mask actual coding regions with N's
    #
    if replace:
        for loc in cluster.loci_list:
            region_length = loc.end_coord - loc.start_coord - (2 * within)
            # if dealing with short sequences
            if region_length < (2 * within):
                # set within to retain minimum sequence length
                this_within = int((region_length - minimum) / 2)
            else:
                # use default if not
                this_within = within
            loc.rel_start_coord = ((loc.start_coord + this_within) -
                                   cluster.global_start_coord)
            loc.rel_end_coord = ((loc.end_coord - this_within) -
                                 cluster.global_start_coord)
            seq_with_ns = str(seq_with_ns[0: loc.rel_start_coord] +
                              str("N" * region_length) +
                              seq_with_ns[loc.rel_end_coord:])

        try:
            # make sure the sequence is proper length, corrected for zero-index
            assert cluster.global_end_coord - cluster.global_start, seq_len
            # make sure replacement didnt change seq length
            assert seq_len, len(seq_with_ns)
        except:
            logger.error("There appears to be an error with the seqeuence " +
                         "coordinate  calculation!")
    # again, plus 1 corrects for 0 index.
    # len("AAAA") = 4 vs AAAA[-1] - AAAA[0] = 3
    logger.info(str("\nexp length {0} \nact length {1}".format(
        cluster.global_end_coord - cluster.global_start_coord + 1, seq_len)))
    # if verbose:
    #     lb = 70  # line break
    #     for i in range(0, int(len(seq_with_ns) / lb)):
    #         print(str(full_seq[i * lb: lb + (i * lb)]))
    #         print(str(seq_with_ns[i * lb: lb + (i * lb)]))
    #         print()
    if not circular:
        seq_id = str(cluster.sequence + "_" + str(cluster.global_start_coord) +
                     ".." + str(cluster.global_end_coord))
    else:  # correct for padding
        seq_id = str(cluster.sequence + "_" + str(cluster.global_start_coord -
                                                  cluster.padding) +
                     ".." + str(cluster.global_end_coord - cluster.padding))

    strands = [x.strand for x in cluster.loci_list]
    # if most are on - strand, return sequence reverse complement
    if revcomp and \
       (sum([x == -1 for x in strands]) > sum([x == 1 for x in strands])):
        logger.info("returning the reverse compliment of the sequence")
        cluster.extractedSeqRecord = SeqRecord(
            Seq(seq_with_ns,
                IUPAC.IUPACAmbiguousDNA()).reverse_complement(),
            id=str(seq_id + "_RC"))
    else:
        cluster.extractedSeqRecord = SeqRecord(
            Seq(seq_with_ns,
                IUPAC.IUPACAmbiguousDNA()),
            id=seq_id)
    return cluster


def prepare_prank_cmd(outdir, combined_fastas, prank_exe,
                      add_args="", outfile_name="best_MSA.fasta",
                      clobber=False, logger=None):
    """returns command line for constructing MSA with
    PRANK and the path to results file
    """
    if logger is None:
        raise ValueError("Must use logger")
    if not os.path.exists(outdir):
        raise FileNotFoundError("output directory not found!")
    prank_cmd = "{0} {1} -d={2} -o={3}".format(
        prank_exe, add_args, combined_fastas,
        os.path.join(outdir, outfile_name))
    logger.debug("PRANK command: \n %s", prank_cmd)
    return (prank_cmd, os.path.join(outdir, outfile_name))


def prepare_mafft_cmd(outdir, combined_fastas, mafft_exe,
                      add_args="", outfile_name="best_MSA",
                      clobber=False, logger=None):
    """returns command line for constructing MSA with
    mafft and the path to results file
    """
    if logger is None:
        raise ValueError("Must use logger")
    if not os.path.exists(outdir):
        raise FileNotFoundError("output directory not found!")
    mafft_cmd = "{0} {1} {2} > {3}".format(
        mafft_exe, add_args, combined_fastas,
        os.path.join(outdir, outfile_name))
    logger.debug("MAFFT command: \n %s", mafft_cmd)
    return (mafft_cmd, os.path.join(outdir, outfile_name))


def calc_Shannon_entropy(matrix):
    """ $j$ has entropy $H(j)$ such that
    $H(j) = -\sum_{i=(A,C,T,G)} p_i(j) \log p_i(j)$
    """
    entropies = []
    for instance in matrix:
        unique = set(instance)
        proportions = {}
        for i in unique:
            proportions[i] = sum([x == i for x in instance]) / len(instance)
        entropy = -sum([prob * (math.log(prob, math.e)) for
                        prob in proportions.values()])
        entropies.append(entropy)
    return entropies


def calc_entropy_msa(msa_path):
    """givn a path to an MSA in FASTA format, this gets the
    $j$ has entropy $H(j)$ such that
    $H(j) = -\sum_{i=(A,C,T,G)} p_i(j) \log p_i(j)$
    return list
    """
    batch_size = 1000  # read seequences in chunks this long
    lengths = []
    with open(msa_path) as fh:
        msa_seqs = list(SeqIO.parse(fh, 'fasta'))
    seq_names = []
    for rec in msa_seqs:
        lengths.append(len(rec))
        seq_names.append(rec.id)
    if not all([i == lengths[0] for i in lengths]):
        raise ValueError("Sequences must all be the same length!")
    entropies = []
    tseq = []
    # save memory by reading in chunks
    for batch in range(0, (math.ceil(lengths[0] / batch_size))):
        # get all sequences into an array
        seq_array = []
        for nseq, record in enumerate(msa_seqs):
            seq_array.append(
                [x for x in record.seq[(batch * batch_size):
                                       ((batch + 1) * batch_size)]])
        # transpose
        tseq_array = list(map(list, zip(*seq_array)))
        tseq.extend(tseq_array)
        entropies.extend(calc_Shannon_entropy(tseq_array))
    # check length of sequence is the same as length of the entropies
    assert len(entropies) == lengths[0]
    return (entropies, seq_names, tseq)


def annotate_msa_conensus(tseq_array, seq_file, barrnap_exe,
                          kingdom="bact",
                          pattern='product=(.+?)$',
                          countdashcov=True,   # include -'s in coverage
                          collapseNs=False,  # include -'s in consensus
                          excludedash=False,
                          logger=None):
    """ returns annotations (as a gfflist),the consensus sequence as a list,
    and named coords  as a list
    TODO: The 'next_best' thing fails is an N is most frequent. Defaults to a T
    """
    if logger is None:
        raise ValueError("Must use logging")
    if excludedash:
        logger.warning("CAUTION: excludedash selected. There is a known " +
                       "bug in the 'next_best' thing fails if an " +
                       "N is most frequent. Defaults to a T")
    consensus = []
    nseqs = len(tseq_array[0])
    logger.info("calc coverage for each of the %i positions", len(tseq_array))
    for position in tseq_array:
        if all([x == position[0] for x in position]):
            if position[0] == '-':
                if collapseNs:
                    continue
                elif not countdashcov:
                    consensus.append([position[0], 0])
            else:
                consensus.append([position[0], nseqs])
        else:
            max_count = 0  # starting count
            nextbest_count = 0
            best_nuc = None
            nextbest_nuc = None
            for nuc in set(position):
                count = sum([nuc == z for z in position])
                # if max count, swap with max to nextbest and update max
                if count > max_count:
                    nextbest_count = max_count
                    max_count = count  # update count if better
                    nextbest_nuc = best_nuc
                    best_nuc = nuc
                else:
                    pass
            # contol whether gaps are allowed in consensus
            if (
                    all([x != '-' for x in position]) and
                    best_nuc == '-' and
                    excludedash):
                # if we dont want n's, choose nextbest
                if nextbest_nuc is None:
                    nextbest_nuc = 't'  # I hate myself for this
                consensus.append([nextbest_nuc, nextbest_count])
            elif best_nuc == '-' and not countdashcov:
                consensus.append([best_nuc, 0])
            else:
                consensus.append([best_nuc, max_count])
    # if any are '-', replace with n's for barrnap
    seq = str(''.join([x[0] for x in consensus])).replace('-', 'n')

    # annotate seq
    with open(seq_file, 'w') as output:
        SeqIO.write(SeqRecord(
            Seq(seq, IUPAC.IUPACAmbiguousDNA())), output, "fasta"),
    barrnap_cmd = "{0} --kingdom {1} {2}".format(barrnap_exe,
                                                 kingdom, seq_file)
    barrnap_gff = subprocess.run(barrnap_cmd,
                                 shell=sys.platform != "win32",
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 check=True)
    results_list = [x.split('\t') for x in
                    barrnap_gff.stdout.decode("utf-8").split("\n")]

    ###  make [name, [start_coord, end_coord]] list
    named_coords = []
    for i in results_list:
        if i[0].startswith("#") or not len(i) == 9:
            logger.debug("skipping gff line: %s", i)
            continue
        m = re.search(pattern, i[8])
        if m:
            found = m.group(1)
        named_coords.append([found, [int(i[3]), int(i[4])]])
    if len(named_coords) == 0:
        raise ValueError(str("Error extracting coords from barrnap gff line" +
                             " %s with pattern %s!" % (str(i), pattern)))
    return (results_list, consensus, named_coords)


def plot_scatter_with_anno(data,
                           names=["Position", "Entropy"],
                           title="Shannon Entropy by Position",
                           consensus_cov=[],
                           anno_list=[],
                           output_prefix="entropy_plot.png"):
    """Given annotation coords [feature, [start, end]],
    consensus cov list ['base', coverage_depth],
    entropy values (list) and consensus sequence
    (same length for consensus_cov and data, no funny business),
    plot out the entropies for each position,
    plot the annotations, and return 0
    """
    if len(consensus_cov) != len(data):
        raise ValueError("data and consensus different lengths!")
    df = pd.DataFrame({names[0]: range(1, len(data) + 1),
                       names[1]: data})  # columns=names)
    df_con = pd.DataFrame(consensus_cov, columns=["base",
                                                  "depth"])
    cov_max_depth = max(df_con['depth'])
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
                                   gridspec_kw={'height_ratios': [4, 1]})
    colors = ['#FF8306', '#FFFB07', '#04FF08', '#06B9FF', '#6505FF', '#FF012F',
              '#FF8306', '#FFFB07', '#04FF08', '#06B9FF', '#6505FF', '#FF012F']
    ax1.set_title(title, y=1.08)
    xmin, xmax = 0, len(data)
    ymin, ymax = -0.1, (max(data) * 1.2)
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
    for index, anno in enumerate(anno_list):
        rect1 = patches.Rectangle(
            (anno[1][0],  # starting x
             ymin),  # starting y
            anno[1][1] - anno[1][0],  # rel x end
            ymax - ymin,  # rel y end
            facecolor=mpl.colors.ColorConverter().to_rgba(
                colors[index], alpha=0.2),
            edgecolor=mpl.colors.ColorConverter().to_rgba(
                colors[index], alpha=0.2))
        rect2 = patches.Rectangle(
            (anno[1][0],  # starting x
             1),  # starting y
            anno[1][1] - anno[1][0],  # rel x end
            cov_max_depth,  # dont -1 beacuse start at 1
            facecolor=mpl.colors.ColorConverter().to_rgba(
                colors[index], alpha=0.2),
            edgecolor=mpl.colors.ColorConverter().to_rgba(
                colors[index], alpha=0.2))
        ax1.add_patch(rect1)
        ax2.add_patch(rect2)
        ax1.text((anno[1][0] + anno[1][1]) / 2,    # x location
                 ymax - 0.48,                      # y location
                 anno[0][0:20],                          # text first 20 char
                 ha='center', color='red', weight='bold', fontsize=10)
    ax1.scatter(x=df["Position"], y=df["Entropy"],
                marker='o', color='black', s=2)
    ax1.set_ylabel('Shannon Entropy')
    ax1.get_yaxis().set_label_coords(-.05, 0.5)
    ax2.set_xlim([xmin, xmax])
    ax2.invert_yaxis()
    ax2.set_ylabel('Consensus Coverage')
    ax2.set_xlabel('Position (bp)')
    ax2.get_yaxis().set_label_coords(-.05, 0.5)
    # ax2.set_ylim([1, cov_max_depth + 1]) #, 1])
    ax2.bar(df_con.index, df_con["depth"],
            width=1, color='darkgrey', linewidth=0, edgecolor='darkgrey')
    # ax2.step(df_con.index, df_con["depth"],
    #          where='mid', color='darkgrey')
    for ax in [ax1, ax2]:
        ax.spines['right'].set_visible(False)
        # ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax1.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')
    ax1.xaxis.set_ticks_position('top')
    ax1.tick_params(axis='y', colors='dimgrey')
    ax2.tick_params(axis='y', colors='dimgrey')
    ax1.tick_params(axis='x', colors='dimgrey')
    ax2.tick_params(axis='x', colors='dimgrey')
    ax1.yaxis.label.set_color('black')
    ax2.yaxis.label.set_color('black')
    ax1.xaxis.label.set_color('black')
    ax2.xaxis.label.set_color('black')
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.set_size_inches(12, 7.5)
    fig.savefig(str(output_prefix + '.png'), dpi=(200))


def get_all_kmers(alph="", length=3):
    """Given an alphabet of charactars, return a list of all permuations
    not actually used, I dont think I'll need this
    """
    mers = [''.join(p) for p in product(alph, repeat=length)]
    return mers


def profile_kmer_occurances(rec_list, alph, k, logger=None):
    """ given a list of seq records, an alphabet, and a value k,
    retrun counts dict of kmer occurances and list of seq names
    """
    counts = defaultdict(list)
    names_list = []
    # part 1: get all kmers from seqs so that we can have equal length lisst
    all_mers = []
    for rec in rec_list:
        all_mers.extend([str(rec.seq).lower()[x: x + k] for
                         x in range(0, len(rec.seq) - k + 1)])
    # unique_mers = list(set().union(all_mers))
    unique_mers = set(all_mers)
    for i in unique_mers:
        counts[i] = []  # initialixe counts dictionary with ker keys
    # part two: count 'em
    for n, rec in enumerate(rec_list):
        logger.info("counting kmer occurances for %s", rec.id)
        names_list.append(rec.id)
        logger.debug("converting to lower")
        string = str(rec.seq).lower()
        logger.debug("getting %imers from seq", k)
        string_mers = [string[x: x + k] for x in range(0, len(string) - k + 1)]
        logger.debug("counting mer occurances")
        # Add counts for those present
        for value in set(string_mers):
            counts[value].append(sum([value == mer for mer in string_mers]))
        # filling in where not found
        for missing in (unique_mers - set(string_mers)):
            counts[missing].append(0)
    return counts, names_list


def plot_pairwise_least_squares(counts, names_list, output_prefix):
    """given a  list of sequence names, a prefix for plot file names,
    and a list of counts from plot_kmer_occurances,
    retruns a pandas df of least squares after plotting heatmaps
    """
    res_list = []
    counts_list = []
    for k, v in counts.items():
        counts_list.append(v)
    # this gives each an index and gets all the pairs
    all_pairs = [[index, value] for index, value in
                 enumerate(product(range(0, len(names_list)), repeat=2))]
    for i in all_pairs:
        this_pairs_diffs = []
        for row in counts_list:
            this_pairs_diffs.append((row[i[1][0]] - row[i[1][1]]) ** 2)
        res_list.append([names_list[i[1][0]],
                         names_list[i[1][1]],
                         sum(this_pairs_diffs)])
    lsdf_wNA = pd.DataFrame(res_list, columns=["locus_1", "locus_2", "sls"])
    wlsdf = lsdf_wNA.pivot(index='locus_1', columns='locus_2', values='sls')
    fig, ax = plt.subplots(1, 1)
    lsdf = lsdf_wNA.fillna(value=0)
    print(wlsdf)
    heatmap = ax.pcolormesh(wlsdf,
                            cmap='Greens')
    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(wlsdf.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(wlsdf.shape[1]) + 0.5, minor=False)
    # # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    plt.xticks(rotation=90)
    # # Set the labels
    ax.set_xticklabels(wlsdf.columns.values, minor=False)
    ax.set_yticklabels(wlsdf.index, minor=False)
    fig.colorbar(heatmap)  # add colorbar key
    plt.tight_layout()  # pad=0, w_pad=5, h_pad=.0)
    fig.set_size_inches(8, 8)
    fig.savefig(str(output_prefix + "heatmap.png"), dpi=(200))
    ####  plot clustered heatmap
    plt.close('all')
    plt.figure(1, figsize=(6, 6))
    h = heatmapcluster(wlsdf.as_matrix(), row_labels=wlsdf.index,
                       col_labels=wlsdf.columns.values,
                       num_row_clusters=2, num_col_clusters=2,
                       label_fontsize=6,
                       xlabel_rotation=-75,
                       cmap=plt.cm.coolwarm,
                       show_colorbar=True,
                       top_dendrogram=True)
    # plt.tight_layout()  # pad=0, w_pad=5, h_pad=.0)
    # fig2.set_size_inches(16, 16)
    plt.savefig(str(output_prefix + "clustered_heatmap.png"), dpi=(200))
    return lsdf_wNA


def main(clusters, genome_records, logger, verbose, within, no_revcomp,
         flanking, replace, output, padding, circular, minimum,
         feature, prefix_name):
    get_rev_comp = no_revcomp is False  # kinda clunky
    for cluster in clusters:  # for each cluster of loci
        # locus_tag_list = cluster[1]
        # get seq record that cluster is  from
        try:
            cluster.SeqRecord = \
                get_genbank_rec_from_multigb(recordID=cluster.sequence,
                                             genbank_records=genome_records)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        # make coord list
        try:
            cluster_with_loci = extract_coords_from_locus(
                cluster=cluster,
                # record=genbank_rec,
                feature=feature,
                # locus_tag_list=[x.locus_tag for x in i.loci_list],
                logger=logger)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        logger.info(str(cluster_with_loci.__dict__))
        if circular:
            cluster_post_pad = pad_genbank_sequence(cluster=cluster_with_loci,
                                                    logger=logger)
        else:
            cluster_post_pad = cluster_with_loci

        #  given coords and a sequnce, extract the region as a SeqRecord
        try:
            cluster_post_stitch =\
                stitch_together_target_regions(cluster=cluster_post_pad,
                                               within=within,
                                               minimum=minimum,
                                               flanking=flanking,
                                               replace=replace,
                                               verbose=False,
                                               logger=logger,
                                               circular=circular,
                                               revcomp=get_rev_comp)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        regions.append(cluster_post_stitch.extractedSeqRecord)
    # after each cluster has been extracted, write out results
    logger.debug(regions)
    output_index = 1
    for i in regions:
        if prefix_name is None:
            filename = str("{0}_region_{1}_{2}.fasta".format(date,
                                                             output_index,
                                                             "riboSnag"))
        else:
            filename = str("{0}_region_{1}_{2}.fasta".format(prefix_name,
                                                             output_index,
                                                             "riboSnag"))
        with open(os.path.join(output, filename), "w") as outfile:
            #TODO make discription work when writing seqrecord
            #TODO move date tag to fasta description?
            # i.description = str("{0}_riboSnag_{1}_flanking_{2}_within".format(
            #                           output_index, args.flanking, args.within))
            SeqIO.write(i, outfile, "fasta")
            outfile.write('\n')
            output_index = output_index + 1
    return regions


if __name__ == "__main__":
    args = get_args()
    output_root = os.path.abspath(os.path.expanduser(args.output))
    # Create output directory only if it does not exist
    try:
        os.makedirs(args.output)
    except FileExistsError:
        # '#' is needed in case streaming output eventually
        print("#Selected output directory %s exists" %
              args.output)
        if not args.clobber:
            print("exiting")
            sys.exit(1)
        else:
            print("# continuing, and risking potential loss of data")
    log_path = os.path.join(output_root,
                            str("{0}_riboSnag_log.txt".format(
                                time.strftime("%Y%m%d%H%M"))))
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=log_path,
                            name=__name__)

    logger.debug("Usage:\n{0}\n".format(str(" ".join([x for x in sys.argv]))))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    # parse cluster file
    try:
        clusters = parse_clustered_loci_file(args.clustered_loci,
                                             gb_filepath=args.genbank_genome,
                                             padding=args.padding,
                                             circular=args.circular,
                                             logger=logger)
    except Exception as e:
        logger.error(e)
        sys.exit(1)
    # parse genbank records
    with open(args.genbank_genome) as fh:
        genome_records = list(SeqIO.parse(fh, 'genbank'))
    regions = []
    logger.info("clusters:")
    logger.info(clusters)
    regions = main(clusters=clusters,
                   genome_records=genome_records,
                   logger=logger,
                   verbose=False, within=args.within,
                   flanking=args.flanking,
                   replace=args.replace,
                   output=args.output,
                   padding=args.padding,
                   circular=args.circular,
                   minimum=args.minimum,
                   prefix_name=args.name,
                   no_revcomp=args.no_revcomp,
                   feature=args.feature)

    # profile_string_kmers(string, size)
    # make MSA and calculate entropy
    if not args.skip_check:
        if args.clobber:
            logger.error("Cannot safely check SMA when --clobber is used!")
            sys.exit(1)
        unaligned_seqs = combine_contigs(contigs_dir=args.output,
                                         pattern="*",
                                         contigs_name="riboSnag_unaligned",
                                         ext=".fasta", verbose=False,
                                         logger=logger)
        if args.msa_tool == "prank":
            if check_installed_tools(executable=args.prank_exe,
                                     hard=False,
                                     logger=logger):
                msa_cmd, results_path = prepare_prank_cmd(
                    outdir=args.output,
                    outfile_name="best_MSA",
                    combined_fastas=unaligned_seqs,
                    prank_exe=args.prank_exe,
                    add_args="",
                    clobber=False, logger=logger)
            else:
                logger.error("Construction of MSA skipped because " +
                             "%s is not a valid executable!", args.prank_exe)
                sys.exit(1)
        elif args.msa_tool == "mafft":
            if check_installed_tools(executable=args.mafft_exe,
                                     hard=False,
                                     logger=logger):
                msa_cmd, results_path = prepare_mafft_cmd(
                    outdir=args.output,
                    outfile_name="best_MSA",
                    combined_fastas=unaligned_seqs,
                    mafft_exe=args.mafft_exe,
                    add_args="",
                    clobber=False, logger=logger)
            else:
                logger.error("Construction of MSA skipped because " +
                             "%s is not a valid executable!", args.mafft_exe)
                sys.exit(1)
        logger.info("Running %s for MSA", args.msa_tool)
        subprocess.run(msa_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)

        seq_entropy, names, tseq = calc_entropy_msa(results_path)
        if args.msa_kmers:
            with open(results_path, 'r') as resfile:
                kmer_seqs = list(SeqIO.parse(resfile, "fasta"))
        else:
            kmer_seqs = regions
        counts, names = profile_kmer_occurances(kmer_seqs,
                                                alph='atcg-',
                                                k=5,
                                                logger=logger)

        mca_df = plot_pairwise_least_squares(
            counts=counts, names_list=names,
            output_prefix=os.path.join(
                args.output,
                "sum_least_squares"))
        # calc_plot_mda(df=mca_df, output_path="entropy_plot.png")
        gff, consensus_cov, annos = annotate_msa_conensus(
            tseq_array=tseq,
            pattern='product=(.+?)$',
            seq_file=os.path.join(
                args.output,
                "test_consensus.fasta"),
            barrnap_exe=args.barrnap_exe,
            kingdom=args.kingdom,
            logger=logger)
        title = str("Shannon Entropy by Position\n" +
                    os.path.basename(
                        os.path.splitext(
                            args.genbank_genome)[0]))

        return_code = plot_scatter_with_anno(
            data=seq_entropy,
            consensus_cov=consensus_cov,
            names=["Position", "Entropy"],
            title=title,
            anno_list=annos,
            output_prefix=os.path.join(
                args.output,
                "entropy_plot"))

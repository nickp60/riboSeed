#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.


# These methods were moved from riboSnag to avoid circular
# import of matplotlib and other horrors

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from .classes import Locus, LociCluster

def parse_clustered_loci_file(filepath, gb_filepath, output_root,
                              circular, padding=1000, logger=None):
    """parses the clusters and returns a list of LociCluster objects

    Given a file from riboSelect or manually created (see specs in README)
    this parses the clusters and returns a list where [0] is sequence name
    and [1] is a list of loci in that cluster

    Args:
        filepath (str): path to loci file form riboSelect
        gb_filepath (str): path to genbank file from riboScan
        output_root (str): output root directory
        circular (bool): whether or not the gneome is circular; if so, the ends of
                         the sequences will be the padded
        padding (int):  how many base pairs to pad width
    Returns:
        (list): list of LociCluster objects
    Raises:
        Exception: one of the various things that can happen with bad files
        ValueError: No clusters found

    """
    assert logger is not None, "logging must be used!"
    clusters = []
    try:
        with open(filepath, "r") as f:
            file_contents = list(f)
    except Exception as e:
        logger.error("Cluster file could not be parsed!")
        raise e
    feature = None
    for line in file_contents:
        try:
            if line.startswith("#$ FEATURE"):
                try:
                    feature = line.split("FEATURE")[1].strip()
                except:
                    raise ValueError("Cannot extract FEATURE from '%s'" % line)
                continue
            elif line.startswith("#") or line.strip() == '':
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
            loci_list.append(Locus(index=i,
                                   locus_tag=loc,
                                   sequence_id=seqname))
        # make and append LociCluster objects
        clusters.append(LociCluster(mappings=[],
                                    output_root=output_root,
                                    sequence_id=seqname,
                                    loci_list=loci_list,
                                    padding=padding,
                                    circular=circular))
        # cluster_index = cluster_index + 1
    ### check feature i;f still none or starts with #$ (ie, no split)
    if feature is None:
        raise ValueError("no feature extracted from coords file! This " +
                         "has been made mandatory 20161108")
    ###
    if len(clusters) == 0:
        raise ValueError("No Clusters Found!!")
    for clu in clusters:
        clu.feat_of_interest = feature
    return clusters


def pad_genbank_sequence(cluster, logger=None, verbose=False):
    """ circular genomes need "padding" to get flanking coordinates near the
    origin

    coords in coords list should be the 1st list item, with the index
    being 0th. Given a genbank record and a coord_list. this returns a seq
    padded on both ends by --padding bp, and returns a coord_list with coords
    adjusted accordingly.  Used to capture regions across origin.
    # as of 20161028, cluster object, not coord_list, is used
    """
    # take care of the coordinates
    if logger:
        logger.debug(str("adjusting coordinates by {0} to account for " +
                        "padding").format(cluster.padding))
    for loc in cluster.loci_list:
        if logger:
            logger.debug("pre-padded")
            logger.debug(str(loc.__dict__))
        start, end = loc.start_coord, loc.end_coord
        loc.start_coord, loc.end_coord = [start + cluster.padding,
                                          end + cluster.padding]
        if logger:
            logger.debug("post-padded")
            logger.debug(str(loc.__dict__))
    #TODO: check for interference with other clusters

    # take care of the sequence
    old_seq = cluster.seq_record.seq
    if cluster.padding > len(old_seq):
        logger.warning("padding cannot be greater than length of " +
                       "sequence! returning original sequence")
        return cluster
    new_seq = str(old_seq[-cluster.padding:]
                  + old_seq
                  + old_seq[0: cluster.padding])
    assert len(new_seq) == (len(old_seq) + (2 * cluster.padding)), \
        "Error within function! new seq should be len of " + \
        "seq plus 2x padding"
    cluster.seq_record = SeqRecord(Seq(new_seq))
    return cluster


def extract_coords_from_locus(cluster, feature="rRNA",
                              logger=None, verbose=False):
    """given a LociCluster object, ammend values
    20161028 returns a LociCluster
    """
    assert logger is not None, "logging must be used!"
    loc_number = 0  # index for hits
    locus_tags = [x.locus_tag for x in cluster.loci_list]
    if verbose:
        logger.debug("Locus tags for cluster %s: %s", cluster.index,
                     " ".join([x for x in locus_tags]))
    for feat in cluster.seq_record.features:
        if not feat.type in feature:
            continue
        if verbose:
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
            logger.debug("Added attributes for %s", this_locus.locus_tag)
            # logger.debug(str(this_locus.__dict__))
            logger.debug(
                "locus_tag: %s; coords: [%i-%i]; sequence_id: %s; " +
                "product: %s", this_locus.locus_tag,
                this_locus.start_coord, this_locus.end_coord,
                this_locus.sequence_id, this_locus.product[0])
            loc_number = loc_number + 1
        else:
            pass
            # logger.error("skipping %s", locus_tag)
    if not loc_number > 0:
        logger.error("no hits found in any record with feature %s! Double " +
                     "check your genbank file", feature)
        raise ValueError


def add_gb_seqrecords_to_cluster_list(cluster_list, gb_filepath):
    """return cluster list with seq_records attached

    We need seqrecords in the lociCluster. So we look through the gb file
    to check for a sequence matching the LociClusters's sequence ID

    Args:
        cluster_list (list): List of LociCluster objects
        gb_filepath (str): oath to genbank file
    Returns:
        (list): the same cluster list, but with their seq records
    Raises:
        None

    """
    # match up seqrecords
    gb_records = SeqIO.index(gb_filepath, 'genbank')
    for clu in cluster_list:
        clu.seq_record = gb_records[clu.sequence_id]
    gb_records.close()
    return cluster_list

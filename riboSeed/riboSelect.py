#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""
Input:
- genbank file
- location:  rRNA (or RRNA, or something similar)

Output:
-text file with sequence name and colon-separated locus_tags,
  one for each region
NC33352.5  ECUMN_23S_3:ECUMN_16S_3:ECUMN_5S_4


USAGE:
 $ python3.5 riboSelect.py genome.gbk --output /path/to/output
"""
import os
import datetime
import argparse
import sys
import jenkspy
from Bio import SeqIO

from .classes import Locus
from .shared_methods import set_up_logging, multisplit


def get_args(test_args=None):  # pragma: no cover
    """get the arguments as a main parser with subparsers
    for named required arguments and optional arguments
    """
    parser = argparse.ArgumentParser(prog="ribo select",
                                     description="This is used to identify" +
                                     " and cluster rRNA regions from a gb " +
                                     "file, returns" +
                                     "a text file with the clusters")
    parser.prog = "ribo select"
    parser.add_argument("genbank_genome", help="Genbank file (WITH SEQUENCE)")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output",
                               help="output directory;"
                               "default: %(default)s",
                               required=True,
                               type=str, dest="output")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-f", "--feature",
                          help="Feature, rRNA or RRNA; default: %(default)s",
                          default='rRNA', type=str)
    optional.add_argument("-s", "--specific_features",
                          help="colon:separated -- specific features"
                          "; default: %(default)s",
                          default='16S:23S:5S', type=str)
    optional.add_argument("--clobber",
                          help="overwrite previous output files: "
                          "default: %(default)s", action='store_true',
                          default=False, dest="clobber")
    # TODO Implement this for work with really fragmented genomes
    # optional.add_argument("--nocluster",
    #                       help="do not bother clustering; treat all " +
    #                       "occurances on a sequence as a cluster" +\
    #                       "default: %(default)s", action='store_true',
    #                       default=False, dest="nocluster")
    optional.add_argument("-c", "--clusters",
                          help="number of rDNA clusters;"
                          "if submitting multiple records, must be a "
                          "colon:separated list whose length matches number "
                          "of genbank records.  Default is inferred from "
                          "specific feature with fewest hits", default='',
                          type=str, dest="clusters")
    optional.add_argument("-v", "--verbosity",
                          dest='verbosity', action="store",
                          default=2, type=int,
                          help="1 = debug(), 2 = info(), 3 = warning(), "
                          "4 = error() and 5 = critical(); "
                          "default: %(default)s")
    optional.add_argument("--debug", dest="debug", action="store_true",
                          help="Enable debug messages", default=False)
    if test_args is None:
        args = parser.parse_args(sys.argv[2:])
    else:
        args = parser.parse_args(test_args)
    return args


def count_feature_hits_per_sequence(all_feature, rec_id_list,
                       specific_features, loci_list, logger=None):
    """count rDNA features in each sequence

    given genome seq records, specific_features, and a
    locus_tag_dict from get_filtered_locus_tag_dict, return two structures:
    - nfeatures_simple {record.id: [5,4,6]}
    SPECIFIC FEATURES MUST BE SORTED

    Args:
        all_feature (bool): if counting all features, exit, return none
        rec_id_list (list): list of genome accessions
        specific_features (list): list of features as strings to observe
        loci_list (list): see get_loci_list_for_features
    Returns:
        (dict): {sequenceID: [n_16s, n_23s]}
        None: if all_features
    Raises:
        None

    """
    logger.info("counting occurances of %s", " ".join(specific_features))
    assert logger is not None, "Must use logging"
    # if looking at all the features, ignore this whole bit
    if all_feature:
        return None
    #  This bit counts the number of hits per specific feature.
    logger.info("counting features")
    nfeat_simple = {}  # makes  {genome : [5, 3]}
    for rec_id in rec_id_list:
        # create the entry
        hit_list = [0 for x in specific_features]
        for locus in loci_list:
            if locus.sequence_id == rec_id:
                for idx, x in enumerate(specific_features):
                    if locus.product == x:
                        hit_list[idx] = hit_list[idx] + 1
                    else:
                        pass
        nfeat_simple[rec_id] = hit_list
    return nfeat_simple


def get_loci_list_for_features(gb_path, feature="rRNA",
                                specific_features="16S:23S",
                                verbose=True, logger=None):
    """get list of loci matching specific features

    Given a genbank path, returns a tuple of loci_list and
    dict of {"sequence" [counts of features]}

    The list of Locus object for each feature of interest
    (ie, all rRNA annotations). dictionary of index:locus_tag id pairs for all
    "feature"  entries.  This then gets clustered.
    requires having locus tag in your genbank file.  Non-negitable.
    should be prokka-friendly, so if you have a gb file with legacy
    annotations, just run through prokka (with an rRNA caller)
    20160922 This was changed to add checking for rRNA type from
    ribosome product annoation, not ht locus tag.
    if specific features is None, return all with type == feature
    Changed from having index used for clustering to using the first coord

    Args:
    gb_path (str): path to genbank file
        feature (str): default="rRNA"; the genbank feature type
        specific_features (str): default="16S:23S"; comma sep list of subunits
                                 to consider
        verbose (bool): default=True
    Returns:
        (tuple):
          (list): loci_list,
          (dict): {accession: [n_16, n_23]}

    Raises:
        SystemExit: Error with the annotations

    """
    assert logger is not None, "must use logging, even if not 'verbose'"
    assert os.path.exists(gb_path), \
        'genbank file not found! '
    if specific_features is not None:
        specific_features = set(specific_features.split(":"))
    loci_list = []
    gb_record_ids = []
    # loop through records
    with open(gb_path, "r") as genome_seq_records:
        for record in SeqIO.parse(genome_seq_records, "genbank"):
            gb_record_ids.append(record.id)
            loc_number = 0  # counter for hits
            logger.debug("scanning %s", record.id)
            for feat in record.features:
                if feat.type in feature:
                    try:
                        locustag = feat.qualifiers.get("locus_tag")[0]
                    except TypeError:
                        if verbose:
                            logger.debug("no locus tag for this feature!")
                        continue
                    product_list = set(multisplit(
                        [",", " ", "-", "_"],
                        feat.qualifiers.get("product")[0]))
                    if verbose:
                        logger.debug("%s:%s", locustag, " ".join(product_list))
                    product_matches = \
                        product_list.intersection(specific_features)
                    if len(product_matches) > 1:
                        logger.error(
                            "multiple hits found in a single product " +
                            "annotation; this is likely indicative of " +
                            "indistinct annotation. Please use with riboScan" +
                            "or barrnap outputted genbank file to avoid this.")
                        sys.exit(1)
                    # if either specific feature is found in product or
                    # only interested in all features, add locus to dict
                    if (specific_features is None or
                        (specific_features is not None and product_matches)):
                        # key is start coord
                        loci_list.append(Locus(
                            index=loc_number,
                            sequence_id=record.id,
                            locus_tag=locustag,
                            feature_type=feat.type,
                            # cause sets cant be indexed
                            product=next(iter(product_matches)),
                            start_coord=feat.location.start.position + 1,
                            end_coord=feat.location.end.position,
                            strand=None))
                        loc_number = loc_number + 1
                    else:
                        if verbose:
                            logger.debug("Not adding this feat to " +
                                         "list: %s", product_list)
                else:
                    if verbose:
                        logger.debug("skipping: %s", feat.type)
            # this is a soft warning, as we want to be able to loop
            # through all records before worrying
            if len(loci_list) < 1:
                logger.info(str("no locus tags found in {0} for {1} " +
                                "features with annotated products matching" +
                                "{2}!\n").format(record.id,
                                                 feature, specific_features))
    if verbose:
        for idx, locus in enumerate(loci_list):
            logger.debug("%s: %s;", idx, str(locus.__dict__))
    # count the occuraces of each feature per genbank record
    nfeat_simple = count_feature_hits_per_sequence(
        all_feature=specific_features is None,
        rec_id_list=gb_record_ids,
        specific_features=specific_features,
        loci_list=loci_list,
        logger=logger)
    return (loci_list, nfeat_simple)


def parse_args_clusters(clusters, nrecs, logger=None):
    """determine number of centers based on cluster arg

    default case, clusters are inferred
    if not, must be equal to the length of genbank records

    Args:
        clusters (str): argparse --clusters; to be coerrced to an int or
                        list of ints
    Returns:
        (int): number of centers to use
    Raises:
        SystemExit: cant coerce to int after splitting

    """

    assert logger is not None, "logging must be used"
    logger.info("parsing --clusters args")
    if clusters != "":
        try:
            centers = [int(x) for x in clusters.split(":")]
            logger.info(str(centers))
        except:
            logger.error("cannot coerce --clusters to integer after " +
                         "splitting on colons!\n")
            sys.exit(1)
    else:
        centers = [0 for x in range(0, nrecs)]
    return centers


def dict_from_jenks(data, centers, logger=None):
    """use jenks natural breaks to split a list of coordinates accordingly

    Jenks natural breaks is used to split a list of coordinates at N centers
    to group neighboring coords

    Args:
        data (list): list of coordinates (int)
        centers (int): number of breaks to use
    Returns:
        (dict): {1: [3000, 4330, 4100],
                 2: [6000, 6200, 6143]}
    Raises:
        None

    """
    assert logger is not None, "must use logging"
    assert centers is not 0, "cannot use 0 center"
    if len(data) <= centers:
        logger.info("making 1 cluster per index!")
        return {str(i + 1): [j] for i, j in enumerate(data)}
    if centers == 1:
        return {"1": data}

    # if centers <= 1:
    #     logger.error("center value must be larger than 1")
    breaks = jenkspy.jenks_breaks(data, nb_class=centers)
    logger.debug("Jenks breaks:")
    logger.debug(breaks)
    newdata = {}
    data2 = data[:]  # we will be removing items, so lets keep track
    for idx, b in enumerate(breaks):
        if idx == 0:
            continue
        for point in sorted(data2):
            if point <= b:
                if str(idx) in newdata:
                    newdata[str(idx)].append(point)
                else:
                    newdata[str(idx)] = [point]
                data2.remove(point)

    final_dict = {}
    for k, v in newdata.items():
        final_dict[k] = [x for x in v]
    return(final_dict)


def main(args, logger=None):
    output_root = os.path.abspath(os.path.expanduser(args.output))
    # Create output directory only if it does not exist
    try:
        os.makedirs(output_root)
    except FileExistsError:
         # leading comment char'#' added for stream output usage
        print("#Selected output directory %s exists" %
              output_root)
        sys.exit(1)
    if logger is None:
        logger = set_up_logging(
            verbosity=args.verbosity,
            outfile=os.path.join(output_root, "riboSelect.log"),
            name=__name__)

    logger.info("Usage:\n%s\n", " ".join([x for x in sys.argv]))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, v)
    date = str(datetime.datetime.now().strftime('%Y%m%d'))

    output_path = os.path.join(output_root,
                               str("riboSelect_grouped_loci.txt"))
    if os.path.splitext(args.genbank_genome)[1] in [".fa", ".fasta"]:
        logger.error("Input must be a genbank file; this appears to " +
                     "be a fasta ")
        sys.exit(1)

    # get genome records into a generator to count them
    genome_records = SeqIO.parse(args.genbank_genome, 'genbank')
    logger.info("Getting number of genbank entries")
    nrecs = sum(1 for x in genome_records)
    logger.info("Input has %i records", nrecs)

    # get list of loci matching feature and optionally specific features
    # also returns nfeat, a dict of feature count by genbank id
    logger.debug(
        str(
            "searching for {0} features containing {1} in the " +
            "product annotation").format(
            args.feature,
            str([x for x in args.specific_features.split(":")])))
    loci_list, nfeat_simple = \
        get_loci_list_for_features(gb_path=args.genbank_genome,
                                   feature=args.feature,
                                   specific_features=args.specific_features,
                                   verbose=args.debug,
                                   logger=logger)

    # default case, clusters are inferred
    # if not, must be equal to the length of genbank records
    centers_per_seq = parse_args_clusters(
        clusters=args.clusters, nrecs=nrecs,
        logger=logger)
    # if unequal lengths, throw error
    # logger.info clusters for accession for user to verify
    if nrecs != len(centers_per_seq):
        logger.error("clusters must be the same length as number" +
                     " of genbank records!\n")
        sys.exit(1)

    #####
    ##### for each genbank record, process, and append any hits to outfile
    #####
    logger.debug("All loci:\n" +
                 "\n".join([str(x.__dict__) for x in loci_list]))
    outlines = []
    with open(args.genbank_genome, "r") as gb:
        for i, rec in enumerate(SeqIO.parse(gb, 'genbank')):
            logger.info("Processing {0}\n".format(rec.id))
            # if user gives clusters, make sure it matches the length:
            if args.clusters:
                logger.info("using %s clusters for %s\n",
                            centers_per_seq[i], rec.id)
            subset = [x for x in loci_list if x.sequence_id == rec.id]
            # skip if that doesnt have any hits
            if len(subset) == 0:
                logger.info("no hits in {0}\n".format(rec.id))
                continue
            logger.debug("Subset loci:\n" +
                         "\n".join([str(x.__dict__) for x in subset]))

            #  find nfeat for this genbank id by subsetting;
            # is this a bad way of doesnt things?

            # logger.debug("centers: {0}".format(centers_per_seq))
            # logger.debug("nfeat_simple: {0}".format(nfeat_simple[i]))
            if nfeat_simple is None and centers_per_seq[i] == 0:
                logger.error(
                    "No specific features submitted, cannot calculate " +
                    " number centers needed for clustering.  Please" +
                    "submit the desired number of clusters with the  " +
                    "--clusters argument!\n")
                sys.exit(1)
            logger.debug("nfeat_simple:\n%s", nfeat_simple)
            rec_nfeat = [v for k, v in nfeat_simple.items() if
                              rec.id in k][0]
            logger.debug("rec_nfeat: {0}".format(rec_nfeat))
            if all([x == 0 for x in rec_nfeat]):
                logger.error("unable to count features!")
                sys.exit(1)
            indexes = [x.start_coord for x in subset]  # get index back from tuple key
            ## if centers[i] is 0, try max and min sequentially; if that fails skip
            if centers_per_seq[i] == 0:
                # if only looking at two features, take the max
                # if the smallest value is not greater than 1, take the max
                # This is a shaky heuristic
                if min(rec_nfeat) <= 1 or len(rec_nfeat) <= 2:
                    current_centers = max(rec_nfeat)
                else:
                    current_centers = min(rec_nfeat)
                if current_centers == 0:
                    logger.info("skipping the clustering for {0}\n".format(i))
                    continue
                # logger.debug("grouping indexes with jenks natural breaks " +
                #              "assuming %d classes", current_centers)
                # logger.debug(indexes)
                indexClusters = dict_from_jenks(
                    data=indexes, centers=current_centers, logger=logger)
            else:
                ## if centers[i] is not 0, use it
                current_centers = centers_per_seq[i]
            # Perform actual clustering
            try:
                # indexClusters should be like { "1": [3,4,6], "2": [66,45,63]}
                logger.debug(
                    "grouping the following indexes with jenks " +
                    "natural breaks assuming %d classes", current_centers)
                logger.debug(indexes)
                indexClusters = dict_from_jenks(
                    data=indexes, centers=current_centers, logger=logger)
            except Exception as e:
                logger.error(e)
                sys.exit(1)
            logger.debug("indexClusters:")
            logger.debug(indexClusters)
            # add output lines to list
            outlines.append("# Generated cluters for {0} on {1}\n".format(
                rec.id, date))
            outlines.append("#$ FEATURE {0}\n".format(args.feature))

            for k, v in sorted(indexClusters.items()):
                outlines.append(
                    "{0} {1}\n".format(
                        rec.id,
                        # ":".join([subset[(rec.id, int(x))][2] for x in v]))
                        ":".join([x.locus_tag for x in subset if \
                                  x.start_coord in v]))
                )

    with open(output_path, "a") as outfile:
        for i in outlines:
            sys.stdout.write(i)
            outfile.write(i)

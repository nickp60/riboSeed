#!/usr/bin/env python3.5
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
import csv
import subprocess
import datetime
import argparse
import sys
import time
import jenkspy
from Bio import SeqIO

from pyutilsnrw.utils3_5 import set_up_logging, multisplit


def get_args():  # pragma: no cover
    """get the arguments as a main parser with subparsers
    for named required arguments and optional arguments
    """
    parser = argparse.ArgumentParser(description="This is used to identify" +
                                     " and cluster rRNA regions from a gb " +
                                     "file, returns" +
                                     "a text file with the clusters")
    parser.add_argument("genbank_genome", help="Genbank file (WITH SEQUENCE)")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output",
                               help="output directory;"
                               "default: %(default)s",
                               default=os.getcwd(),
                               type=str, dest="output")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-f", "--feature",
                          help="Feature, rRNA or RRNA; default: %(default)s",
                          default='rRNA', type=str)
    optional.add_argument("-s", "--specific_features",
                          help="colon:separated -- specific features"
                          "; default: %(default)s",
                          default='16S:23S:5S', type=str)
    optional.add_argument("--keep_temps",
                          help="view intermediate clustering files"
                          "default: %(default)s",
                          action='store_true',
                          default=False, dest="keep_temps")
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
    args = parser.parse_args()
    return args


def count_feature_hits(all_feature, gb_path,
                       specific_features, locus_tag_dict, logger=None):
    """ given genome seq records, specific_features, and a
    locus_tag_dict from get_filtered_locus_tag_dict, return two structures:
    - nfeatures_occur {record.id, [[16s, 5],[23s, 4],[5s,6]]}
    - nfeatures_simple {record.id [record.id, [5,4,6]}
    """
    logger.info("counting occurances of %s", " ".join(specific_features))
    assert logger is not None, "Must use logging"
    #  This bit counts the number of hits per specific feature.
    logger.info("counting features")
    if not all_feature:
        nfeatures_occur = {}  # makes  {genome : [['18S', 5],['28S',3]]}
        nfeat_simple = {}  # makes  {genome : [5, 3]}
        with open(gb_path, "r") as genome_seq_records:
            for record in SeqIO.parse(genome_seq_records, "genbank"):
                logger.debug("counting hits in %s", record.id)
                hit_list = []  # [specific feature, count] list
                hit_list_simple = []  # [count] list
                subset = {k: v for k, v in locus_tag_dict.items()
                          if record.id in v}
                if len(subset) == 0:
                    continue
                logger.debug(subset)
                for i in specific_features:
                    hits = 0
                    for k, v in subset.items():
                        # hint: v[-1] should be the product annotation
                        if any([i == x for x in v[-1]]):
                            hits = hits + 1
                        else:
                            pass
                    hit_list.append([i, hits])
                    hit_list_simple.append(hits)
                nfeatures_occur[record.id] = (hit_list)
                nfeat_simple[record.id] = hit_list_simple
    else:
        nfeatures_occur, nfeat_simple = None, None
    return(nfeatures_occur, nfeat_simple)


def get_filtered_locus_tag_dict(gb_path, nrecs, feature="rRNA",
                                specific_features="16S:23S",
                                verbose=True, logger=None):
    """ Given a LIST (as of 20160927) of genbank records,
    returns dictionary of index:locus_tag id pairs for all
    "feature"  entries.  This then gets clustered.
    requires having locus tag in your genbank file.  Non-negitable.
    should be prokka-friendly, so if you have a gb file with legacy
    annotations, just run through prokka (with an rRNA caller)
    20160922 This was changed to add checking for rRNA type from
    ribosome product annoation, not ht locus tag.
    if specific features is None, return all with type == feature
    Changed from having index used for clustering to using the first coord
    """
    assert os.path.exists(gb_path), \
        'genbank file not found! '
    # if specific_features is None:
    #     all_feature = True
    # else:
    #     all_feature = False
    if specific_features is not None:
        specific_features = specific_features.split(":")
    locus_tag_dict = {}  # recipient structure
    preunique_feats = []

    # loop through records
    with open(gb_path, "r") as genome_seq_records:
        # genome_records = list(SeqIO.parse(fh, 'genbank'))
        for record in SeqIO.parse(genome_seq_records, "genbank"):
            loc_number = 0  # counter that for index of hits; this is clustered
            logger.debug("scanning %s", record.id)
            for feat in record.features:
                if feat.type in feature:
                    try:
                        locustag = feat.qualifiers.get("locus_tag")[0]
                    except TypeError:
                        if verbose:
                            logger.debug("no locus tag for this feature!")
                        continue
                    product_list = multisplit(
                        [",", " ", "-", "_"],
                        feat.qualifiers.get("product")[0])
                    coords = [feat.location.start.position + 1,
                              feat.location.end.position]
                    if verbose:
                        logger.debug(product_list)
                        logger.debug(coords)
                    # if either specific feature is found in product or
                    # only interested in all features, add locus to dict
                    if (specific_features is None or
                        (specific_features is not None and
                         any([x in specific_features for x in product_list]))):
                        # key is start coord
                        preunique_feats.extend([x for x in specific_features if
                                                x in product_list])
                        locus_tag_dict[(record.id, coords[0])] = [loc_number,
                                                                  record.id,
                                                                  locustag,
                                                                  feat.type,
                                                                  product_list]
                    else:
                        if verbose:
                            logger.debug("Not adding this feat to " +
                                         "list: %s", product_list)
                        pass

                    loc_number = loc_number + 1  # increment idx after feature
                else:
                    if verbose:
                        logger.debug("skipping: %s", feat.type)
                    loc_number = loc_number + 1  # increment idx after feature
            # this is a soft warning, as we want to be able to loop
            # through all records before worrying
            if len(locus_tag_dict) < 1 and logger:
                logger.info(str("no locus tags found in {0} for {1} " +
                                "features with annotated products matching" +
                                "{2}!\n").format(record.id,
                                                 feature, specific_features))
    # locus_tag_dict = locus_tag_dict
    if verbose and logger:
        for key in sorted(locus_tag_dict):
            logger.debug("%s: %s;", key, locus_tag_dict[key])

    # count the occuraces of each feature per genbank record
    nfeatures_occur, \
        nfeat_simple = count_feature_hits(
            all_feature=specific_features is None,
            gb_path=gb_path,
            specific_features=set(preunique_feats),
            locus_tag_dict=locus_tag_dict,
            logger=logger)
    return(locus_tag_dict, nfeatures_occur, nfeat_simple)


def parse_args_clusters(clusters, nrecs, logger=None):
    """
    """
    # default case, clusters are inferred
    # if not, must be equal to the length of genbank records
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
    if logger:
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


if __name__ == "__main__":
    args = get_args()
    output_root = os.path.abspath(os.path.expanduser(args.output))
    # Create output directory only if it does not exist
    try:
        os.makedirs(output_root)
    except FileExistsError:
         # leading comment char'#' added for stream output usage
        print("#Selected output directory %s exists" %
              output_root)
        sys.exit(1)
    logger = set_up_logging(
        verbosity=args.verbosity,
        outfile=os.path.join(output_root, "riboSelect.log"),
        name=__name__)

    # log = sys.stderr.write  # to keep streaming clean if this goes that route
    logger.info("Current usage:\n{0}\n".format(" ".join(sys.argv[1:])))
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
    lociDict, nfeat, nfeat_simple = \
        get_filtered_locus_tag_dict(gb_path=args.genbank_genome,
                                    feature=args.feature,
                                    specific_features=args.specific_features,
                                    verbose=args.debug,
                                    nrecs=nrecs,
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
    logger.debug("All loci:")
    for k, v in lociDict.items():
        logger.debug(str(k) + "\t" + str(v))
    outlines = []
    with open(args.genbank_genome, "r") as gb:
        for i, rec in enumerate(SeqIO.parse(args.genbank_genome, 'genbank')):
            logger.info("Processing {0}\n".format(rec.id))
            # if user gives clusters, make sure it matches the length:
            if args.clusters:
                logger.info("using %s clusters for %s\n",
                            centers_per_seq[i], rec.id)
            # get subset of lociDict for that id
            subset = {key: value for key, value in lociDict.items() if
                      rec.id == value[1]}
            # skip if that doesnt have any hits
            if len(subset) == 0:
                logger.info("no hits in {0}\n".format(rec.id))
                continue
            logger.debug("Subset loci:")
            for k, v in subset.items():
                logger.debug(str(k) + "\t" + str(v))
            logger.info("hits in {0}\n".format(rec.id))

            #  find nfeat for this genbank id by subsetting;
            # is this a bad way of doesnt things?

            # logger.debug("centers: {0}".format(centers_per_seq))
            # logger.debug("nfeat_simple: {0}".format(nfeat_simple[i]))

            if nfeat_simple is None and centers_per_seq[i] == 0:
                logger.error("No specific features submitted, cannot calculate " +
                             " number centers needed for clustering.  Please" +
                             "submit the desired number of clusters with the  " +
                             "--clusters argument!\n")
                sys.exit(1)
            logger.debug("nfeat_simple")
            # for k, v in nfeat_simple.items():
            #     print(k)
            #     print(v)
            rec_nfeat = list({k: v for k, v in nfeat_simple.items() if
                              rec.id in k}.values())[0]
            logger.debug("rec_nfeat: {0}".format(rec_nfeat))
            if all([x == 0 for x in rec_nfeat]):
                logger.error("unable to count features!")
                sys.exit(1)
            indexes = [x[1] for x in list(subset)]  # get index back from tuple key
            ## if centers[i] is 0, try max and min sequentially; if that fails skip
            if centers_per_seq[i] == 0:
                # if only looking at two features, take the max
                # if the smallest value is not greater than 1, take the max
                # This is a shakey heuistic
                if min(rec_nfeat) <= 1 or len(rec_nfeat) <= 2:
                    current_centers = max(rec_nfeat)
                else:
                    current_centers = min(rec_nfeat)
                if current_centers == 0:
                    logger.info("skipping the clustering for {0}\n".format(i))
                    continue
                logger.debug("grouping indexes with jenks natural breaks " +
                             "assuming %d classes", current_centers)
                logger.debug(indexes)
                indexClusters = dict_from_jenks(
                    data=indexes, centers=current_centers, logger=logger)
            else:
                ## if centers[i] is not 0, use it
                current_centers = centers_per_seq[i]
            # Perform actual clustering
            try:
                # indexClusters should be like { "1": [3,4,6], "2": [66,45,63]}
                logger.debug("grouping indexes with jenks natural breaks " +
                             "assuming %d classes", current_centers)
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

            for k, v in indexClusters.items():
                outlines.append(
                    "{0} {1}\n".format(
                        rec.id,
                        ":".join([subset[(rec.id, int(x))][2] for x in v]))
                )

    with open(output_path, "a") as outfile:
        for i in outlines:
            sys.stdout.write(i)
            outfile.write(i)

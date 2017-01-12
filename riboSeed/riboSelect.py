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
from Bio import SeqIO

from pyutilsnrw.utils3_5 import set_up_logging, multisplit


def get_args():  # pragma: no cover
    """get the arguments as a main parser with subparsers
    for named required arguments and optional arguments
    """
    parser = argparse.ArgumentParser(description="This is used to extract" +
                                     " rRNA regions from a gb file, returns" +
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
                          "colon:separated list that matches number "
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


def count_feature_hits(all_feature, genome_seq_records,
                       specific_features, locus_tag_dict):
    """ given genome seq records, specific_features, and a
    locus_tag_dict from get_filtered_locus_tag_dict, return two structures:
    - nfeatures_occur {record.id, [[16s, 5],[23s, 4],[5s,6]]}
    - nfeatures_simple {record.id [record.id, [5,4,6]}
    """
    #  This bit counts the number of hits per specific feature.
    if not all_feature:
        nfeatures_occur = {}  # makes  {genome : [['18S', 5],['28S',3]]}
        nfeat_simple = {}  # makes  {genome : [5, 3]}
        for record in genome_seq_records:
            hit_list = []  # [specific feature, count] list
            hit_list_simple = []  # [count] list
            for i in specific_features:
                hits = 0
                subset = {k: v for k, v in locus_tag_dict.items()
                          if record.id in v}
                for k, v in subset.items():
                    # hint: v[-1] should be the product annotation
                    if any([i in x for x in v[-1]]):
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


def get_filtered_locus_tag_dict(genome_seq_records, feature="rRNA",
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
    assert isinstance(genome_seq_records, list), \
        'must pass list of genomes to function, even if single genome '
    # if specific_features is None:
    #     all_feature = True
    # else:
    #     all_feature = False
    if specific_features is not None:
        specific_features = specific_features.split(":")
    locus_tag_dict = {}  # recipient structure
    # loop through records
    for record in genome_seq_records:
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
                product_list = multisplit([",", " ", "-", "_"],
                                          feat.qualifiers.get("product")[0])
                coords = [feat.location.start.position + 1,
                          feat.location.end.position]
                logger.debug(product_list)
                logger.debug(coords)
                # if either specific feature is found in product or
                # only interested in all features, add locus to dict
                if (specific_features is None or
                    (specific_features is not None and
                     any([x in specific_features for x in product_list]))):
                    # key is start coord
                    locus_tag_dict[(record.id, coords[0])] = [loc_number,
                                                              record.id,
                                                              locustag,
                                                              feat.type,
                                                              product_list]
                else:
                    logger.debug("Not adding this feat to " +
                                 "list: %s", product_list)
                    pass
                loc_number = loc_number + 1  # increment index after feature
            else:
                if verbose:
                    logger.debug("skipping: %s", feat.type)
                loc_number = loc_number + 1  # increment index after feature
        # this is a soft warning, as we want to be able to loop
        # through all records before worrying
        if len(locus_tag_dict) < 1 and logger:
            logger.info(str("no locus tags found in {0} for " +
                            "{1} features with annotated products matching " +
                            "{2}!\n").format(record.id,
                                             feature, specific_features))
    locus_tag_dict = locus_tag_dict
    if verbose and logger:
        for key in sorted(locus_tag_dict):
            logger.debug("%s: %s;", key, locus_tag_dict[key])

    # count the occuraces of each feature per genbank record
    nfeatures_occur, \
        nfeat_simple = count_feature_hits(
            all_feature=specific_features is None,
            genome_seq_records=genome_seq_records,
            specific_features=specific_features,
            locus_tag_dict=locus_tag_dict)
    return(locus_tag_dict, nfeatures_occur, nfeat_simple)


def pure_python_kmeans(data, centers=3, kind=int, DEBUG=True):
    """giveb 1d list of numberic data and number of centers, returns a
    csv with the data and cluster, and LP's disapointment
    """
    # handle single instance "clusters" when data is same or less than centers
    if len(data) <= centers:
        # {"index": [coord]}
        return {str(i): [j] for i, j in enumerate(data)}

    with open(os.path.join(os.getcwd(),
                           "pure_python_kmeans_list.csv"), "w") as f:
        for i in data:
            f.write("".join([str(i), "\n"]))
    rcmds = ["# Generated by riboSelect.py on {0}".format(time.asctime()),
             "centers <- {0}".format(centers),
             str("data <- read.csv('pure_python_kmeans_list.csv', " +
                 "header=F, col.names=c('index'))"),
             "set.seed(27)",
             "km <- kmeans(data[,1], nstart=100,iter.max=100,centers=centers)",
             "data[,2] <- km[1]",
             "data <- data[order(data[,1]),  ]",
             str("write.table(data, 'pure_python_kmeans_list.csv', " +
                 "sep=',', row.names=F, col.names=F)")
             ]
    with open(os.path.join(os.getcwd(), "km_script.R"), "w") as f:
        for i in rcmds:
            f.write(i + "\n")
    subprocess.run("Rscript km_script.R", shell=sys.platform != 'win32',
                   check=True)
    with open("pure_python_kmeans_list.csv", mode='r') as infile:
        reader = csv.reader(infile)
        indexClusterDict = {}
        for row in reader:
            try:
                if row[1] in indexClusterDict:
                    indexClusterDict[row[1]].append(row[0])
                else:
                    indexClusterDict[row[1]] = [row[0]]
            except:
                raise ImportError("error constructing dictionary from csv; "
                                  "possibly due to type casting? adjust the "
                                  " 'kind' arg to string if in doubt")
    if not DEBUG:
        os.remove(os.path.join(os.getcwd(), "pure_python_kmeans_list.csv"))
        os.remove(os.path.join(os.getcwd(), "km_script.R"))
    cast_dict = {}
    if kind not in [int, str, float]:
        raise TypeError("can only attempt coersion to int, str, and float; " +
                        "cant coerce to {0}.".format(kind))
    for k, v in indexClusterDict.items():
        cast_dict[k] = [kind(x) for x in v]
    return cast_dict


if __name__ == "__main__":
    args = get_args()
    output_root = os.path.abspath(os.path.expanduser(args.output))
    # Create output directory only if it does not exist
    try:
        os.makedirs(args.output)
    except FileExistsError:
         # leaving comment char'#' added for stream output usage
        print("#Selected output directory %s exists" %
              args.output)
        if not args.clobber:
            print("exiting")
            sys.exit(1)
        else:
            print("# continuing, and risking potential loss of data")
    logger = set_up_logging(
        verbosity=args.verbosity,
        outfile=str("%s_riboSelect_log.txt" %
                    os.path.join(output_root,
                                 time.strftime("%Y%m%d%H%M"))),
        name=__name__)

    # log = sys.stderr.write  # to keep streaming clean if this goes that route
    logger.info("Current usage:\n{0}\n".format(" ".join(sys.argv[1:])))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, v)
    date = str(datetime.datetime.now().strftime('%Y%m%d'))

    # Check if output file exists; if so, remove it
    output_path = os.path.join(args.output,
                               str("riboSelect_grouped_loci.txt"))
    if os.path.exists(output_path):
        if args.clobber:
            logger.info("removing existing output file\n")
            os.remove(output_path)
        else:
            logger.error("Existing output file found!")
            sys.exit(1)

    # get genome records into a list
    with open(args.genbank_genome) as fh:
        genome_records = list(SeqIO.parse(fh, 'genbank'))
    print("Loaded %d records" % len(genome_records))

    # get list of loci matching feature and optionally specific features
    # also returns nfeat, a dict of feature count by genbank id
    logger.debug(
        str(
            "searching {0} for {1} features containing {2} in the " +
            "product annotation").format(
            str([x.id for x in genome_records]),
            args.feature,
            str([x for x in args.specific_features.split(":")])))
    lociDict, nfeat, nfeat_simple = \
        get_filtered_locus_tag_dict(genome_seq_records=genome_records,
                                    feature=args.feature,
                                    specific_features=args.specific_features,
                                    verbose=args.debug,
                                    logger=logger)

    # default case, clusters are inferred
    # if not, must be equal to the length of genbank records
    if args.clusters != "":
        try:
            centers = [int(x) for x in args.clusters.split(":")]
            logger.info(str(centers))
        except:
            logger.error("cannot coerce --clusters to integer after " +
                         "splitting on colons!\n")
            sys.exit(1)
    else:
        centers = [0 for x in genome_records]

    # if unequal lengths, throw error
    # logger.info clusters for accession for user to verify
    if len(genome_records) != len(centers):
        logger.error("centers must be the same length as number" +
                     " of genbank records!\n")
        sys.exit(1)

    #####
    ##### for each genbank record, process, and append any hits to outfile
    #####
    logger.debug("All loci:")
    for k, v in lociDict.items():
        logger.debug(str(k) + "\t" + str(v))
    for i in range(0, len(genome_records)):

        logger.info("Processing {0}\n".format(genome_records[i].id))
        # if user gives clusters, make sure it matches the length:
        if args.clusters:
            logger.info("using %s clusters for %s\n",
                        centers[i], genome_records[i].id)
        # get subset of lociDict for that id
        subset = {key: value for key, value in lociDict.items() if
                  genome_records[i].id in value}
        # skip if that doesnt have any hits
        logger.debug("Subset loci:")
        for k, v in subset.items():
            logger.debug(str(k) + "\t" + str(v))

        if len(subset) == 0:
            logger.info("no hits in {0}\n".format(genome_records[i].id))
            continue
        logger.info("hits in {0}\n".format(genome_records[i].id))

        #  find nfeat for this genbank id by subsetting;
        # is this a bad way of doesnt things?
        logger.debug("centers: {0}".format(centers))
        logger.debug("nfeat_simple: {0}".format(nfeat_simple))
        if nfeat_simple is None and centers[i] == 0:
            logger.error("No specific features submitted, cannot calculate " +
                         " number centers needed for clustering.  Please" +
                         "submit the desired number of clusters with the  " +
                         "--clusters argument!\n")
            sys.exit(1)
        rec_nfeat = list({k: v for k, v in nfeat_simple.items() if
                          genome_records[i].id in k}.values())[0]
        logger.debug("rec_nfeat: {0}".format(rec_nfeat))
        indexes = [x[1] for x in list(subset)]  # get index back from tuple key
        ## if centers[i] is 0, try max and min sequentially; if that fails skip
        if centers[i] == 0:
            if min(rec_nfeat) == 0:
                current_centers = max(rec_nfeat)
            else:
                current_centers = min(rec_nfeat)
            if current_centers == 0:
                logger.info("skipping the clustering for {0}\n".format(i))
                continue
            indexClusters = pure_python_kmeans(indexes,
                                               centers=current_centers,
                                               DEBUG=args.keep_temps)
        else:
            ## if centers[i] is not 0, use it
            current_centers = centers[i]
        # Perform actual clustering
        try:
            # indexClusters should be like { "1": [3,4,6], "2": [66,45,63]}
            indexClusters = pure_python_kmeans(indexes,
                                               centers=current_centers,
                                               DEBUG=args.keep_temps)
        except Exception as e:
            logger.error(e)
            sys.exit(1)

        with open(output_path, "a") as outfile:
            outfile.write("# Generated cluters for {0} on {1}\n".format(
                i, date))
            outfile.write("#$ FEATURE {0}\n".format(args.feature))
            sys.stdout.write("#$ FEATURE {0}\n".format(args.feature))
            sys.stdout.write("# Generated cluters for {0} on {1}\n".format(
                i, date))
            for k, v in indexClusters.items():
                # for each k:v, this replaces the index in v with the locus tag
                # from subset, and writes it out in the way that plays nice
                # riboSeed
                outstr = str(genome_records[i].id + " " +
                             str(":".join([subset[(genome_records[i].id,
                                                   int(x))][2] for x in v])) +
                             '\n')
                outfile.write(outstr)
                # this should be the only thing going to stdout.
                sys.stdout.write(outstr)

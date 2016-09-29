#!/usr/bin/env python
"""
version 0.8.6
Minor version changes:
  - this should work better to filter things by specific feature,
     and to go contig by contig
#TODO:
- set up logging
- Make this less awful
- This would be cool to make this stream friendly
Input:
- genbank file
- location:  rRNA (or RRNA, or something similar)

Output:
-text file with colon separated locus_tags, one for each region
  ECUMN_23S_3:ECUMN_16S_3:ECUMN_5S_4


USAGE:
 $ python riboSnag.py genome.gbk --features rRNA --output /path/to/output
"""
import re
import os
import csv
import subprocess
import datetime
import argparse
import sys
from Bio import SeqIO
import time
from pyutilsnrw.utils3_5 import get_genbank_record, check_single_scaffold


def get_args(DEBUG=False):
    parser = argparse.ArgumentParser(description="This is used to extract" +
                                     " rRNA regions from a gb file, returns" +
                                     "a text file with the clusters")
    parser.add_argument("genbank_genome", help="Genbank file (WITH SEQUENCE)")
    parser.add_argument("-f", "--feature", help="Feature, rRNA or RRNA; " +
                        "default: %(default)s",
                        default='rRNA', type=str)
    parser.add_argument("-s", "--specific_features", help="colon:separated" +
                        " -- specific features\
                         ; default: %(default)s",
                        default='16S:23S:5S', type=str)
    parser.add_argument("-o", "--output", help="output directory;" +
                        "default: %(default)s", default=os.getcwd(),
                        type=str, dest="output")
    parser.add_argument("-c", "--clusters", help="number of rDNA clusters;" +
                        "can be a colon:separated list that matches number " +
                        "of genbank records"
                        "default is inferred: %(default)s", default='',
                        type=str, dest="clusters")
    # parser.add_argument( "--per_contig",
    #                     help="if genome is not small or nearly completed, " +
    #                     "run with --per_contigs. This finds clusters " +
    #                     "on a per-contig  basis as opposed to globally " +
    #                     "searching for loci; default: %(default)s",
    #                     default=False, action="store_true")  # , dest="per_contig")
    args = parser.parse_args()
    return(args)


def multisplit(delimiters, string, maxsplit=0):
    """from SO.
    """
    import re
    regexPattern = '|'.join(map(re.escape, delimiters))
    return re.split(regexPattern, string, maxsplit)


def get_filtered_locus_tag_dict(genome_seq_records, feature="rRNA",
                                specific_features="16S:23S",
                                verbose=True, logger=None):
    """ Given a LIST (as of 20160927) or genbank records,
    returns dictionary of index:locus_tag id pairs for all
    "feature"  entries.  This then gets clustered.
    requires having locus tag in your genbank file.  Non-negitable.
    should be prokka-friendly, so if you have a gb file with legacy
    annotations, just run through prokka (with an rRNA caller)
    20160922 This was changed to add checking for rRNA type from
    ribosome product annoation, not ht locus tag.
    if specific features is None, return all with type == feature
    """
    if not isinstance(genome_seq_records, list):
        raise("Error! this function can only accept a list of records" +
              "simply put your genbank record in brackets if you only " +
              "have a single record")
    if specific_features is None:
        just_feature = True
    else:
        just_feature = False
        specific_features = specific_features.split(":")
    if verbose and logger:
        log_status = logger.info
    elif verbose:
        log_status = sys.stderr.write
    else:
        log_status = print
    locus_tag_dict = {}  # recipient structure
    # loop through records
    for record in genome_seq_records:
        loc_number = 0  # counter
        for feat in record.features:
            try:
                locustag = feat.qualifiers.get("locus_tag")[0]
                product_list = multisplit([",", " ", "-", "_"],
                                          feat.qualifiers.get("product")[0])
                if not just_feature:
                    if feat.type in feature and \
                       any(x in specific_features for x in product_list):
                        locus_tag_dict[loc_number] = [loc_number,
                                                      record.id,
                                                      locustag,
                                                      feat.type,
                                                      product_list]
                    else:
                        pass
                else:
                    if feat.type in feature:
                        locus_tag_dict[loc_number] = [record.id,
                                                      locustag,
                                                      feat.type,
                                                      product_list]
                    else:
                        pass
            except TypeError:
                pass
            loc_number = loc_number + 1  # increment index after each feature
        if len(locus_tag_dict) < 1:
            log_status(str("no locus tags found in {0} for " +
                           "{1} features with annotated products matching" +
                           "{2}!").format(record.id,
                                          feature, specific_features))

    filtered = locus_tag_dict
    if verbose:
        for key in sorted(filtered):
                log_status("%s: %s;" % (key, filtered[key]))

    lociDict = filtered
    nfeatures_occur = {}  # [0 for x in specific_features] n features per gb
    nfeat_simple = {}
    for record in genome_seq_records:
        hit_list = []
        hit_list_simple = []
        for i in specific_features:
            hits = 0
            subset = {k: v for k, v in lociDict.items() if record.id in v}
            for k, v in subset.items():
                if any([i in x for x in v[-1]]):
                    hits = hits + 1
                else:
                    pass
            hit_list.append([i, hits])
            hit_list_simple.append(hits)
        nfeatures_occur[record.id] = (hit_list)
        nfeat_simple[record.id] = hit_list_simple
    return(filtered, nfeatures_occur, nfeat_simple)




def pure_python_kmeans(data, group_by=None, centers=3, DEBUG=True):
    """giveb 1d list of numberic data and number of centers, returns a
    csv with the data and cluster, and LP's disapointment
    """
    with open(os.path.join(os.getcwd(), "list.csv"), "w") as f:
        for i in data:
            f.write("".join([str(i), "\n"]))
    rcmds = ["# Generated by riboSelect.py on {0}".format(time.asctime()),
             "centers <- {0}".format(centers),
             "data <- read.csv('list.csv', header=F, col.names=c('index'))",
             "km <- kmeans(data[,1], nstart=100,iter.max=100,centers=centers)",
             "data[,2] <- km[1]",
             "data <- data[order(data[,1]),  ]",
             "write.csv(data, 'list.csv',  row.names=F)"]
    with open(os.path.join(os.getcwd(), "km_script.R"), "w") as f:
        for i in rcmds:
                f.write(i)
                f.write('\n')
    subprocess.run("Rscript km_script.R", shell=sys.platform != 'win32',
                   check=True)
    with open('list.csv', mode='r') as infile:
        reader = csv.reader(infile)
        next(reader, None)  # skip the headers
        indexClusterDict = dict((rows[0], rows[1]) for rows in reader)

    if not DEBUG:
        os.remove(os.path.join(os.getcwd(), "list.csv"))
        os.remove(os.path.join(os.getcwd(), "km_script.R"))
    return(indexClusterDict)


if __name__ == "__main__":
    args = get_args(DEBUG=False)
    print("Current usage:")
    print(sys.argv[1:])
    # if args.specific_features is not None:
    #     specific_features = args.specific_features.split(":")
    # else:
    #     specific_features = None
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    # if check_single_scaffold(args.genbank_genome):
    #     print("You really should only give this genbank files with one " +
    #           "scaffold for now.")
    #     sys.exit(1)
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    # if args.per_contig:
    output_path = os.path.join(args.output,
                               str(date + "_riboSelect_grouped_loci.txt"))
    if os.path.exists(output_path):
        print("removing existing output file")
        os.remove(output_path)

    # get genome records into a list
    genome_records = get_genbank_record(args.genbank_genome,
                                        first_only=False)
    # get list of loci matching feature and optionally specific features
    # also returns nfeat, a dict of feature count by genbank id
    lociDict, nfeat, nfeat_simple = \
        get_filtered_locus_tag_dict(genome_seq_records=genome_records,
                                    feature=args.feature,
                                    specific_features=args.specific_features,
                                    verbose=True)
    print(lociDict)
    # default case, clusters are inferred
    # if not, must be equal to the length of genbank records
    if args.clusters != "":
        try:
            centers = [int(x) for x in args.clusters.split(":")]
            print(centers)
        except:
            print("cannot coerce --clusters to integer!")
            sys.exit(1)
    else:
        centers = [0 for x in genome_records]
    # if unequal lengths, throw error
    if len(genome_records) != len(centers):
        print("centers must be the same length as number" +
              " of genbank records!")
        sys.exit(1)
    # print clusters for accession for user to verify
    for i in range(0, len(genome_records)):  # for each genbank id
        print("using {0} clusters for {1}".format(
            centers[i], genome_records[i].id))
        # get subset of lociDict for that id
        subset = {key: value for key, value in lociDict.items() if \
                  genome_records[i].id in value }
        # skip if that doesnt have any hits
        if len(subset) == 0:
            print("no hits in {0}".format(genome_records[i].id))
            continue
        #  find nfeat for this genbank id by subsetting;
        # is this a bad way of doesnt things?
    rec_nfeat  = list({k: v for k, v in nfeat_simple.items() if \
                       genome_records[i].id in k }.values())[0]
    if centers[i] == 0:
        indexClusterDict = pure_python_kmeans(lociDict.keys(),
                                              centers=min(rec_nfeat))
    else:
        indexClusterDict =  pure_python_kmeans(lociDict.keys(),
                                               centers=centers[i])
        ## csv to dict
        # with open('list.csv', mode='r') as infile:
        #     reader = csv.reader(infile)
        #     next(reader, None)  # skip the headers
        #     indexClusterDict = dict((rows[0], rows[1]) for rows in reader)
    print(indexClusterDict)
    clusteredDict = {}
    for k, v in indexClusterDict.items():
        clusteredDict.setdefault(v, []).append(
            [x for x in subset[int(k)]])
    with open(os.path.join(args.output,
                           str(date + "_riboSelect_grouped_loci.txt")),
              "a") as outfile:
        for k, v in clusteredDict.items():
            outstr = v[0][1] + " " + ":".join(x[2] for x in v) + '\n'
            outfile.write(outstr)
            sys.stdout.write(outstr)


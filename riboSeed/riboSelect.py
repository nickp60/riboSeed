#!/usr/bin/env python
"""
version 0.8.4
Minor version changes:
  - there was a bug where specific_features were being looked for
    in the locus tagl;  this was a corner case for the coli strain, but
     largly irrelevant others
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
    args = parser.parse_args()
    return(args)


def check_single_scaffold(input_genome_path, logger=None):
    """Test for single scaffold
    """
    counter = -1  # if all goes well, returns 0, else returns 1 or more or -1
    for line in open(input_genome_path, "r"):
        if re.search("ORIGIN", line) is not None:
            counter = counter + 1
    return(counter)


def get_genbank_record(input_genome_path, verbose=True, logger=None):
    """reads the FIRST record only from a genbank file;
    will probably only work for first scaffold
    """
    if verbose and logger:
        log_status = logger.info
    elif verbose:
        log_status = sys.stderr.write
    else:
        pass

    print("Reading genbank file...")
    with open(input_genome_path) as input_genome_handle:
        genome_seq_record = next(SeqIO.parse(input_genome_handle, "genbank"))
    if genome_seq_record.seq[0: 100] == str("N" * 100):
        log_status("Careful: the first 100 nucleotides are N's; " +
                   "did you download the full .gb file?")
    return(genome_seq_record)


def get_filtered_locus_tag_dict(genome_seq_record, feature="rRNA",
                                specific_features="16s:23s:5s",
                                verbose=True, logger=None):
    """returns dictionary of index:locus_tag id pairs for all
    "feature"  entries.  This then gets clustered.
    requires having locus tag in your genbank file.  Non-negitable.
    should be prokka-friendly, so if you have a gb file with legacy
    annotations, just run through prokka (with an rRNA caller)
    20160922 This was changed to add checking for rRNA type from
    ribosome product annoation, not ht locus tag.
    """
    if verbose and logger:
        log_status = logger.info
    elif verbose:
        log_status = sys.stderr.write
    else:
        pass
    loc_number = 0  # counter
    locus_tag_dict = {}  # recipient structure
    for feat in genome_seq_record.features:
        try:
            locustag = feat.qualifiers.get("locus_tag")[0]
            product = feat.qualifiers.get("product")[0]
            locus_tag_dict[loc_number] = [locustag, feat.type, product]
            loc_number = loc_number + 1
        except TypeError:
            pass
    if len(locus_tag_dict) < 1:
            raise ValueError("no locus tags found!")
    if verbose:
        log_status("filtering by feature of interest")
    filtered = {k: v for k, v in locus_tag_dict.items() if v[1] == feature.strip()}
    if len(filtered) == 0:
        log_status("ERROR! no {0} found in locus_tag_dict; rRNA's must have " +
                   "locus tags".format(feature))
        sys.exit(1)
    if verbose:
        for key in sorted(filtered):
                log_status("%s: %s;" % (key, filtered[key]))
    ###
    lociDict = filtered
    nfeatures_occur = []  #  [0 for x in specific_features]
    for i in specific_features:
        hits = 0
        for k, v in lociDict.items():
            if any([i in x for x in v]):
                hits = hits +1
            else:
                pass
        nfeatures_occur.append(hits)
    print(nfeatures_occur)
    for i in range(0, len(specific_features)):
        if nfeatures_occur[i] == 0:
            log_status(str("no features found! check that your file contains" +
                           " {0}; case-sensitive.  rRNA's must have locus " +
                           "tags!").format(specific_features[i]))
            sys.exit(1)
    log_status(str(" occuraces of each specific feature: {0}; using " +
                   " {1} clusters").format(nfeatures_occur,
                                           min(nfeatures_occur)))

    ###
    return(filtered, nfeatures_occur)


def pure_python_kmeans(data, centers=3):
    """giveb 1d list of numberic data and number of centers, returns a
    csv with the data and cluster, and LP's disapointment
    """
    with open(os.path.join(os.getcwd(), "list.csv"), "w") as f:
        for i in data:
            f.write("".join([str(i), "\n"]))
    rcmds = ["# Generated by riboSelect.py on {0}".format(time.asctime()),
             "centers <- {0}".format(centers),
             "data <- read.csv('list.csv', header=F, col.names=c('index'))",
             "km <- kmeans(data[,1], nstart=100, iter.max=100, centers=centers)",
             "data[,2] <- km[1]",
             "data <- data[order(data[,1]),  ]",
             "write.csv(data, 'list.csv',  row.names=F)"]
    with open(os.path.join(os.getcwd(), "km_script.R"), "w") as f:
        for i in rcmds:
                f.write(i)
                f.write('\n')
    subprocess.run("Rscript km_script.R", shell=sys.platform != 'win32',
                   check=True)


if __name__ == "__main__":
    args = get_args(DEBUG=False)
    print("Current usage:")
    print(sys.argv[1:])
    specific_features = args.specific_features.split(":")
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    if check_single_scaffold(args.genbank_genome):
        print("You really should only give this genbank files with one " +
              "scaffold for now.")
        sys.exit(1)
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    genome_record = get_genbank_record(args.genbank_genome)
    lociDict, nfeat = \
        get_filtered_locus_tag_dict(genome_seq_record=genome_record,
                                    feature=args.feature,
                                    specific_features=specific_features,
                                    verbose=True)
    print(lociDict)
    # for each specific feature, get number of occurances
    pure_python_kmeans(lociDict.keys(), centers=min(nfeat))
    ## csv to dict
    with open('list.csv', mode='r') as infile:
        reader = csv.reader(infile)
        next(reader, None)  # skip the headers
        indexClusterDict = dict((rows[0], rows[1]) for rows in reader)
    ##
    # print(lociDict)
    # print(indexClusterDict)
    clusteredDict = {}
    for k, v in indexClusterDict.items():
        clusteredDict.setdefault(v, []).append(lociDict[int(k)])
    # print(clusteredDict)
    with open(os.path.join(args.output,
                           str(date +
                               "_riboSelect_grouped_loci.txt")),
              "w") as outfile:
        for k, v in clusteredDict.items():
            outfile.write(str(":".join(str(x[0]) for x in v) + '\n'))
            sys.stdout.write(str("\t".join(str(x[0]) for x in v) + '\n'))

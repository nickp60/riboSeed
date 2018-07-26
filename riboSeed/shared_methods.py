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
import logging
import sys
import re
import os
import glob
import multiprocessing
import subprocess
import hashlib

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from .classes import Locus, LociCluster

logger = logging.getLogger('root')

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


def set_up_logging(verbosity, outfile, name):
    """
    Set up logging a la pyani, with
    a little help from:
    https://aykutakin.wordpress.com/2013/08/06/
        logging-to-console-and-file-in-python/
    requires logging, os, sys, time
    and some coloring
    logs debug level to file, and [verbosity] level to stderr
    return a logger object
    """
    logger = logging.getLogger('root')
    if (verbosity * 10) not in range(10, 60, 10):
        raise ValueError('Invalid log level: %s' % verbosity)
    # logging.basicConfig(level=logging.DEBUG)
    logger.setLevel(logging.DEBUG)
    logging.addLevelName(logging.DEBUG,  "DBUG")
    # create console handler and set level to given verbosity
    console_err = logging.StreamHandler(sys.stderr)
    console_err.setLevel(level=(verbosity * 10))
    console_err_format = logging.Formatter(
        str("%(asctime)s " + "%(levelname).4s" +" %(message)s"),
        "%H:%M:%S")
    console_err.setFormatter(console_err_format)
    # set some pretty colors, shorten names of loggers to keep lines aligned
    # sadly, this doesn't jive with log files (casue they cant display the col
    # logging.addLevelName(logging.DEBUG, "\u001b[30m%s\033[1;0m" % "..")
    # logging.addLevelName(logging.INFO,  "\u001b[32m%s\033[1;0m" % "--")
    # logging.addLevelName(logging.WARNING, "\u001b[33m%s\033[1;0m" % "!!")
    # logging.addLevelName(logging.ERROR, "\u001b[31m%s\033[1;0m" % "xx")
    # logging.addLevelName(logging.CRITICAL, "\u001b[31m%s\033[1;0m" % "XX")
    logger.addHandler(console_err)
    # create debug file handler and set level to debug
    try:
        logfile_handler = logging.FileHandler(outfile, "w")
        logfile_handler.setLevel(logging.DEBUG)
        logfile_handler_formatter = \
            logging.Formatter("%(asctime)s - %(levelname).4s - %(message)s")
        logfile_handler.setFormatter(logfile_handler_formatter)
        logger.addHandler(logfile_handler)
    except:
        logger.error("Could not open {0} for logging".format(outfile))
        sys.exit(1)
    logger.debug("Initializing logger")
    logger.debug("logging at level {0}".format(verbosity))
    return logger


def make_barrnap_cmd(infasta, outgff, exe, thresh, kingdom, evalue=1e-06, threads=1):
    """ construct a system call to barrnap

    Args:
        infasta (str): path to fasta file
        outgff (str): name of gff file to create
        exe (str): path to system executable for barrnap
        thresh (float): Proportional length threshold to reject prediction
        kingdom (str): bac, euk, arc, or mito
        threads (int): number of threads to use
    Returns:
        (str): command

    """
    assert thresh > 0 and thresh < 1, "Thresh must be between 0 and 1!"
    if exe.endswith("py"):
        # ensure running python barrnap uses >3.5
        pyexe = str(sys.executable + " ")
    else:
        pyexe = ""
    cmd = "{0}{1} -k {2} {3} --reject {4} --threads {5} --evalue {6} > {7}".format(
        pyexe,
        exe,
        kingdom,
        infasta,
        thresh,
        threads,
        evalue,
        outgff
        )
    return cmd


def test_barrnap_ok(exe):
    """Ensure barrnap isnt missing perl dependencies if not using py barrnap

    Args:
        exe (str): path to system executable for barrnap
    Raises:
        SystemExit: if we get an error when running barrnap
    Returns:

    """
    cmd = "{0} -h".format(exe)
    result = subprocess.run(cmd,
                            shell=sys.platform != "win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=False)
    if result.returncode == 127:
        print("It appears that you have barrnap installed, but it lacks " +
              "certain perl dependecies. Either install those, or install " +
              "python barrnap from https://github.com/nickp60/barrnap-python")
        sys.exit(1)


def check_version_from_cmd(
        exe,
        cmd, line,
        pattern=r"^__version__ = '(?P<version>[^']+)'$",
        where='stderr',
        min_version="0.0.0", logger=None,
        coerce_two_digit=False):
    """the guts have been stolen from pyani; returns version
    from an system call that should return a version string.
    Hacky, but better than nothing.
    line arg is 1-indexed
    .strip() is called on match to remove whitspaces
    20170920 changed to remove shutil.which call.
    That should be done outside of this funciton
    """
    assert logger is not None, "must use logging"
    from distutils.version import StrictVersion
    result = subprocess.run("{0} {1}".format(exe, cmd),
                             # is this a securiy risk?
                            shell=sys.platform != "win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=False)
    logger.debug(result)
    try:
        if where == 'stderr':
            printout = result.stderr.decode("utf-8").split("\n")
        elif where == 'stdout':
            printout = result.stdout.decode("utf-8").split("\n")
        else:
            raise ValueError("where option can only be 'stderr' or 'stdout'")
    except Exception as e:
        raise e
    logger.debug(printout)

    this_version = None
    try:
        m = re.search(pattern, printout[line - 1])
    except IndexError as e:
        raise e
    if m:
        this_version = m.group('version').strip()
    logger.debug("this_version: %s", this_version)
    if coerce_two_digit:
        this_version = "0.{0}".format(this_version)
        logger.debug("coerced this_version: %s", this_version)
    if this_version is None:
        raise ValueError("No verison was captured with pattern" +
                         "{0}".format(pattern))
    try:
        if StrictVersion(this_version) < StrictVersion(min_version):
            raise ValueError("{0} version {1} must be greater than {2}".format(
                cmd, this_version, min_version))
    except Exception as e:
        raise e
    return(this_version)


def combine_contigs(contigs_dir, pattern="*",
                    contigs_name="riboSeedContigs_aggregated",
                    ext=".fasta", verbose=False, logger=None):
    """changed over to biopython
    combine all *ext files in dir, return path to concatenated file
    requires Bio.SeqIO, glob, os
    """
    # make sure it ends in a path sep
    if contigs_dir[-1] != os.path.sep:
        contigs_dir = str(contigs_dir + os.path.sep)
    output = os.path.join(contigs_dir, str(contigs_name + ext))
    fastas = glob.glob(str(contigs_dir + pattern + ext))
    fastas.sort()
    if verbose:
        print(str("combining the following files matching pattern " +
                  "{0}:{1}".format(pattern, " ".join(fastas))))
    if logger:
        logger.debug(str("combining the following files matching pattern " +
                        "{0}:{1}".format(pattern, " ".join(fastas))))
    if len(fastas) == 0:
        if logger:
            logger.error("No matching files to combine found" +
                         " in {0}!".format(contigs_dir))
        sys.exit(1)
    with open(output, 'w') as w_file:
        for filen in fastas:
            with open(filen, 'r') as o_file:
                seq_records = SeqIO.parse(o_file, 'fasta')
                SeqIO.write(seq_records, w_file, 'fasta')

    return output


def get_number_mapped(bam, samtools_exe, logger=None):
    """use samtools flagstats to retrieve total mapped reads as a diagnostic
    returns a string to be printed, the 4th line of flagstat
    requires subprocess, sys, logger (SAMtools)
    """
    flagstatcmd = str("{0} flagstat {1}").format(samtools_exe, bam)
    flagstats = subprocess.run(flagstatcmd,
                               shell=sys.platform != "win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
    try:
        printout = flagstats.stdout.decode("utf-8").split("\n")[4]
    # TODO test for none mapped
    except IndexError:
        if logger:
            logger.error("Error reading {0}".format(bam))
        sys.exit(1)
    return printout


def keep_only_first_contig(ref, newname="contig1"):
    # TODO make a biopython version
    """
    given a multi fasta from SPAdes, extract first entry,
    rename "NODE_1" with newname, overwrite file
    requires os, re
    fasta_sequences = SeqIO.parse(open(ref),'fasta')
    with open(output_file) as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, fasta.seq.tostring()
            new_sequence = some_function(sequence)
            write_fasta(out_file)

    """
    temp = os.path.join(os.path.dirname(ref), "temp.fasta")
    with open(temp, "w") as outfile:
        with open(ref, "r") as file_handle:
            lines = file_handle.readlines()
            if lines[0][0] != ">":
                raise ValueError(str("Error  with spades output contig!" +
                                     " Not a valid fasta!"))
            new_header = str(re.sub("NODE_\d*", newname, str(lines[0])))
            outfile.write(new_header)
            for line in range(1, len(lines)):
                if lines[line][0] == ">":
                    break  # stop after first entry
                outfile.write(str(lines[line]))
    os.remove(ref)
    os.rename(temp, ref)


def get_fasta_lengths(fasta):
    """ given a fasta file, return list of sequence lengths as list
    """
    len_list = []
    with open(fasta, "rt") as data:
        for read in SeqIO.parse(data, "fasta"):
            len_list.append(len(read))
    return len_list


def file_len(fname):
    """http://stackoverflow.com/questions/845058/
    how-to-get-line-count-cheaply-in-python
    """
    if os.path.splitext(fname)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    with open_fun(fname) as f:
        for i, l in enumerate(f):
            index = i
            pass
    return index + 1


def multisplit(delimiters, string, maxsplit=0):
    """from SO. takes a list of delimiters and a string, and
    returns a list of split string
    """
    import re
    assert isinstance(delimiters, list)
    regexPattern = '|'.join(map(re.escape, delimiters))
    return re.split(regexPattern, string, maxsplit)


def md5(fname, string=False):
    """straight from quantumsoup
    http://stackoverflow.com/questions/3431825/
        generating-an-md5-checksum-of-a-file
    updated 20160916 for strings; if true, fname can be a string
    """
    if string:
        return hashlib.md5(fname.encode('utf-8')).hexdigest()
    else:
        hash_md5 = hashlib.md5()
        with open(fname, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
            return hash_md5.hexdigest()

def get_genbank_record(input_genome_path, logger=None, first_only=True):
    """Returns a list of records in a GenBank file"""
    if logger:
        logger.info("Reading genbank file...")
    with open(input_genome_path) as input_genome_handle:
        records = list(SeqIO.parse(input_genome_handle, "genbank"))
    if first_only:
        return [records[0]]
    else:
        return records

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
import glob
import re
import shutil
import subprocess
import os
import traceback

from .shared_methods import set_up_logging, make_barrnap_cmd, \
    test_barrnap_ok
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_args(test_args=None):  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(prog="ribo scan",
        description="Given a directory of one or more chromosomes as fasta " +
        "files, this facilitates reannotation of rDNA regions with Barrnap " +
        " and outputs all sequences as a single, annotated genbank file",
        # usage='xuxxxxxu',
        add_help=False)  # to allow for custom help
    parser.prog = "ribo scan"
    parser.add_argument("contigs", action="store",
                        help="either a (multi)fasta or a directory " +
                        "containing one or more chromosomal " +
                        "sequences in fasta format")

    # taking a hint from http://stackoverflow.com/questions/24180527
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-k", "--kingdom", dest='kingdom',
                          action="store",
                          choices=["bac", "euk", "arc", "mito"],
                          help="whether to look for eukaryotic, archaeal, or" +
                          " bacterial rDNA; " +
                          "default: %(default)s", default="bac",
                          type=str)
    optional.add_argument("-t", "--id_thresh", dest='id_thresh',
                          action="store", type=float,
                          help="partial rRNA hits below this threshold will " +
                          "be ignored. " +
                          "default: %(default)s", default=0.5,
                          )
    optional.add_argument("-b", "--barrnap_exe", dest='barrnap_exe',
                          action="store",
                          help="path to barrnap executable; " +
                          "default: %(default)s", default="barrnap",
                          type=str)
    optional.add_argument("-n", "--name", dest='name',
                          action="store",
                          # default="contig",
                          help="name to give the contig files; "
                          "default: infered from file",
                          type=str)
    optional.add_argument("-c", "--cores", dest='cores',
                          action="store",
                          help="number of threads/cores to use; " +
                          "default: %(default)s", default=2,
                          type=int)
    optional.add_argument("-s", "--seqret_exe", dest='seqret_exe',
                          action="store",
                          help="path to seqret executable, usually " +
                          "installed with emboss; " +
                          "default: %(default)s", default="seqret",
                          type=str)
    optional.add_argument("-m", "--min_length", dest='min_length',
                          action="store",
                          help="skip the scaffold if its shorter than this " +
                          "default: %(default)s", default=0,
                          type=int)
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    if test_args is None:
        args = parser.parse_args(sys.argv[2:])
    else:
        args = parser.parse_args(test_args)
    return args


def last_exception():
    """ Return last exception as a string, or use in logging.
    stolen verbatim from pyani
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


def parse_fasta_header(first_line):
    """ return accession from first line of a fasta file as a string

    Args:
        first_line (str): The first line of a file
    Raise:
        ValueError: if first_line doesnt start with a ">"
        ValueError: if regex to extract accession fails to get a group
    Returns:
        str: name of the accession

    """
    if not first_line.startswith(">"):
        raise ValueError("A valid fasta must start with a '>'!")
    if first_line.startswith(">gi"):
        accession = re.search(r">gi\|.*?\|.*?\|(.*?)\|",
                              first_line).groups()[0]
    else:
        try:  # this handles normal flavored fastas with '>seqid descr'
            accession = re.search(r">(.*?) .*", first_line).groups()[0]
        except AttributeError:  # this handles headers with only '>seqid\n'
            accession = first_line.strip().replace(">", "")
    if accession is None:
        raise ValueError(
            str("unable to extract accession from first line in file: \n" +
                first_line))
    return accession


def add_locus_tags_to_gff(gff, acc):
    """ Add fake locus tags to gff file in place.

    Args:
        gff (str): path to gff file
        acc (str): accession to label the loci with
    Returns:
        retrun path to new gff file (the old one, but with _tagged at the end

    """
    LOCUS = 0
    gff_list = []
    with open(gff, 'r') as g:
        for idx, line in enumerate(g):
            line = line.strip()
            if idx == 0:
                gff_list.append(line)
            else:
                gff_list.append(
                    re.sub("Name=",
                           "locus_tag={0}_{1};Name=".format(acc, LOCUS),
                           line))
                LOCUS = LOCUS + 1

    new_gff = str(os.path.splitext(gff)[0] + "_tagged.gff")
    with open(new_gff, 'w') as outgff:
        for l in gff_list:
            # print(str(l))
            outgff.write(str(l + "\n"))
    return(new_gff)


def combine_gbs(finalgb, gb_list):
    """ pretty self explanatory - it combines files

    it skips lines starting with "##"

    Args:
        finalgb (str): path to ourput file
        gb_list (str): list of other file paths to be combined
    Returns:
        None: doesn't return something

    """
    with open(finalgb, 'w') as outfile:
        for idx, fname in enumerate(gb_list):
            with open(fname) as infile:
                for line in infile:
                    if line.startswith("##") and idx != 0:
                        continue
                    outfile.write(line)


def append_accession_and_version(accession, ingb, finalgb):
    """ add accession, version and GI to a gb file

    Biopython requires a version and an accession to parse correctly.
    This adds it

    Args:
        accession (str): Accession string to add
        ingb (str): path to initial genbank file
        finalgb (str): path to resulting genbank file
    Returns:
        None
    Raises:
        None

    """
    with open(finalgb, 'w') as outfile:
        with open(ingb) as infile:
            for idx, line in enumerate(infile):
                if idx == 1:
                    repline = "ACCESSION   {0}\nVERSION     {0}{1}{2}".format(
                        accession,
                        "  GI:0000000\n",
                        line)
                    outfile.write(repline)
                else:
                    outfile.write(line)


def make_seqret_cmd(exe, outgb, ingff, infasta):
    """construct system call to seqret to make genbank file

    Construct a genbank file from a gff and fasta file

    Args:
        exe (str): path to seqret executable or just the name (seqret)
        outgb (str): path to resulting genbank file
        ingff (str): path to gff(3?) file from barrnap
        infasta (str): path to fasta file
    Returns:
        (str): command
    Raises:
        None

    """
    assert shutil.which(exe) is not None, "seqret executable not found!"
    cmd = str(
        "{0} -sequence {1} -feature -fformat gff3 -fopenfile {2} " +
        "-osformat genbank -auto -outseq {3}"
    ).format(
        shutil.which(exe),
        infasta,
        ingff,
        outgb)
    return cmd


def checkSingleFasta(fasta, logger):
    """raise systen exit if multiple entries in fasta.

    This is mostly depreciated, but left in place juuuust in case

    Args:
        fasta (type): path to file
    Returns:
        None
    Raises:
        SystemExit: thrown if fasta has multiple entries

    """
    with open(fasta, 'r') as f:
        counter = 0
        for rec in SeqIO.parse(f, "fasta"):
            counter = counter + 1
            if counter > 1:
                logger.error(
                    "Fasta has multiple entries! " +
                    "use splitMultifasta.sh to create single-entry " +
                    "fastas from it. exiting...")
                sys.exit(1)


def getFastas(inp, output_root, name, logger):
    """return list of fasta files

    If the input is a multifasta, a new directory is created to contain
    single sequence files.  given an input that could either
    be a single fasta, multifasta,
    or a directory, retrun the appropriate files

    Args:return list of fasta files
        inp (str): either the path to a file or a directory
        output_root: path to where the contigs dir will be created, if needed
        name (str): accession for resulting contig files or sequences
    Returns:
        (list): paths to sequences
    Raises:
        SystemExit: if file is empty or path doesnt exist
        SystemExit: if no files end in *fasta or *fa

    """
    if not os.path.isdir(os.path.expanduser(inp)):
        if not os.path.isfile(os.path.expanduser(inp)):
            logger.error("'%s' is not a valid directory or file!",
                         os.path.expanduser(inp))
            sys.exit(1)
        else:
            if not os.path.getsize(os.path.expanduser(inp)) > 0:
                logger.error("'%s' looks like an empty file!",
                             os.path.expanduser(inp))
                sys.exit(1)
            logger.info("spliting multifasta into multiple fastas " +
                        "for easier processing")
            splitMultifasta(multi=inp, output=output_root, name=name,
                            logger=logger)

            fastas = glob.glob(os.path.join(output_root, "contigs",
                                            "*.fa"))
    else:
        fastas = glob.glob(os.path.join(os.path.expanduser(inp),
                                        "*.fa"))
    if len(fastas) == 0:
        logger.info("No fasta files found with extention '.fa'. " +
                       "searching for '*.fasta' files instead...")
        fastas = glob.glob(os.path.join(os.path.expanduser(inp),
                                        "*.fasta"))
        if len(fastas) == 0:  # still
            logger.error("No files in %s with extention '*fasta'! Exiting",
                         inp)
            sys.exit(1)
    # unify the output of multifasta
    combine_gbs(finalgb=os.path.join(output_root, "scannedScaffolds.fa"),
                gb_list=fastas)
    return(fastas)


def splitMultifasta(multi, output, name, dirname="contigs", logger=None):
    """create new files containing a single fasta entry each

    regex stolen from SO
    name is the name of the file, output is the parent dir for the output dir

    Args:
        multi (str): path to multifasta
        output (str): path to output directory
        name (str): accession for renaming multifasta
        dirname (str): name of directory to fold split fastas
    Returns:
        None
    Raises:
        None

    """
    idlist = []
    assert logger is not None, "must use logging"
    os.makedirs(os.path.join(output, dirname))
    source_fname = os.path.splitext(os.path.basename(multi))[0]
    with open(os.path.expanduser(multi), "r") as mf:
        for idx, rec in enumerate(SeqIO.parse(mf, "fasta")):
            logger.debug("record %d: %s", idx, rec.id)
            # logger.debug(rec)
            if name is not None:
                fname = name
            else:
                fname = rec.id
            # if fname in idlist:
            fname = fname + "_" + str(idx)
            if not re.match("^[a-zA-Z0-9_\.]*$", fname):
                logger.warning("Problem with header %s", fname)
                logger.warning("file header contains special characters! " +
                               "Valid characters are in set [a-zA-Z0-9_]; " +
                               "rename with the --name option to prevent " +
                               "later issues. Renaming as source filename " +
                               "with index, sorr..." )
                fname = source_fname + "_" + str(idx)
            with open(os.path.join(output, dirname,
                                   fname + ".fa"), "w") as outf:
                renamed_rec = SeqRecord(rec.seq, id=fname,
                                        description="from riboScan")
                SeqIO.write(renamed_rec, outf, "fasta")
            idlist.append(fname)
    logger.debug("rewrote the following records")
    logger.debug("\n".join(idlist))


def main(args, logger=None):
    # allow user to give relative paths
    output_root = os.path.abspath(os.path.expanduser(args.output))
    sys_barrnap = shutil.which(args.barrnap_exe)
    assert sys_barrnap is not None, \
        "barrnap executable not found!"
    test_barrnap_ok(exe=sys_barrnap)
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory already exists; exiting...")
        sys.exit(1)
    t0 = time.time()
    if logger is None:
        logger = set_up_logging(
            verbosity=args.verbosity,
            outfile=os.path.join(output_root, "riboScan.log"),
            name=__name__)
    logger.info("Usage:\n%s\n", " ".join([x for x in sys.argv]))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, str(v))

    output_file = os.path.join(output_root, "scannedScaffolds.gb")
    output_file_gff = os.path.join(output_root, "scannedScaffolds.gff")
    ##  get and check list of input files
    fastas = getFastas(inp=args.contigs, output_root=output_root,
                       name=args.name, logger=logger)
    gb_list = []
    gff_list = []  # for tagged gffs, that is
    for fasta in sorted(fastas):
        # check for multiple entry fastas.
        # This can still happen if mutlifastas are in dir input
        checkSingleFasta(fasta, logger=logger)
        with open(fasta, 'r') as f:
            rec = next(SeqIO.parse(f, "fasta"))
        accession = rec.id
        if args.min_length > len(rec.seq):
            logger.debug("Skipping %s: must be %d bp, set by --min_length",
                         accession, args.min_length)
            continue
        logger.info("Accession for %s: %s", fasta, accession)
        barrnap_cmd = make_barrnap_cmd(
            infasta=fasta,
            outgff=os.path.join(output_root, "{0}.gff".format(accession)),
            exe=sys_barrnap,
            threads=args.cores,
            thresh=args.id_thresh,
            kingdom=args.kingdom)
        logger.info("running barrnap cmd: %s", barrnap_cmd)
        subprocess.run(barrnap_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        logger.info("grooming gff")
        tagged_gff = add_locus_tags_to_gff(
            gff=os.path.join(output_root, "{0}.gff".format(accession)),
            acc=accession)
        unfinished_gb = os.path.join(
            output_root, "{0}_pre.gb".format(accession))
        seqret_cmd = make_seqret_cmd(
            exe=args.seqret_exe,
            outgb=os.path.join(output_root, "{0}_pre.gb".format(accession)),
            ingff=tagged_gff, infasta=fasta)
        subprocess.run(seqret_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        append_accession_and_version(
            accession=accession,
            ingb=unfinished_gb,
            finalgb=os.path.join(output_root, "{0}.gb".format(accession)))
        gb_list.append(os.path.join(output_root, "{0}.gb".format(accession)))
        gff_list.append(tagged_gff)
    ### no more fastas
    combine_gbs(finalgb=output_file, gb_list=gb_list)
    # this is why you give things intellegent method names
    combine_gbs(finalgb=output_file_gff, gb_list=gff_list)
    # Report that we've finished
    logger.info("Done: %s", time.asctime())
    logger.info("combined scaffolds can be found here: %s", output_file)
    logger.info("Time taken: %.3fs" % (time.time() - t0))

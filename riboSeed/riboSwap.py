#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
"""

import argparse
import sys
import time
import logging
import os
import traceback

from Bio import SeqIO

from .shared_methods import set_up_logging

# --------------------------- methods --------------------------- #


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        prog="ribo swap",
        description="Given de novo and de fere novo contigs files, a " +
        "misjoined de fere novo contig name,  and a colon:separated " +
        "list of de novo contig names, replace the offending contig with " +
        "the de novo contig(s) ",
        add_help=False)  # to allow for custom help
    parser.prog = "ribo swap"
    parser.add_argument("de_novo_file", action="store",
                        help="multifasta containing de novo contigs")

    parser.add_argument("de_fere_novo_file", action="store",
                        help="multifasta containing de fere novo contigs")

    parser.add_argument("bad_contig", action="store",
                        help="name of the bad contig")
    parser.add_argument("good_contigs", action="store",
                        help="colon separated good contigs for replacement")

    # taking a hint from http://stackoverflow.com/questions/24180527
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    # had to make this faux "optional" parse so that the named required ones
    # above get listed first
    optional = parser.add_argument_group('optional arguments')
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
    args = parser.parse_args(sys.argv[2:])
    return args


def last_exception():
    """ Returns last exception as a string, or use in logging.
    stolen verbatim from pyani
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


def remove_bad_contig(infile, outfile, bad_name,
                      logger=None):
    """uruk-hai voice: "Find the bad contig! Gah! Find the bad contig!"
    loop through the input file, and write out if the bad name
    isnt in the record id.
    return a path to the output file, or raise an error if duplicates found
    """
    FOUND = False
    with open(infile, "r") as incontigs:
        with open(outfile, "a") as outcontigs:
            for i in SeqIO.parse(incontigs, "fasta"):
                if bad_name in i.id:
                    if FOUND:
                        logger.error("mutiple contigs with that name found!")
                        raise ValueError
                    else:
                        logger.info("found bad contig: %s; removing", i.id)
                        FOUND = True
                else:
                    SeqIO.write(i, outcontigs, "fasta")
    if not FOUND:
        logger.warning("No contig found matching the supplied name: %s",
                       bad_name)
    return(outfile)


def append_replacement_contigs(infile, outfile, name_list, logger=None):
    """
    loop through the input file (de novo), and write out to the de fere novo
    fasta without the bad contig. if namelist isnt in the record id.
    return a path to the output file, or raise an error if duplicates found
    """
    found_list = []  # this is used to check for dups late on

    with open(infile, "r") as de_novo:
        with open(outfile, "a") as clean_de_fere:
            for i in SeqIO.parse(de_novo, "fasta"):
                for name in name_list:
                    if name not in i.id:
                        pass
                    else:
                        if name in found_list:
                            logger.error("mutiple contigs found with name %s!",
                                         name)
                            raise ValueError
                        else:
                            logger.info("appending %s to output file", i.id)
                            i.id = str("SWAPPED_" + i.id)
                            i.name, i.description = '', ''  # just change id
                            SeqIO.write(i, clean_de_fere, "fasta")
                        found_list.append(name)
            clean_de_fere.write("\n")
    if len(found_list) != len(name_list):
        logger.error("Not all replacement contigs found!")
        raise ValueError
    return(outfile)


def main(args, logger=None):
    # allow user to give relative paths
    output_root = os.path.abspath(os.path.expanduser(args.output))
    try:
        os.makedirs(output_root, exist_ok=False)
    except OSError:
        print("Output directory already exists; exiting...")
        sys.exit(1)
    t0 = time.time()
    log_path = os.path.join(output_root, "riboSwap.log")
    if logger is None:
        logger = set_up_logging(verbosity=args.verbosity,
                                outfile=log_path,
                                name=__name__)
    # # log version of riboSeed, commandline options, and all settings

    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))

    output_file = os.path.join(
        output_root,
        "{0}_swapped.fasta".format(
            os.path.splitext(os.path.basename(args.de_fere_novo_file))[0]))
    if os.path.exists(output_file):
        raise OSError("existing output file")
    rep_names = args.good_contigs.strip().split(":")
    ###
    try:
        outf = remove_bad_contig(infile=args.de_fere_novo_file,
                                 outfile=output_file,
                                 bad_name=args.bad_contig,
                                 logger=logger)
    except:
        logger.error("something bad happened when trying to remove the bad " +
                     "contig.  Something real bad.")
        sys.exit(1)

    try:
        append_replacement_contigs(infile=args.de_novo_file, outfile=outf,
                                   name_list=rep_names, logger=logger)
    except:
        logger.error("something bad happened when trying to append the good " +
                     "contig(s).  Something real bad.")
        sys.exit(1)
    ###

    # Report that we've finished
    logger.info("Done: %s", time.asctime())
    logger.info("file with bad contig replace can be found here: %s", outf)
    logger.info("Time taken: %.3fs" % (time.time() - t0))

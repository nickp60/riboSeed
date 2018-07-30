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
import os
import sys
import subprocess
import argparse
import multiprocessing
import glob

from Bio import SeqIO
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from .shared_methods import set_up_logging, combine_contigs


def get_args(test_args=None):  # pragma: no cover
    parser = argparse.ArgumentParser(prog="ribo score",
        description="This does some simple blasting to detect correctness " +
        "of riboSeed results")
    parser.prog = "ribo score"
    parser.add_argument("indir",
                        help="dir containing a genbank file, assembly files" +
                        "as fastas. Usually the 'mauve' dir in the riboSeed " +
                        "results")
    parser.add_argument("-o", "--output", dest='output',
                        help="directory in which to place the output files",
                        default=None)
    parser.add_argument("-l", "--flanking_length",
                        help="length of flanking regions, in bp; " +
                        "default: %(default)s",
                        default=1000, type=int, dest="flanking")
    parser.add_argument("-p", "--min_percent", dest="min_percent",
                        help="minimum percent identity",
                        default=97, type=int)
    parser.add_argument("-f", "--assembly_ext", dest="assembly_ext",
                        help="extenssion of reference, usually fasta",
                        default="fasta", type=str)
    parser.add_argument("-g", "--ref_ext", dest="ref_ext",
                        help="extension of reference, usually .gb",
                        default="gb", type=str)
    parser.add_argument("-F", "--blast_Full", dest="blast_full",
                        help="if true, blast full sequences along with " +
                        "just the flanking. Interpretation is not " +
                        "implemented currently as false positives cant " +
                        "be detected this way",
                        default=False, action="store_true")
    parser.add_argument("-v", "--verbosity", dest='verbosity',
                        action="store",
                        default=2, type=int, choices=[1, 2, 3, 4, 5],
                        help="Logger writes debug to file in output dir; " +
                        "this sets verbosity level sent to stderr. " +
                        " 1 = debug(), 2 = info(), 3 = warning(), " +
                        "4 = error() and 5 = critical(); " +
                        "default: %(default)s")
    # parser.add_argument("-t", "--blast_type",
    #                     help="blastn or tblastx", default="tblastx")
    if test_args is None:
        args = parser.parse_args(sys.argv[2:])
    else:
        args = parser.parse_args(test_args)
    return(args)


def make_nuc_nuc_recip_blast_cmds(
        query_list, output, subject_file=None, logger=None):
    """given a file, make a blast cmd, and return path to output csv
    """
    assert logger is not None, "must use logging"
    blast_cmds = []
    blast_outputs = []
    recip_blast_outputs = []
    for f in query_list:
        # run forward, nuc aganst prot, blast
        output_path_tab = str(
            os.path.join(output,
                         os.path.splitext(os.path.basename(f))[0] +
                         "_vs_ref.tab"))
        blast_cline = NcbiblastnCommandline(query=f,
                                            subject=subject_file,
                                            # evalue=.001,
                                            outfmt=6, out=output_path_tab)
        add_params = str(" -num_threads 1 -num_alignments 50")
        blast_command = str(str(blast_cline) + add_params)
        blast_cmds.append(blast_command)
        blast_outputs.append(output_path_tab)
        # run reverse, prot against nuc, blast
        recip_output_path_tab = os.path.join(
            output,
            "ref_vs_" + os.path.splitext(os.path.basename(f))[0] + ".tab")
        recip_blast_cline = NcbiblastnCommandline(
            query=subject_file,
            subject=f,
            # evalue=.001,
            outfmt=6, out=recip_output_path_tab)
        recip_blast_command = str(str(recip_blast_cline) + add_params)
        blast_cmds.append(recip_blast_command)
        recip_blast_outputs.append(recip_output_path_tab)

    return(blast_cmds, blast_outputs, recip_blast_outputs)


def merge_outfiles(filelist, outfile):
    """
    """
    # only grab .tab files, ie, the blast output
    filelist = [i for i in filelist if i.split(".")[-1:] == ['tab']]
    if len(filelist) == 1:
        # print("only one file found! no merging needed")
        return(filelist)
    else:
        # print("merging all the blast results to %s" % outfile)
        nfiles = len(filelist)
        fout = open(outfile, "a")
        # first file:
        with open(filelist[0]) as firstf:
            for line in firstf:
                fout.write(line)
        #  now the rest:
        for num in range(1, nfiles):
            with open(filelist[num]) as otherf:
                for line in otherf:
                    fout.write(line)
        fout.close()
    return(outfile)


def BLAST_tab_to_df(path):
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score"]
    with open(path) as tab:
        raw_csv_results = pd.read_csv(
            tab, comment="#", sep="\t", names=colnames)
    return raw_csv_results


def filter_recip_BLAST_df(df1, df2, min_percent, min_lens, logger=None):
    """ results from pd.read_csv with default BLAST output 6 columns

    returns a df
    """
    assert logger is not None, "must use a logger"
    logger.debug("shape of blast results")
    logger.debug("shape of recip blast results")
    # df1['genome'] = df1.query_id.str.split('_').str.get(0)
    # df2['genome'] = df2.subject_id.str.split('_').str.get(0)
    df1['genome'] = df1.query_id
    df2['genome'] = df2.subject_id
    logger.debug(df1.shape)
    logger.debug(df2.shape)
    # recip structure
    filtered = pd.DataFrame(columns=df1.columns)
    unq_subject = df1.subject_id.unique()
    unq_query = df1.genome.unique()
    recip_hits = []
    nonrecip_hits = []
    for gene in unq_subject:
        for genome in unq_query:
            logger.debug("Checking %s in %s for reciprocity" % (gene, genome))
            tempdf1 = df1.loc[(df1["subject_id"] == gene) &
                              (df1["genome"] == genome), ]
            tempdf2 = df2.loc[(df2["query_id"] == gene) &
                              (df2["genome"] == genome), ]
            if tempdf1.empty or tempdf2.empty:
                logger.info("skipping %s in %s", gene, genome)
            else:
                subset1 = tempdf1.loc[
                    (tempdf1["identity_perc"] > min_percent)
                ]
                subset2 = tempdf2.loc[
                    (tempdf2["identity_perc"] > min_percent)
                ]
                logger.debug("grouped df shape: ")
                logger.debug(tempdf1.shape)
                logger.debug("grouped df2 shape: " )
                logger.debug(tempdf2.shape)
                if subset1.empty or subset2.empty:
                    logger.info("No reciprocol hits for %s in %s", gene, genome)
                    logger.debug(tempdf1)
                    logger.debug(tempdf2)
                    nonrecip_hits.append([gene, genome])
                else:
                    if subset1.iloc[0]["query_id"] == subset2.iloc[0]["subject_id"]:
                        recip_hits.append([gene, genome])
                        logger.debug("Reciprocol hits for %s in %s!", gene, genome)
                        if subset1.iloc[0]["alignment_length"] >= \
                           (min_lens[subset1.iloc[0]["query_id"]] - 0):
                            filtered = filtered.append(subset1)
                            logger.info("%s in %s passed min len test!", gene, genome)
                        else:
                            pass

                    else:
                        nonrecip_hits.append([gene, genome])
                        logger.debug("No reciprocol hits for %s in %s",
                                     gene, genome)

            # logger.debug(subset.shape)
    logger.debug("Non-reciprocal genes:")
    logger.debug(nonrecip_hits)
    logger.debug("Reciprocal genes:")
    logger.debug(recip_hits)
    logger.debug("filtered shape:")
    logger.debug(filtered.shape)
    return(filtered)


def checkBlastForMisjoin(df, fasta, ref_lens, BUF, flanking, logger=None):
    """ results from pd.read_csv with default BLAST output 6 columns
    returns a df
    """
    logger.debug("length of references:")
    logger.debug(ref_lens)
    df['name'] = df.query_id.str.replace("_upstream", "").str.replace("_downstream", "")
    # df['name2'] = df.name.str.replace("_downstream", "")
    df['query_name'] = df['name'].str.split('flanking').str.get(0)
    df['query_len'] = df.query_id.str.split('length_').str.get(1).str.split("_").str.get(0)
    where = []
    for i, row in df.iterrows():
        where.append("down" if "downstream" in row['query_id'] else "up")
    df['where'] = where
    assert logger is not None, "must use a logger"
    # print(ref_lens)
    print("\n")
    queries = df.query_name.unique()
    # subjects = df.subject_id.unique()
    sdf = df.loc[(df["alignment_length"] > (flanking * 0.9) - BUF)]
    naughty_nice_list = []
    for query in queries:
        logger.debug("checking hits for %s", query)
        tempdf = sdf.loc[(df["query_name"] == query)]
        for i, row in tempdf.iterrows():
            # print("outer row")
            subject_start = None
            # if both start the same (around 1), we have the first hit
            if row["s_start"] - 1 < BUF and abs(row["q_start"] - 1) < BUF:
                subject_start = row["subject_id"]
                ref_len = ref_lens[row["subject_id"]]
                logger.debug("checking %s and %s, len %d",
                             query, subject_start, ref_len)
                # print(tempdf)
                foundpair = False
                for i, innerrow in tempdf.iterrows():
                    subject_len = ref_lens[innerrow["subject_id"]]
                    subject_end = innerrow["subject_id"]
                    # if hit extends to end of reference
                    logger.debug("subject len: %s", subject_len)
                    logger.debug(innerrow)
                    logger.debug(abs(innerrow["s_end"] - subject_len))
                    if (abs(innerrow["s_end"] - subject_len)) < BUF:
                        # if same contig
                        if subject_start == subject_end:
                            naughty_nice_list.append(
                                [fasta, "good", query, subject_start, subject_end])
                            foundpair = True
                        else:
                            naughty_nice_list.append(
                                [fasta, "bad", query, subject_start, subject_end]
                            )
                            foundpair = True
                if not foundpair:
                    naughty_nice_list.append(
                        [fasta, "?", query, subject_start, "?"])
    print("Results for %s:" % fasta)
    for line in naughty_nice_list:
        print("\t".join(line))
    print("\n")
    return(naughty_nice_list)


def write_results(df, fasta_name, outfile, logger=None):
    #% parse output
    assert logger is not None, "must use a logger"
    logger.debug("writing out the results")
    with open(outfile, "a") as outf:
        outf.write("# {0} \n".format(fasta_name))
        df.to_csv(outf)


def parseDirContents(dirname, ref_ext, assembly_ext):
    """retursn a tuple (ref, [assembly1, assembly2, etc])
    """
    return (glob.glob(dirname + "*" + ref_ext)[0],
            glob.glob(dirname + "*" + assembly_ext))


def getScanCmd(ref, outroot, other_args):
    """ returns (cmd, path/to/dir/)
    """
    if other_args != "":
        other_args = " " + other_args  # pad with space for easier testing
    if ref.endswith(".gb"):
        return (None, ref)

    resulting_gb = os.path.join(outroot, "scan", "scannedScaffolds.gb")
    return ("ribo scan {0} --min_length 5000 -o {1}{2}".format(
        ref,
        os.path.join(outroot, "scan"),
        other_args
    ), resulting_gb)


def getSelectCmd(gb, outroot, other_args):
    resulting_clusters = os.path.join(outroot, "select",
                                      "riboSelect_grouped_loci.txt")
    if other_args != "":
        other_args = " " + other_args  # pad with space for easier testing
    return ("ribo select {0} -o {1}{2}".format(
        gb,
        os.path.join(outroot, "select"),
        other_args
    ), resulting_clusters)


def getSnagCmd(scangb, cluster, flank, outroot, other_args=""):
    if other_args != "":
        other_args = " " + other_args  # pad with space for easier testing
    return ("ribo snag {0} {1} -l {2} --just_extract -o {3}{4}".format(
        scangb,
        cluster,
        flank,
        os.path.join(outroot, "snag"),
        other_args
    ), os.path.join(outroot, "snag"))


def check_scan_select_snag_retruncodes(subreturns, logger):
    if subreturns[0].returncode != 0:
        logger.error("error with riboScan! Check the riboScan log files")
        sys.exit(1)
    if subreturns[1].returncode != 0:
        logger.error("error with riboSelect! Check the riboSelect log files")
        sys.exit(1)
    if subreturns[2].returncode != 0:
        logger.info("error with riboSnag! This often happens if " +
                    "the assembly doesnt reconstruct any rDNAs.")
        # note the lack of sys exit


def main(args, logger=None):
    if args.output is None:
        args.output = os.path.dirname(
            os.path.join(args.indir, "")
        ) + "_riboScored"
    output_root = os.path.abspath(os.path.expanduser(args.output))
    if not os.path.isdir(output_root):
        sys.stderr.write("creating output directory %s\n" % output_root)
        os.makedirs(output_root)
    else:
        sys.stderr.write("Output Directory already exists!\n")
        sys.exit(1)
    log_path = os.path.join(output_root, "riboScore.log")
    if logger is None:
        logger = set_up_logging(verbosity=args.verbosity,
                                outfile=log_path,
                                name=__name__)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    if not os.path.isdir(os.path.join(args.indir, "")) or len(
            os.listdir(os.path.join(args.indir, ""))) == 0:
        logger.error("input directory doesnt exist or is empty! Exiting...")
        sys.exit(1)
    gb, fastas = parseDirContents(dirname=os.path.join(args.indir, ""),
                                  ref_ext=args.ref_ext,
                                  assembly_ext=args.assembly_ext)

    # snags from reference
    bs_dir1 = os.path.join(output_root, "bridgeSeeds_ref")
    scancmd1, scangb1 = getScanCmd(ref=gb, outroot=bs_dir1, other_args="")
    selectcmd1, cluster1 = getSelectCmd(gb=scangb1, outroot=bs_dir1,
                                        other_args="-s 16S:23S")
    snagcmd1, snagdir1 = getSnagCmd(scangb=scangb1, cluster=cluster1,
                                    flank=args.flanking,
                                    outroot=bs_dir1,
                                    other_args="")
    logger.info(
        "Running riboScan, riboSelect, and riboSnag on reference: %s", gb)
    report_list = []
    for i in [scancmd1, selectcmd1, snagcmd1]:
        if i is None:
            continue
        logger.debug(i)
        subprocess.run(
            [i],
            shell=sys.platform != "win32",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True)
    for index, fasta in enumerate(fastas):
        logger.debug("processing %s", fasta)
        this_root = os.path.join(
            args.output, os.path.splitext(os.path.basename(fasta))[0])
        bs_dir2 = os.path.join(this_root, "bridgeSeeds_contigs")
        os.makedirs(bs_dir2)
        # snags from assembly
        scancmd2, scangb2 = getScanCmd(ref=fasta, outroot=bs_dir2,
                                       other_args='')
        selectcmd2, cluster2 = getSelectCmd(gb=scangb2, outroot=bs_dir2,
                                            other_args="-s 16S:23S")
        snagcmd2, snagdir2 = getSnagCmd(scangb=scangb2, cluster=cluster2,
                                        flank=args.flanking,
                                        outroot=bs_dir2)
        logger.info(
            "Running riboScan, riboSelect, and riboSnag on " +
            "%s, assembly %d of %d",
            fasta, index + 1, len(fastas))
        returncodes = []
        for i in [scancmd2, selectcmd2, snagcmd2]:
            if i is None:
                continue
            logger.debug(i)
            returncodes.append(subprocess.run(
                [i],
                shell=sys.platform != "win32",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False))  # we check later due to likely de novo failure
        check_scan_select_snag_retruncodes(
            subreturns=returncodes, logger=logger)

        ref_snags = sorted(glob.glob(
            snagdir1 + "/*_riboSnag.fasta"))

        if args.blast_full:
            full_blast_results = os.path.join(this_root, "BLAST")
            os.makedirs(full_blast_results)
            combined_full_snags = combine_contigs(
                contigs_dir=snagdir2,
                pattern="*riboSnag",
                contigs_name="combinedSnags",
                logger=logger)
            commands, paths_to_outputs, paths_to_recip_outputs = \
                make_nuc_nuc_recip_blast_cmds(
                    query_list=ref_snags,
                    subject_file=combined_full_snags,
                    output=full_blast_results,
                    logger=logger)
        else:
            commands = []

        contig_snags = sorted(glob.glob(
            os.path.join(snagdir2, "") +
            "*_riboSnag.fasta"))
        contig_snags_flanking = sorted(glob.glob(
            os.path.join(snagdir2, "flanking_regions_output", "") +
            "*_riboSnag_flanking_regions.fasta"))
        logger.debug(contig_snags)
        logger.debug(contig_snags_flanking)
        # combine the assembly contigs
        if len(contig_snags) == 0:
            report_list.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                os.path.abspath(os.path.expanduser(args.indir)),  # 0
                os.path.basename(fasta),  # 1
                len(ref_snags),  # 2
                0,  # 3
                0,  # 4
                0  # 5
            ))
            continue
        combined_flanking_snags = combine_contigs(
            contigs_dir=os.path.join(
                snagdir2, "flanking_regions_output", ""),
            pattern="*riboSnag_flanking_regions",
            contigs_name="combinedSnagFlanking",
            logger=logger)

        ref_snag_dict = {}
        contig_snag_dict = {}
        for snag in ref_snags:
            rec = SeqIO.read(snag, "fasta")
            ref_snag_dict[rec.id] = len(rec.seq)
        for snag in contig_snags:
            rec = SeqIO.read(snag, "fasta")
            contig_snag_dict[rec.id] = len(rec.seq)
        logger.debug(ref_snag_dict)
        logger.debug(contig_snag_dict)
        flanking_blast_results = os.path.join(this_root, "BLAST_flanking")
        os.makedirs(flanking_blast_results)
        f_commands, f_paths_to_outputs, f_paths_to_recip_outputs = \
            make_nuc_nuc_recip_blast_cmds(
                query_list=ref_snags,
                subject_file=combined_flanking_snags,
                output=flanking_blast_results,
                logger=logger)
        # check for existing blast results
        pool = multiprocessing.Pool()
        logger.debug("Running the following commands in parallel " +
                     "(this could take a while):")
        logger.debug("\n" + "\n".join([x for x in commands + f_commands]))
        logger.info("Running BLAST commands")
        results = [
            pool.apply_async(subprocess.run,
                             (cmd,),
                             {"shell": sys.platform != "win32",
                              "stdout": subprocess.PIPE,
                              "stderr": subprocess.PIPE,
                              "check": True})
            for cmd in commands + f_commands]
        pool.close()
        pool.join()
        reslist = []
        reslist.append([r.get() for r in results])
        logger.info("Parsing BLAST results")
        if args.blast_full:
            merged_tab = merge_outfiles(
                filelist=paths_to_outputs,
                outfile=os.path.join(this_root, "merged_results.tab"))
            recip_merged_tab = merge_outfiles(
                filelist=paths_to_recip_outputs,
                outfile=os.path.join(this_root, "recip_merged_results.tab"))
            resultsdf = BLAST_tab_to_df(merged_tab)
            recip_resultsdf = BLAST_tab_to_df(recip_merged_tab)
            filtered_hits = filter_recip_BLAST_df(
                df1=resultsdf,
                df2=recip_resultsdf,
                min_lens=ref_snag_dict,
                min_percent=args.min_percent,
                logger=logger)
            write_results(
                outfile=os.path.join(output_root,
                                     "riboScore_hits_fulllength.txt"),
                fasta_name=fasta,
                df=filtered_hits, logger=logger)

        f_merged_tab = merge_outfiles(
            filelist=f_paths_to_outputs,
            outfile=os.path.join(
                this_root, "merged_flanking_results.tab"))
        f_recip_merged_tab = merge_outfiles(
            filelist=f_paths_to_recip_outputs,
            outfile=os.path.join(
                this_root, "recip_merged_flanking_results.tab"))
        f_resultsdf = BLAST_tab_to_df(f_merged_tab)
        f_recip_resultsdf = BLAST_tab_to_df(f_recip_merged_tab)
        # 5 columns: [fasta, good/bad/?, query, startseq, end_seq]
        flanking_hits = checkBlastForMisjoin(
            fasta=fasta,
            df=f_recip_resultsdf,
            ref_lens=ref_snag_dict,
            flanking=args.flanking,
            BUF=50, logger=logger)
        with open(os.path.join(output_root, "riboScore_hits.txt"), "a") as f:
            for line in flanking_hits:
                f.write("\t".join(line) + "\n")
        good_hits = 0 + sum([1 for x in flanking_hits if x[1] == "good"])
        ambig_hits = 0 + sum([1 for x in flanking_hits if x[1] == "?"])
        bad_hits = 0 + sum([1 for x in flanking_hits if x[1] == "bad"])
        report_list.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
            os.path.abspath(os.path.expanduser(args.indir)),  # 0
            os.path.basename(fasta),  # 1
            len(ref_snags),  # 2
            good_hits,  # 3
            ambig_hits,  # 4
            bad_hits  # 5
        ))
        logger.debug(report_list)
    with open(os.path.join(output_root, "riboScore_report.txt"), "a") as r:
        for line in report_list:
            r.write(line)

#!/usr/bin/env python3
"""
"""
import os
import sys
import shutil
import datetime
import subprocess
import argparse
import multiprocessing
import logging
import time
import glob

from Bio import SeqIO
import pandas as pd
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from pyutilsnrw.utils3_5 import set_up_logging, combine_contigs

BUF=10

def get_args():
    parser = argparse.ArgumentParser(
        description="This does some simple reciprocol blasting to get region of interest")
    parser.add_argument("indir",
                        help="dir containing a genbank file and other file")
    parser.add_argument("-o", "--output", dest='output',
                        help="directory in which to place the output files",
                        default=os.path.join(os.getcwd(), "simpleOrtho"))
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
    args = parser.parse_args()
    return(args)


def make_nuc_nuc_recip_blast_cmds(
        query_list, date,
        output, subject_file=None, logger=None):
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


def merge_outfiles(filelist, outfile_name):
    """
    """
    # only grab .tab files, ie, the blast output
    filelist = [i for i in filelist if i.split(".")[-1:] == ['tab']]
    if len(filelist) == 1:
        # print("only one file found! no merging needed")
        return(filelist)
    else:
        # print("merging all the blast results to %s" % outfile_name)
        nfiles = len(filelist)
        fout = open(outfile_name, "a")
        # first file:
        for line in open(filelist[0]):
            fout.write(line)
        #  now the rest:
        for num in range(1, nfiles):
            f = open(filelist[num])
            for line in f:
                fout.write(line)
            f.close()  # not really needed
        fout.close()
    return(outfile_name)


def BLAST_tab_to_df(path):
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score"]
    raw_csv_results = pd.read_csv(
        open(path), comment="#", sep="\t", names=colnames)
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
                    # (tempdf1["bit_score"] == tempdf1["bit_score"].max())
                ]
                # (tempdf1["alignement_l"] == tempdf1["bit_score"].max())]
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
                    # logger.debug(tempdf1)
                    # logger.debug("tempdf2")
                    # logger.debug(tempdf2)
                    # logger.debug("subset1")
                    # logger.debug(subset1)
                    # logger.debug("subset2")
                    # logger.debug(subset2)
                    if subset1.iloc[0]["query_id"] == subset2.iloc[0]["subject_id"]:
                        recip_hits.append([gene, genome])
                        logger.debug("Reciprocol hits for %s in %s!", gene, genome)
                        # print(subset1.iloc[0]["query_id"])
                        # print(min_lens[subset1.iloc[0]["query_id"]])
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


# def checkBlastForMisjoin(df, contig_lens, min_percent, logger=None):
#     """ results from pd.read_csv with default BLAST output 6 columns
#     pineapple
#     returns a df
#     """
#     assert logger is not None, "must use a logger"
#     queries = df.query_id.unique()
#     for query in queries:
#         tempdf = df.loc[(df["query_id"] == query)]
#         contig_len = contig_lens[query]
#         subject_start = None
#         # print(contig_len)
#         for i, row in tempdf.iterrows():
#             # if both start the same
#             if abs(row["s_start"] - 1) < BUF and abs(row["q_start"] - 1) < BUF and row["alignment_length"] > args.flanking:
#                 subject_start = row["subject_id"]
#                 for j, innerrow in tempdf.iterrows():
#                     if innerrow["subject_id"] == subject_start:
#                         continue
#                     if abs(contig_len - innerrow["q_end"]) < BUF and innerrow["alignment_length"] > args.flanking:
#                         print ("YOP")
#                         print(subject_start)

def checkBlastForMisjoin(df, contig_lens, min_percent, logger=None):
    """ results from pd.read_csv with default BLAST output 6 columns
    pineapple
    returns a df
    """
    assert logger is not None, "must use a logger"
    print("contig_lens")
    print(contig_lens)
    queries = df.query_id.unique()
    for query in queries:
        tempdf = df.loc[(df["query_id"] == query)]
        contig_len = contig_lens[query]
        subject_start = None
        print(contig_len)
        for i, row in tempdf.iterrows():
            # if both start the same
            if abs(row["s_start"] - 1) < BUF and abs(row["q_start"] - 1) < BUF and row["alignment_length"] > args.flanking:
                subject_start = row["subject_id"]
                logger.info("checkign %s and %s, len %d",query,  subject_start, contig_len)
                for j, innerrow in tempdf.iterrows():
                    if abs(contig_len - innerrow["q_end"]) < BUF and innerrow["alignment_length"] > args.flanking:
                        if innerrow["subject_id"] == subject_start:
                            print(subject_start)
                            print ("POP")
                        else:
                            print("YOP")
                            print(subject_start)
        # if percent_identity <  threshold &
        # :
        #     continue



#     # df1['genome'] = df1.query_id.str.split('_').str.get(0)
#     # df2['genome'] = df2.subject_id.str.split('_').str.get(0)
#     df1['genome'] = df1.query_id
#     df2['genome'] = df2.subject_id
#     logger.debug(df1.shape)
#     logger.debug(df2.shape)
#     # recip structure
#     filtered = pd.DataFrame(columns=df1.columns)
#     unq_subject = df1.subject_id.unique()
#     unq_query = df1.genome.unique()
#     recip_hits = []
#     nonrecip_hits = []
#     for gene in unq_subject:
#         # for genome in unq_query:
#             logger.debug("Checking %s in %s for reciprocity" % (gene, genome))
#             tempdf1 = df1.loc[(df1["subject_id"] == gene) &
#                               (df1["genome"] == genome), ]
#             tempdf2 = df2.loc[(df2["query_id"] == gene) &
#                               (df2["genome"] == genome), ]
#             if tempdf1.empty or tempdf2.empty:
#                 logger.info("skipping %s in %s", gene, genome)
#             else:
#                 subset1 = tempdf1.loc[
#                     (tempdf1["identity_perc"] > min_percent)
#                     # (tempdf1["bit_score"] == tempdf1["bit_score"].max())
#                 ]
#                 # (tempdf1["alignement_l"] == tempdf1["bit_score"].max())]
#                 subset2 = tempdf2.loc[
#                     (tempdf2["identity_perc"] > min_percent)
#                 ]
#                 logger.debug("grouped df shape: ")
#                 logger.debug(tempdf1.shape)
#                 logger.debug("grouped df2 shape: " )
#                 logger.debug(tempdf2.shape)
#                 if subset1.empty or subset2.empty:
#                     logger.info("No reciprocol hits for %s in %s", gene, genome)
#                     logger.debug(tempdf1)
#                     logger.debug(tempdf2)
#                     nonrecip_hits.append([gene, genome])
#                 else:
#                     # logger.debug(tempdf1)
#                     # logger.debug("tempdf2")
#                     # logger.debug(tempdf2)
#                     logger.debug("subset1")
#                     logger.debug(subset1)
#                     logger.debug("subset2")
#                     logger.debug(subset2)
#                     if subset1.iloc[0]["query_id"] == subset2.iloc[0]["subject_id"]:
#                         recip_hits.append([gene, genome])
#                         logger.debug("Reciprocol hits for %s in %s!", gene, genome)
#                         # print(subset1.iloc[0]["query_id"])
#                         # print(min_lens[subset1.iloc[0]["query_id"]])
#                         if subset1.iloc[0]["alignment_length"] >= \
#                            (min_lens[subset1.iloc[0]["query_id"]] - 0):
#                             filtered = filtered.append(subset1)
#                             logger.info("%s in %s passed min len test!", gene, genome)
#                         else:
#                             pass

#                     else:
#                         nonrecip_hits.append([gene, genome])
#                         logger.debug("No reciprocol hits for %s in %s", gene, genome)

#             # logger.debug(subset.shape)
#     logger.debug("Non-reciprocal genes:")
#     logger.debug(nonrecip_hits)
#     logger.debug("Reciprocal genes:")
#     logger.debug(recip_hits)
#     logger.debug("filtered shape:")
#     logger.debug(filtered.shape)
#     return(filtered)


def write_results(df, fasta_name, outfile, logger=None):
    #% parse output
    assert logger is not None, "must use a logger"
    logger.debug("cleaning up the csv output")
    with open(outfile, "a") as outf:
        outf.write("# {0} \n".format(fasta_name))
        df.to_csv(outf)


def parseDirContents(dirname, ref_ext, assembly_ext):
    """retursn a tuple (ref, [assembly1, assembly2, etc])
    """
    return (glob.glob(dirname + "*" + ref_ext)[0],
            glob.glob(dirname + "*" + assembly_ext))


def getScanCmd(ref, outroot):
    """ returns (cmd, path/to/dir/)
    """
    # print(__file__)
    if ref.endswith(".gb"):
        return (None, ref)

    resulting_gb = os.path.join(outroot, "scan", "scannedScaffolds.gb")
    return ("{0} {1} {2} -o {3}".format(
        sys.executable,
        os.path.join(
            os.path.dirname(__file__),
            "riboScan.py"),
        ref,
        os.path.join(outroot, "scan")
    ), resulting_gb)


def getSelectCmd(gb, outroot):
    resulting_clusters = os.path.join(outroot, "select",
                                      "riboSelect_grouped_loci.txt")
    return ("{0} {1} {2} -s 16S:23S -o {3}".format(
        sys.executable,
        os.path.join(
            os.path.dirname(__file__),
            "riboSelect.py"),
        gb,
        os.path.join(outroot, "select")
    ), resulting_clusters)


def getSnagCmd(scangb, cluster, flank, outroot):
    return ("{0} {1} {2} {3} -l {4} --just_extract -o {5}".format(
        sys.executable,
        os.path.join(
            os.path.dirname(__file__),
            "riboSnag.py"),
        scangb, cluster,
        flank,
        os.path.join(outroot, "snag")
    ), os.path.join(outroot, "snag"))


def main(args):
    EXISTING_DIR = False
    output_root = os.path.abspath(os.path.expanduser(args.output))
    if not os.path.isdir(output_root):
        sys.stderr.write("creating output directory %s\n" % output_root)
        os.makedirs(output_root)
    else:
        sys.stderr.write("Output Directory already exists!\n")
        EXISTING_DIR = True
        # sys.exit(1)
    log_path = os.path.join(output_root,
                            "riboScore_log.txt")
    logger = set_up_logging(verbosity=args.verbosity,
                            outfile=log_path,
                            name=__name__)

    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    gb, fastas = parseDirContents(dirname=os.path.join(args.indir, ""),
                                  ref_ext=args.ref_ext,
                                  assembly_ext=args.assembly_ext)
    for fasta in fastas:
    # fasta = fastas[0]
        this_root = os.path.join(
            args.output, os.path.splitext(os.path.basename(fasta))[0])
        bs_dir1 = os.path.join(this_root, "bridgeSeeds_ref")
        bs_dir2 = os.path.join(this_root, "bridgeSeeds_contigs")
        os.makedirs(bs_dir1)
        os.makedirs(bs_dir2)
        # snags from reference
        scancmd1, scangb1 = getScanCmd(ref=gb, outroot=bs_dir1)
        selectcmd1, cluster1 = getSelectCmd(gb=scangb1, outroot=bs_dir1)
        snagcmd1, snagdir1 = getSnagCmd(scangb=scangb1, cluster=cluster1,
                                        flank=args.flanking,
                                        outroot=bs_dir1)
        # snags from assembly
        scancmd2, scangb2 = getScanCmd(ref=fasta, outroot=bs_dir2)
        selectcmd2, cluster2 = getSelectCmd(gb=scangb2, outroot=bs_dir2)
        snagcmd2, snagdir2 = getSnagCmd(scangb=scangb2, cluster=cluster2,
                                        flank=args.flanking,
                                        outroot=bs_dir2)
        for i in [scancmd1, selectcmd1, snagcmd1,
                  scancmd2, selectcmd2, snagcmd2]:
            if i is None:
                continue
            logger.debug(i)
            subprocess.run([i],
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)

        ref_snags = sorted(glob.glob(snagdir1 + "/*_riboSnag.fasta"))

        contig_snags = sorted(glob.glob(snagdir2 + "/*_riboSnag.fasta"))
        # print(contig_snags)
        # contig_len_str = [
        #     os.path.basename(x).replace(".fasta", "").split(
        #         "_RC")[0].split("_")[-1].split("..")
        #     for x in contig_snags]
        # print(contig_len_str)
        # contig_lens = [int(x[1]) - int(x[0]) for x in contig_len_str]
        # print(contig_lens)
        combined_snags = combine_contigs(  # combine all the assembly contigs
            contigs_dir=snagdir2,
            pattern="*riboSnag",
            contigs_name="combinedSnags",
            logger=logger)

        ref_snag_dict = {}
        contig_snag_dict = {}
        for snag in ref_snags:
            rec = SeqIO.read(snag, "fasta")
            ref_snag_dict[rec.id] = len(rec.seq)
        for snag in contig_snags:
            rec = SeqIO.read(snag, "fasta")
            contig_snag_dict[rec.id] = len(rec.seq)
        # print(ref_snags)
        combined_df = pd.DataFrame()
        # for i, snag in enumerate(snags):
        blast_results = os.path.join(this_root, "BLAST")

        os.makedirs(blast_results)
        commands, paths_to_outputs, paths_to_recip_outputs = \
            make_nuc_nuc_recip_blast_cmds(
                query_list=ref_snags,
                subject_file=combined_snags,
                output=blast_results, date=date,
                logger=logger)
        # check for existing blast results
        if not all([os.path.isfile(x) for x in paths_to_outputs]):
            if EXISTING_DIR:
                logger.error("existing output dir found, but not all the " +
                             "needed blast results were found. Cannot use " +
                             "this directory")
                sys.exit(1)
            pool = multiprocessing.Pool()
            logger.debug("Running the following commands in parallel " +
                         "(this could take a while):")
            logger.debug("\n" + "\n".join([x for x in commands]))
            results = [
                pool.apply_async(subprocess.run,
                                 (cmd,),
                                 {"shell": sys.platform != "win32",
                                  "stdout": subprocess.PIPE,
                                  "stderr": subprocess.PIPE,
                                  "check": True})
                for cmd in commands]
            pool.close()
            pool.join()
            reslist = []
            reslist.append([r.get() for r in results])
        else:
            pass

        merged_tab = os.path.join(this_root,
                                  "merged_results.tab")
        recip_merged_tab = os.path.join(this_root,
                                        "recip_merged_results.tab")
        merge_outfiles(filelist=paths_to_outputs,
                       outfile_name=merged_tab)
        merge_outfiles(filelist=paths_to_recip_outputs,
                       outfile_name=recip_merged_tab)
        resultsdf = BLAST_tab_to_df(merged_tab)
        recip_resultsdf = BLAST_tab_to_df(recip_merged_tab)
        filtered_hits = filter_recip_BLAST_df(
            df1=resultsdf,
            df2=recip_resultsdf,
            min_lens=ref_snag_dict,
            min_percent=args.min_percent,
            logger=logger)
        print(ref_snag_dict)
        checkBlastForMisjoin(df=recip_resultsdf,
                             contig_lens=contig_snag_dict,
                             min_percent=44, logger=logger)
        # combined_df = combined_df.append(filtered_hits)
        # print(combined_df)
        write_results(
            outfile=os.path.join(output_root, "riboScore_hits.txt"),
            fasta_name=fasta,
            df=filtered_hits, logger=logger)

if __name__ == '__main__':
    assert ((sys.version_info[0] == 3) and (sys.version_info[1] >= 5)), \
        "Must use python3.5 or higher!"
    args = get_args()
    main(args=args)

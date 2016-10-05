# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
The Goal of this is to have a unified place to put the useful
python 3.5 functions or templates

how I got the fastq file
# seqtk sample -s 27 ~/GitHub/FA/pseudochromosome/data/20150803_Abram1/ \
    reads/3123-1_1_trimmed.fastq .0005

bam file was from a riboseed mapping; md5: 939fbf2c282091aec0dfa278b05e94ec

mapped bam was made from bam file with the following command
 samtools view -Bh -F 4 /home/nicholas/GitHub/FB/Ecoli_comparative_genomics/
    scripts/riboSeed_pipeline/batch_coli_unpaired/map/
    mapping_20160906_region_7_riboSnag/
    test_smalt4_20160906_region_7_riboSnagS.bam >
     ~/GitHub/pyutilsnrw/tests/test_mapped.sam
md5: 27944249bf064ba54576be83053e82b0

"""
__version__ = "0.0.3"
import time
import sys
import shutil
import logging
import subprocess
import os
import unittest
import hashlib
import glob
import argparse
sys.dont_write_bytecode = True

from pyutilsnrw.utils3_5 import make_output_prefix, check_installed_tools,\
    copy_file, get_ave_read_len_from_fastq, get_number_mapped,\
    extract_mapped_and_mappedmates, keep_only_first_contig, md5,\
    combine_contigs, clean_temp_dir

from riboSeed.riboseed import  check_smalt_full_install,\
    map_to_ref_smalt, convert_bams_to_fastq



@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among otherthings wont run if you try this" +
                 " with less than python 3.5")
class utils3_5TestCase(unittest.TestCase):
    def setUp(self):
        self.ref_gb
        self.ref_fasta
        self.ref_Ffastq
        self.ref_rfastq


        pass

    def test_map_to_ref_smalt(ref, ref_genome, fastq_read1, fastq_read2,
                              distance_results,
                              map_results_prefix, cores, samtools_exe,
                              smalt_exe, fastq_readS="",
                              read_len=100, step=3, k=5,
                              scoring="match=1,subst=-4,gapopen=-4,gapext=-3"):

def convert_bams_to_fastq(map_results_prefix, fastq_results_prefix,
                          keep_unmapped):

get_filtered_locus_tag_dict(genome_seq_record, feature="rRNA",
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




    def tearDown(self):
        pass

if __name__ == '__main__':
    args = get_args()
    curdir = os.getcwd()
    # samtools_exe = args.samtools_exe
    samtools_exe = "samtools"
    check_smalt_full_install(smalt_exe, logger=None)
    unittest.main()

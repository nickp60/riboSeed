#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.


# This config file was generated Wed Sep 13 17:01:28 2017


#------------------------#
##  Required programs   ##
#------------------------#
# executable for barrnap
BARRNAP_EXE : /home/nicholas/barrnap/bin/barrnap

# executable for seqret
SEQRET_EXE : /home/nicholas/bin/seqret

# executable for spades.py
SPADES_EXE : /home/nicholas/bin/SPAdes-3.9.1-Linux/bin/spades.py

# executable for bwa
BWA_EXE : /usr/bin/bwa

# executable for samtools
SAMTOOLS_EXE : /usr/bin/samtools

# executable for blastn
BLAST_EXE : /usr/bin/blastn

#------------------------#
##  Optional programs   ##
#------------------------#
# executable for quast.py
QUAST_EXE : /home/nicholas/bin/quast-4.5/quast.py

# executable for smalt
SMALT_EXE : /usr/local/bin/smalt

#------------------------#
##   Mauve (Optional)   ##
#------------------------#
# Muave is used by riboSketch for orienting contigs and 
# determining synteny.  riboSketch uses the Mauve.jar java program,
# but usually only mauveAligner (in the platform specific subdir)
# is actually added to the path after installation.
# executable for mauveAligner
MAUVE_ALIGNER_EXE : /home/nicholas/mauve_snapshot_2015-02-13/linux-x64/mauveAligner

# executable for Mauve.jar
MAUVE_JAR : /home/nicholas/mauve_snapshot_2015-02-13/Mauve.jar

#------------------------#
##  riboScan Parameters ##
#------------------------#
# --id_thresh: partial rRNA hits below this threshold will be ignored
SCAN_ID_THRESH : 0.5

# --name: stem to give the files generated by riboScan  
SCAN_CONTIG_NAME : null

# --min_length: skip annotating if contig shorter than this
SCAN_MIN_LENGTH : 0

# -v: verbosity for riboScan
SCAN_VERBOSITY : 2

#------------------------#
## riboSelect Parameters #
#------------------------#
# --feature: which annotations to pay attention to; 
# barrnap uses 'rRNA', but others may use 'RRNA', etc
SELECT_FEATURE : 'rRNA'

# -v: verbosity for riboSelect
SELECT_VERBOSITY : 2

#------------------------#
##  riboSeed Parameters ##
#------------------------#
# 
# --method_for_map: which mapper to use, (default bwa)
SEED_MAP_METHOD : bwa

# 
# --score_min: If using smalt, this sets the '-m' param; 
# default with smalt is inferred from 
# read length. If using BWA, reads mapping with AS
# score lower than this will be rejected
# ; default with BWA is half of read length
SEED_SCORE_MIN : null

# 
# if initial SPAdes assembly largest contig 
# is not at least as long as --min_assembly_len, 
# reject. Set this to the length of the seed 
# sequence; if it is not achieved, seeding across 
# regions will likely fail; default: %(default)s
SEED_MIN_ASSEMBLY_LENGTH : 6000

# 
# if assembled contig is smaller than --min_assembly_len, contig
#will still be included in assembly; default: inferred
SEED_INCLUDE_SHORTS : false

# 
# if --subtract reads already used in previous
# round of subassembly will not be included in 
# subsequent rounds.  This can lead to problems 
# with multiple mapping and inflated coverage.
SEED_SUBTRACT : false

# 
# if --skip_control, no de novo 
# assembly will be done; default: %(default)s
SEED_SKIP_CONTROL : false

# 
# ignore: reference will not be used in 
# subassembly. trusted: SPAdes will use the seed
#  sequences as a --trusted-contig; untrusted: 
# SPAdes will treat as --untrusted-contig. 
# infer: if mapping 
#percentage over 80%: 'trusted', else 'untrusted'
# See SPAdes docs for details.  default: infer
SEED_REF_AS_CONTIG : infer

# 
# If set, iterations will continue until 
# contigs reach this length, or max iterations (
#set by --iterations) have been completed. Set as 
#fraction of original seed length by giving a 
#decimal between 0 and 5, or set as an absolute 
#number of base pairs by giving an integer greater
# than 50. Not used by default
SEED_TARGET_LEN : null

# 
# Submit custom parameters to mapper. 
#And by mapper, I mean bwa, cause we dont support 
#this option for SMALT, sorry. 
#This requires knowledge of your chosen mapper's 
#optional arguments. Proceed with caution!
SEED_MAPPER_ARGS : -L 0,0 -U 0 -a

# 
# If mapping with SMALT, 
#submit custom smalt scoring via smalt -S 
#scorespec option
SEED_SMALT_SCORING : match=1,subst=-4,gapopen=-4,gapext=-3

# 
# -v: verbosity for riboSeed
SEED_VERBOSITY : 2

#------------------------#
## riboSketch Parameters #
#------------------------#
# 
# -f extension of assemblies, usually fasta
SKETCH_ASSEMBLY_EXT : 'fasta'

# 
# -g extension of reference, usually gb
SKETCH_REF_EXT : 'gb'

# 
# -v: verbosity for riboSketch
SKETCH_VERBOSITY : 2

#------------------------#
## riboScore Parameters ##
#------------------------#
# NOTE: other than verbosity, there are no parameters that are 
# unique to riboScore, as all the arguments have den defined in 
# previous programs above. See run_riboSeed.py for more info as to 
# which args are used.
# -v: verbosity for riboSketch
SCORE_VERBOSITY : 2

#------------------------#
##    Run Parameters    ##
#------------------------#
# NOTE: This section will only be filled in when run by run_riboSeed
# as this contains info specific to a particular run.  This feature 
# is aimed at improving reproducibility.  See the docs for details 
# about the args

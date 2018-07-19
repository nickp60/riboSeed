#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Copyright 2017, National University of Ireland and The James Hutton Insitute
# Author: Nicholas Waters
#
# This code is part of the riboSeed package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

import pkg_resources
import sys
import os
import shutil
import subprocess

helpstring = """
Welcome to the runme! Here we test the integration of several parts of the
riboSeed pipeline.  First, `ribo run` is performed on the included test
dataset.  Then, essentially the same thing is done, but calling the
individual subcommands (`ribo scan`, `ribo select`, etc)

If all goes well, no errors should occur, and you should essentially have
two "identical" riboSeed assemblies (although due to random assignments
of mapping duplicates, the nature of error correction, etc, I can't
guarantee that you will get the exact same result

Have fun!
"""


resource_package = pkg_resources.Requirement.parse("riboSeed")  # Could be any module/package name
print(resource_package)
resource_path_fasta = '/'.join(('riboSeed',
                                'integration_data', 'concatenated_seq.fasta'))
resource_path_reffasta = '/'.join(('riboSeed',
                                   'integration_data', 'NC_000913.3.fasta'))
resource_path_1 = '/'.join(('riboSeed',
                            'integration_data', 'test_reads1.fq'))
resource_path_2 = '/'.join(('riboSeed',
                            'integration_data', 'test_reads2.fq'))

print(resource_path_fasta)
fasta = pkg_resources.resource_filename(resource_package, resource_path_fasta)
reffasta = pkg_resources.resource_filename(resource_package,
                                           resource_path_reffasta)
fastq1 = pkg_resources.resource_filename(resource_package, resource_path_1)
fastq2 = pkg_resources.resource_filename(resource_package, resource_path_2)
# fasta_path = pkg_resources.resource_string("/", resource_path)
print(fasta)
print(fastq1)
print(fastq2)

for i in ["blastn", "spades.py", "bwa", "mafft",
          "samtools", "seqret", "barrnap"]:
    assert shutil.which(i) is not None, "{0} executable not found in PATH!".format(i)

root_dir = os.path.join(os.getcwd(), "integration_test_results")
os.makedirs(root_dir, exist_ok=True)
cmds = [
    # riboScan
    "ribo scan {0}{1} -o {2} -k bac".format(
        os.path.join(os.path.dirname(fasta), ""),
        os.path.basename(fasta),
        os.path.join(root_dir, "scan", "")),
    # riboSelect
    "ribo select {0}scannedScaffolds.gb -o {1}".format(
        os.path.join(root_dir, "scan", ""),
        os.path.join(root_dir, "select","")),
    # riboSnag
    "ribo snag {0}scannedScaffolds.gb {1}riboSelect_grouped_loci.txt -o {2}".format(
        os.path.join(root_dir, "scan", ""),
        os.path.join(root_dir, "select", ""),
        os.path.join(root_dir, "snag", "")),
    # riboSeed
    "ribo seed -r {0}scannedScaffolds.gb {1}riboSelect_grouped_loci.txt -o {2} -F {3} -R {4} --serialize -v 1 ".format(
        os.path.join(root_dir, "scan", ""),
        os.path.join(root_dir, "select", ""),
        os.path.join(root_dir, "seed", ""),
        fastq1,
        fastq2)
]

ribo_run_cmd = \
    "ribo run -r {0} -o {1} -F {2} -R {3} --serialize -v 1 --stages stack score spec".format(
        fasta,
        os.path.join(root_dir, "run"),
        fastq1,
        fastq2)
print("running " + ribo_run_cmd)

subprocess.run([ribo_run_cmd],
               shell=sys.platform != "win32",
               stdout=subprocess.PIPE,
               stderr=subprocess.PIPE,
               check=True)
print("finished running integration test with ribo run!")


print("running integration tests with individual modules")
for i in cmds:
    print(i)
    subprocess.run([i],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
print("finished running integration test on parts!")

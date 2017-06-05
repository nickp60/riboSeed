#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import pkg_resources
import sys
import os
import shutil
import subprocess
# sys.path.append(os.path.join('..', 'riboSeed'))
# from riboSeed import riboSnag as rs

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

for i in ["blastn", "spades.py", "quast.py", "bwa", "mafft",
          "samtools", "seqret", "barrnap"]:
    assert shutil.which(i) is not None, "{0} executable not found in PATH!".format(i)

root_dir = os.path.join(os.getcwd(), "integration_test_results")
os.makedirs(root_dir, exist_ok=True)
cmds = [
    # riboScan
    "{0} {1} {2}{3} -o {4} -k bac".format(
        sys.executable,
        shutil.which("riboScan.py"),
        os.path.join(os.path.dirname(fasta), ""),
        os.path.basename(fasta), os.path.join(root_dir, "scan", "")),
    # riboSelect
    "{0} {1} {2}scannedScaffolds.gb -o {3}".format(
        sys.executable,
        shutil.which("riboSelect.py"),
        os.path.join(root_dir, "scan", ""),
        os.path.join(root_dir, "select","")),
    # riboSnag
    "{0} {1} {2}scannedScaffolds.gb {3}riboSelect_grouped_loci.txt -o {4} ".format(
        sys.executable,
        shutil.which("riboSnag.py"),
        os.path.join(root_dir, "scan", ""),
        os.path.join(root_dir, "select", ""),
        os.path.join(root_dir, "snag", "")),
    # riboSeed
    "{0} {1} -r {2} {3}riboSelect_grouped_loci.txt -o {4} -F {5} -R {6} --serialize ".format(
        sys.executable,
        shutil.which("riboSeed.py"),
        reffasta,
        os.path.join(root_dir, "select", ""),
        os.path.join(root_dir, "seed", ""),
        fastq1,
        fastq2)
]
for i in cmds:
    print(i)
    subprocess.run([i],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
print("finished running integration test!")

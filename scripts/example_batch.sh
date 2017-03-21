#!/bin/bash
# version 0.1.0
#$ -cwd
#$ -j yes
#$ -V

echo "Hostname: $HOSTNAME"
echo "Task started `date`"

set -e # exit on error!
# Runs a riboseed assembly with default settings

# output scanScaffolds_combined.gb in current directory
echo 'USAGE: /path/to/genome.gb path/to/genome.fasta path/to/read1 path/to/read2 /path/to/outdir/ n_iterations n_flanking n_cores'

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] ;
then
    echo "All mandatory arguments: genome.gb, genome.fasta, read1, read2, output_dir, iterations, flanking_width, virtenve_exe"
    exit 1
fi
echo "example_script.sh PATH: \n $PATH"
if [ -d "$4" ]; then
    echo "output dir exists!"
    echo "play things safe, make a new directory for results"
    exit 1
fi
FASTA="$1"
READ1="$2"
READ2="$3"
OUTDIR="$4"
ITERATIONS="$5"
FLANK="$6"
VENVEXE="$7"
## make  dirs
mkdir ${OUTDIR}
## enter virtual environment
source ${VENVEXE}
## riboScan
riboScan.py ${FASTA} -o ${OUTDIR}/scan/ -v 1
## RiboSelect
riboSelect.py ${OUTDIR}/scan/scannedScaffolds.gb -o ${OUTDIR}select/ -v 1
## RiboSeed
riboSeed.py ${OUTDIR}select/riboSelect_grouped_loci.txt -F ${READ1} -R ${READ2} -l ${FLANK} -r ${OUTDIR}/scan/scannedScaffolds.gb -i ${ITERATIONS} -o ${OUTDIR}seed/ -v 1 -c 4 -t 2

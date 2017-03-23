#!/bin/bash
# version 0.2.0
#  If submitting jobs to a cluster with SGE, set the variables
#  below and the environment variables starting with #$.  Otherwise, you
#  can submit arguments to commandline
#$ -cwd
#$ -j yes
#$ -V
#$ -N name
#$ -pe mpi 4

echo "Hostname: $HOSTNAME"
echo "Task started `date`"
echo "if running without qsub:"
echo "USAGE\n sge_batch.sh genome.gb reads_f.fastq reads_r.fastq /path/to/outdir/ n_iterations, n_flanking, venv_exe"
echo "sge_batch.sh PATH: \n $PATH"

set -e # exit on error!

# set these variables
FASTA="$1"
READ1="$2"
READ2="$3"
ITERATIONS="$4"
FLANK="$5"
VENVEXE="$6"
NAME="$7"
echo "$PATH"

OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_${FLANK}bp_${ITERATIONS}x_${NAME}/"

# Runs a riboseed assembly with default settings

# output
if [ -d "$OUTDIR" ]; then
    echo "output dir exists!"
    echo "play things safe, make a new directory for results"
    exit 1
fi

## make  dirs
mkdir ${OUTDIR}
## enter virtual environment
source ${VENVEXE}
## riboScan
riboScan.py ${FASTA} -o ${OUTDIR}scan/ -v 1 -n ${NAME}
## RiboSelect
riboSelect.py ${OUTDIR}/scan/scannedScaffolds.gb -o ${OUTDIR}select/ -v 1
## RiboSeed
riboSeed.py ${OUTDIR}select/riboSelect_grouped_loci.txt -F ${READ1} -R ${READ2} -l ${FLANK} -r ${OUTDIR}scan/scannedScaffolds.gb -i ${ITERATIONS} -o ${OUTDIR}seed/ -v 1 -c 4 -t 2 -n ${NAME}


echo "Task ended `date`"

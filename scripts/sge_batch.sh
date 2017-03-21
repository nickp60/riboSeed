#!/bin/bash
# version 0.1.0
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

FASTA="$1"
READ1="$2"
READ2="$3"
OUTDIR="$4"
ITERATIONS="$5"
FLANK="$6"
VENVEXE="$7"
echo "$PATH"
# include any commands and/or other
# bash scripting instructions here
./scripts/example_batch.sh $FASTA $READ1 $READ2 $OUTDIR $ITERATIONS $FLANK $VENVEXE

echo "Task ended `date`"

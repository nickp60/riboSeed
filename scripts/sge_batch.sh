#!/bin/bash
# version 0.0.31
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
echo "USAGE\n sge_batch.sh genome.gb genome.fasta reads_f.fastq reads_r.fastq /path/to/outdir/ n_iterations, n_flanking, venv_exe"
echo "sge_batch.sh PATH: \n $PATH"

GB="$1"
FASTA="$2"
READ1="$3"
READ2="$4"
OUTDIR="$5"
ITERATIONS="$6"
FLANK="$7"
VENVEXE="$8"
echo "$PATH"
# include any commands and/or other
# bash scripting instructions here
./scripts/example_batch.sh $GB $FASTA $READ1 $READ2 $OUTDIR $ITERATIONS $FLANK $VENVEXE

echo "Task ended `date`"

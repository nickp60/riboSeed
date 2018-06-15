#!/bin/bash
# version 0.0.4
#  If submitting jobs to a cluster with SGE, set the variables
#  below and the environment variables starting with #$.  Otherwise, you
#  can submit arguments to commandline
#$ -cwd
#$ -j yes
#$ -V
#$ -N simKleb
#$ -pe mpi 8

if [ -z "$1" ]
then
    echo "3 required arguments: forward reads, reverse reads, name. \n USAGE: mini_assembly.sh reads1 reads2 NAME"
    exit 1
else
    NAME="$3"
fi



OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_mini_assembly_${NAME}/"
mkdir ${OUTDIR}
mkdir ${OUTDIR}/tmp/


seqtk sample -s100 $1 10000 > ${OUTDIR}/tmp/reads1.fq
seqtk sample -s100 $2 10000 > ${OUTDIR}/tmp/reads2.fq

spades.py -1 ${OUTDIR}/tmp/reads1.fq -2 ${OUTDIR}/tmp/reads2.fq -o ${OUTDIR}/spades/

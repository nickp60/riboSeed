#!/bin/bash

if [ -z "$3" ]
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

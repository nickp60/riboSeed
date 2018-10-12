#!/bin/bash
# This runs  a mini assembly of a subset of reads to use with the average nucleotide identity analysis
if [ -z "$3" ]
then
    echo "3 required arguments: forward reads, reverse reads, outdir. \n USAGE: mini_assembly.sh reads1 reads2 /path/to/outdir"
    exit 1
else
    OUTDIR="$3"
fi


mkdir ${OUTDIR}/downsampled_reads/

# seqtk sample -s100 $1 10000 > ${OUTDIR}/tmp/reads1.fq
# seqtk sample -s100 $2 10000 > ${OUTDIR}/tmp/reads2.fq
{
    seqtk sample -s100 $1 100000 > ${OUTDIR}/downsampled_reads/reads1.fq
    seqtk sample -s100 $2 100000 > ${OUTDIR}/downsampled_reads/reads2.fq
} || {
    seqtk sample -s100 $1 .1    > ${OUTDIR}/downsampled_reads/reads1.fq
    seqtk sample -s100 $2 .1    > ${OUTDIR}/downsampled_reads/reads2.fq
}

# spades.py -1 ${OUTDIR}/downsampled_reads/reads1.fq -2 ${OUTDIR}/downsampled_reads/reads2.fq -o ${OUTDIR}/assembly/ -t 4 -m 64
mkdir $OUTDIR/assembly/
skesa --fastq ${OUTDIR}/downsampled_reads/reads1.fq ${OUTDIR}/downsampled_reads/reads2.fq --contigs_out ${OUTDIR}/assembly/contigs.fasta

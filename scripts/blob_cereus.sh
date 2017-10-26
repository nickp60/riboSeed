#!/bin/bash
# This script is what we used to determine the degree of contamination present
# in the GAGE-B HiSeq B. cereus dataset.
# 2017-10-23
# Nick Waters

REF=~/Downloads/hi/cereus/AE017194.1.fasta
READS1=~/Downloads/hi/cereus/trimmed/insert_180_1__cov250x.fastq
READS2=~/Downloads/hi/cereus/trimmed/insert_180_2__cov250x.fastq
OUTDIR="./blob_cereus_output/"
LOG=${OUTDIR}blob.log

mkdir $OUTDIR
echo "index and map reads to reference" 
bwa index $REF >> $LOG 2>&1
BWASAM=${OUTDIR}mapped_reads.sam
bwa mem $REF $READS1 $READS2 > $BWASAM
echo " extract reads mapping" 
READS1_MAPPED=${OUTDIR}mapped_1.fastq
READS2_MAPPED=${OUTDIR}mapped_2.fastq
READS1_UNMAPPED=${OUTDIR}unmapped_1.fastq
READS2_UNMAPPED=${OUTDIR}unmapped_2.fastq
samtools flagstat $BWASAM >> $LOG 2>&1
samtools fastq -F 12 $BWASAM -1 $READS1_MAPPED -2 $READS2_MAPPED >> $LOG 2>&1
samtools fastq -f 12 $BWASAM -1 $READS1_UNMAPPED -2 $READS2_UNMAPPED >> $LOG 2>&1
echo " assembling raw reads and subset reads with metaspades"
SPADES_ALL=${OUTDIR}metaspades_all/
SPADES_MAPPED=${OUTDIR}metaspades_mapped/
SPADES_UNMAPPED=${OUTDIR}metaspades_unmapped/
metaspades.py --pe1-1 $READS1 --pe1-2 $READS2 -o $SPADES_ALL
metaspades.py --pe1-1 $READS1_MAPPED --pe1-2 $READS2_MAPPED -o $SPADES_MAPPED
metaspades.py --pe1-1 $READS1_UNMAPPED --pe1-2 $READS2_UNMAPPED -o $SPADES_UNMAPPED
BDB=~/BLASTDB/nt
echo "Running blast for mapped"
HITS_MAPPED=${OUTDIR}mapped_blastn_nt.out
HITS_UNMAPPED=${OUTDIR}unmapped_blastn_nt.out
HITS_ALL=${OUTDIR}all_blastn_nt.out
blastn -db $BDB -query ${SPADES_MAPPED}contigs.fasta -out $HITS_MAPPED -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads 12
echo "Running blast for unmapped"
blastn -db $BDB -query ${SPADES_UNMAPPED}contigs.fasta -out $HITS_UNMAPPED -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads 12
echo "Running blast for all"
blastn -db $BDB -query ${SPADES_ALL}contigs.fasta -out $HITS_ALL -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -num_threads 12

echo "running blobtools cmds"
BLOB_ALL=${OUTDIR}blob_all/
BLOB_MAPPED=${OUTDIR}blob_mapped/
BLOB_UNMAPPED=${OUTDIR}blob_unmapped/
mkdir $BLOB_ALL
mkdir $BLOB_MAPPED
mkdir $BLOB_UNMAPPED
blobtools create -i ${SPADES_MAPPED}contigs.fasta  -y spades -o $BLOB_MAPPED -t $HITS_MAPPED  
blobtools create -i ${SPADES_UNMAPPED}contigs.fasta  -y spades -o $BLOB_UNMAPPED -t $HITS_UNMAPPED  
blobtools create -i ${SPADES_ALL}contigs.fasta  -y spades -o $BLOB_ALL -t $HITS_ALL  

echo "Generating plots at the Species level" 
blobtools blobplot -i ${BLOB_MAPPED}blobDB.json -o ${BLOB_MAPPED}mapped --rank species
blobtools blobplot -i ${BLOB_UNMAPPED}blobDB.json -o ${BLOB_UNMAPPED}unmapped --rank species
blobtools blobplot -i ${BLOB_ALL}blobDB.json -o ${BLOB_ALL}all --rank species

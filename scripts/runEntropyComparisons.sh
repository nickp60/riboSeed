#!/bin/bash

# if [ -d "assembly_statistics.txt" ]
# then
#     wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
# fi

while read refacc refother genus species strain gene geneacc x1
do
    echo "processing $genus $species  $strain $geneacc"
    # # get our reference protein sequence for the gene neighboring the first rDNA
    # curl -o ${geneacc}.fasta  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$geneacc&rettype=fasta&retmode=text"
    # # fetch 25 other genomes to compare ours to
    # Rscript ../manuscript_results/entropy/selectGenomes.R ./assembly_summary.txt ${genus}_${species}_${strain}_${geneacc}_output/ $genus $species ;
    echo "$refacc $refother"

done <./test.csv

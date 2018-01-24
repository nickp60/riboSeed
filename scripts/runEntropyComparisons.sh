#!/bin/bash

# if [ -d "assembly_statistics.txt" ]
# then
#     wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
# fi


while read refacc refother genus species strain gene geneacc
do
    echo "Downloading $genus $species $strain genomes..."
    thisdir="${genus}_${species}_${strain}"
    mkdir $thisdir
    cd $thisdir
    # get our reference protein sequence for the gene neighboring the first rDNA
    curl -o ${geneacc}.fasta  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$geneacc&rettype=fasta&retmode=text"
    # fetch 25 other genomes to compare ours to
    echo "Running R script"
    Rscript ../../scripts/getCompleteGenomeSubset.R ../assembly_summary.txt genome_selection $genus $species
    # unzip all the new genomes
    gunzip genome_selection/genomes/*
    cd ../
done < ./entropy_manifest.tab

# while read refacc refother genus species strain gene geneacc
# do
#     echo "Processing $genus $species $strain $geneacc ..."
#     # get the interesting regions
#     ~/GitHub/open_utils/orthoML/simpleOrtho.py ./${geneacc}.fasta ./genome_selection/genomes/ > regions.txt
#     # combine all genomes to create a single file
#     cat ./genome_selection/genomes/*.fna > combined_genomes.fasta
#     # use extractregions.py to pull out those regions, butwith 10kb on each side.  this all ends up in one file with each of the 25 loci
#     while read line
#     do
# 	echo $line | ~/GitHub/open_utils/extractRegion/extractRegion.py ./combined_genomes.fasta -e 10000
#     done < regions.txt  > ${gene}_regions.fa
#     # use riboScan and riboSelect to annotate and group them into operons
#     ribo scan ${gene}_regions.fa -o ./scan_${gene}/ --name ${gene} -v 1
#     ribo select ./scan_${gene}/scannedScaffolds.gb -o ./select_${gene}/ -v 1
#     ribo snag ./scan_${gene}/scannedScaffolds.gb ./select_${gene}/riboSelect_grouped_loci.txt -o ./snag_${gene}/ --msa_tool mafft -v 1 --title "$ ${genus}$ $ ${species}$ ${gene}"
#     cd ../
# done < ./entropy_manifest.tab

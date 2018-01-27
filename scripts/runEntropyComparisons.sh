#!/bin/bash

# if [ -d "assembly_statistics.txt" ]
# then
#     wget -o bacteria.txt "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
#     wget -o archea.txt "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt"
#     cat bacteria.txt > assembly_summary_raw.txt
#     # get rid of the headers for archea for better combining of files
#     cat assembly_summary_archaea.txt | grep "^#" --invert-match >> assembly_summary.txt
# fi

set -e
OLD_IFS=$IFS
IFS=$'\t'

while read refacc refother genus species strain gene geneacc phylum order class family sra compstrain
do
    echo "Downloading $genus $species $strain genomes... and reference $refacc, $refother"
done  < ./entropy_manifest.tab


while read refacc refother genus species strain gene geneacc phylum order class family sra compstrain
do
    thisdir="${genus}_${species}_${strain}"
    if [ -d "$thisdir" ]
    then
	continue
    fi
    mkdir $thisdir
    cd $thisdir
    # get our reference protein sequence for the gene neighboring the first rDNA
    curl -o ${geneacc}.fasta  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$geneacc&rettype=fasta&retmode=text"
    # fetch 25 other genomes to compare ours to
    echo "Running R script"
    Rscript ../../scripts/getCompleteGenomeSubset.R ../assembly_summary.txt genome_selection $genus $species
    # unzip all the new genomes
    gunzip genome_selection/genomes/*
    # get the genome of interest as well
    mkdir ref
    get_genomes.py -q $refacc -o ./ref/
    IFS=","

    if [ "$refother" == "NA" ]
    then
	echo "No additional reference chromosomes."
    else
	for refgenome in $refother
	do
	    get_genomes.py -q $refgenome -o ./ref/
	done
    fi

    IFS=$'\t'
    cd ../
done < ./entropy_manifest.tab



while read refacc refother genus species strain gene geneacc phylum order class family sra compstrain
do
    echo "Processing $genus $species $strain $geneacc ..."
    thisdir="${genus}_${species}_${strain}"
    # mkdir $thisdir
    # lets skip this one if we have already worked on it
    cd $thisdir
    if [ -d "simpleOrtho" ]
    then
	cd ../
	continue
    fi
    # get the interesting regions
    ~/GitHub/open_utils/orthoML/simpleOrtho.py ./${geneacc}.fasta ./genome_selection/genomes/ > regions.txt
    # select only the unique ones
    sort regions.txt | uniq  > unique_regions.txt

    # combine all genomes to create a single file
    cat ./genome_selection/genomes/*.fna > combined_genomes.fasta
    # use extractregions.py to pull out those regions, butwith 10kb on each side.  this all ends up in one file with each of the 25 loci
    while read line
    do
	echo $line | ~/GitHub/open_utils/extractRegion/extractRegion.py ./combined_genomes.fasta -e 10000 -v 2
    done < unique_regions.txt > ${gene}_regions.fa
    # how many strains are we processing?
    files=(./genome_selection/genomes/*.fna)
    count="${#files[@]}"
    # count=$(find  ./genome_selection/genomes/ -type f -name '*.fna')
    # use riboScan and riboSelect to annotate and group them into operons
    ribo scan ${gene}_regions.fa -o ./scan_${gene}/ --name ${gene} -v 2 --min_length 200
    ribo select ./scan_${gene}/scannedScaffolds.gb -o ./select_${gene}/ -v 2
    title="$ ${genus} $ $ ${species} $ rDNAs near ${gene} (N=${count})"
    echo "$title"
    ribo snag ./scan_${gene}/scannedScaffolds.gb ./select_${gene}/riboSelect_grouped_loci.txt -o ./snag_${gene}/ --msa_tool mafft -v 3 --title "$title" --skip_kmers --skip_blast

    # dont forget to process our genome of interest!
    ribo scan ./ref/ -o ./scan_${refacc}/ --name ${refacc} -v 2 --min_length 200
    ribo select ./scan_${refacc}/scannedScaffolds.gb -o ./select_${refacc}/ -v 2
    title="$ ${genus} $ $ ${species} $ ${strain}"
    echo "$title"
    ribo snag ./scan_${refacc}/scannedScaffolds.gb ./select_${refacc}/riboSelect_grouped_loci.txt -o ./snag_${refacc}/ --msa_tool mafft -v 3 --title "rDNAs within ${title}" --skip_kmers --skip_blast
    cd ../
done < ./entropy_manifest.tab


# copy and rename results:
mkdir results_figs
counter=1
while read refacc refother genus species strain gene geneacc phylum order class family sra compstrain
do
    echo "Processing $genus $species $strain $geneacc ..."
    thisdir="${genus}_${species}_${strain}"
    if [ -d "${thisdir}/snag_${gene}" ]
    then

	cp ${thisdir}/snag_${refacc}/entropy_plot.png results_figs/${counter}_genome.png
	cp ${thisdir}/snag_${gene}/entropy_plot.png results_figs/${counter}_gene.png
	counter=$((counter + 1))
    else
	echo "no results found"
    fi
done < ./entropy_manifest.tab

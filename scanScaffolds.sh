#!/bin/bash
# version 0.0.1
# Many scaffolded assemblies from NCBI will not have rRNA annotated.  In order make things work well with riboSeed, point this script at a dirctory into which you have all the contigs.
#
# argument 1 is directory containing contig fastas with the prefix given by arg 2
#    ( must end with sep "/")
# argument 2 is suffix, such as .fa or .fna, etc
# argument 3 is the desired output directory
# output scanScaffolds_combined.gb in current directory
echo 'USAGE: /path/to/contigs/dir/ *ext /path/to/outdir/'
NFILES=$(ls $1*$2  -1 | wc -l) # count files that you will process
THISFILE=1 # counter to increment
mkdir "$3"  # make destination directory
for i in "$1"*"$2";  # for each file matching extension
do
BASENAME=$(basename "$i")  # get the basename of the file without path
# outputs to stderr to make progress easier to see when piping output to log
(>&2 echo "Processing ${BASENAME}, item ${THISFILE} out of ${NFILES}") # status
#echo "Processing ${BASENAME}, item ${THISFILE} out of ${NFILES}"
outdir="$3"${BASENAME}_barrnap
sed 's/^[^ ]*[|]\([^|]*\)[|] .*$/>\1/' ${i} > ${outdir}_renamed${2}

~/barrnap/bin/barrnap -kingdom euk ${outdir}_renamed${2} --reject 0.1 > ${outdir}.gff
# add dumb locus tags
LOCUS=1
while read j; do
    newline=`echo $j | sed -e "s/Name=/locus\_tag=${BASENAME}${LOCUS}\;Name=/"`
    for k in {1..8}; do # for each of the first 8 spaces, convert to tab
    	newline=`echo -e "$newline" | sed -e "s/ /\t/"`
    done
    echo -e "$newline" >> ${outdir}_renamed.gff # -e allows for escaping tabs
    LOCUS=$((LOCUS+1))
done < ${outdir}.gff


# merge fasta and gff3 back to a genbank
seqret -sequence ${outdir}_renamed${2} -feature -fformat gff3 -fopenfile ${outdir}_renamed.gff -osformat genbank -auto  -outseq ${outdir}.gb
THISFILE=$((THISFILE+1))
done
# combine results
cat "$3"*.gb > scanScaffolds_combined.gb

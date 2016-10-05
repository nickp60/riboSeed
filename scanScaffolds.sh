#!/bin/bash
# version 0.0.1
# Requires barrnap and seqret in PATH

# Many assemblies, genomes from NCBI will not have rRNA annotated.
#In order make things work well with riboSeed,
#point this script at a dirctory into which you have all the contigs.
#
# argument 1 is directory containing contig fastas with the prefix given by arg 2
#    ( must end with sep "/")
# argument 2 is suffix, such as .fa or .fna, etc
# argument 3 is the desired output directory
# argument 4 is kingdom: bac  euk mito arc
# argument 5 is threshold [float 0-1], where 1 is 100% identity
# output scanScaffolds_combined.gb in current directory
echo 'USAGE: /path/to/contigs/dir/ *ext /path/to/outdir/ kingdom threshold'
echo 'example: $ barrnap'
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3"]
then
    echo "mandatory arguments: dir, extension, and output_dir"
    exit 1
fi

if [ -d "$3" ]; then
    echo "output dir exists!"
    echo "play things safe, make a new directory for results"
    exit 1
fi
NFILES=$(ls $1*$2  -1 | wc -l) # count files that you will process
## check args
if "$4"!="euk" and "$4"!="bac"
then
    echo "invalid kindom argument! Cannot give threshold without a kingdom"
    exit 1
fi

###
### handle optional arguments
if [ -z "$5" ] and [ -z "$4" ]
then
    echo "using default threshold of .5 identity"
    THRESH=.5
    echo "using default kingdom: euk"
    KINGDOM="euk"
elif [ -z "$5" ]
then
    echo "using default threshold of .5 identity"
    KINGDOM="$4"
    THRESH=.5
else
    KINGDOM="$4"
    THRESH="$5"
fi
###
###

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

barrnap -kingdom "$KINGDOM" ${outdir}_renamed${2} --reject "$THRESH" > ${outdir}.gff
# add dumb locus tags
LOCUS=0
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

cat "$3"*.gb > ${3}scanScaffolds_combined.gb

#!/bin/bash
# version 0.0.4
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

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
then
    echo "mandatory arguments: dir, extension, and output_dir"
    exit 1
fi

if [ -d "$3" ]; then
    echo "output dir exists!"
    echo "play things safe, make a new directory for results"
    exit 1
fi
NFILES=$( ls -1 $1*$2 | wc -l) # count files that you will process

## check args
if [ "$4" != "euk" ] && [ "$4" != "bac" ]
then
    if [ -z "$4" ]
    then
	echo "no kingdom or threshold arguments; using defaults"
    else
	echo "invalid kindom argument! Must be either 'bac' or 'euk'. Cannot give threshold without a kingdom"
	exit 1
    fi
fi

###
### handle optional arguments
if [ -z "$5" ] && [ -z "$4" ]
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
    #  verify threshold is float between 0 and 1
    if [ $(bc <<< "0.0 < $5" ) -eq 1 ] && [ $(bc <<< "$5 < 1.0" ) -eq 1 ]
    then
	KINGDOM="$4"
	THRESH="$5"
    else
        echo "Threshold must be between 0 and 1!"
        exit 1
    fi
fi
###
###

THISFILE=1 # counter to increment
mkdir "$3"  # make destination directory


###############################################################################
for i in "$1"*"$2";  # for each file matching extension
do
BASENAME=$(basename "$i")  # get the basename of the file without path
# outputs to stderr to make progress easier to see when piping output to log
(>&2 echo "Processing ${BASENAME}, item ${THISFILE} out of ${NFILES}") # status
#echo "Processing ${BASENAME}, item ${THISFILE} out of ${NFILES}"
outdir="$3"${BASENAME}_barrnap
sed 's/^[^ ]*[|]\([^|]*\)[|] .*$/>\1/' ${i} > ${outdir}_renamed${2}

#get accession name
ACCNAME=$(grep "^>" ${outdir}_renamed${2} |sed -e 's/>//' -e 's/[[:blank:]].*$//' )

barrnap -kingdom "$KINGDOM" ${outdir}_renamed${2} --reject "$THRESH" > ${outdir}.gff
echo "making locus tags and adding to .gff"
# add dumb locus tags
LOCUS=0
while read j; do
    newline=`echo $j | sed -e "s/Name=/locus\_tag=${BASENAME}${LOCUS}\;Name=/"`
    for k in {1..8}; do # for each of the first 8 spaces, convert to tab
    	newline=`echo -e "$newline" | sed -e "s/ /$(printf '\t')/"`
    done
    echo -e "$newline" >> ${outdir}_renamed.gff # -e allows for escaping tabs
    LOCUS=$((LOCUS+1))
done < ${outdir}.gff

OUTGB=$(echo ${outdir}.gb)

echo "merge fasta and gff3 back to a genbank"
cmd="seqret -sequence ${outdir}_renamed${2} -feature -fformat gff3 \
            -fopenfile ${outdir}_renamed.gff -osformat genbank -auto \
            -outseq ${OUTGB}"
echo ${cmd}
${cmd}

# add version and accession cause biopython complains otherwise
# fixed mac bug when overwriting by adding .bu to i
echo "adding the VERSION and ACCESSION tags to .gb"
if [ `uname` == Darwin ]
then
    echo "using OSX sed"
    sed -i'.bu' -e "1 a\\
    VERSION     $ACCNAME" $OUTGB
    sed -i'.bu' -e "1 a\\
    ACCESSION   $ACCNAME" $OUTGB
elif [ `uname` == Linux ]
then
    echo "using Linux sed"
    sed -i.bu "1 aVERSION     $ACCNAME" $OUTGB
    sed -i.bu "1 aACCESSION   $ACCNAME" $OUTGB
else
    echo "Unsupported OS! Only made for linux and Uni Exiting..."
    exit 1
fi

if [ $THISFILE == 1 ];
then # write
cat $OUTGB > ${3}scannedScaffolds.gb
else # append
cat $OUTGB >> ${3}scannedScaffolds.gb
fi

THISFILE=$((THISFILE+1))
###############################################################################
done

echo "Combined $((THISFILE-1)) gb file(s)"
echo "Results gb file: ${3}scannedScaffolds.gb"

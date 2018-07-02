#!/bin/bash
set -e
#  Takes 5 mandatory args -
# -e experiement name
# -f Forward Reads
# -r Reverse Reads
# -n number of strains
# -o organism name

# And 2 optional
# -p prokaryotes.txt file
# -g directory of genomes of interest

# hacky way to get the path tho the scripts dir
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"


PROKFILE=""
NSTRAINS=""
GENOMESDIR="./genomes_for_run_ANI/"
USAGE="USAGE:run_ani.sh  -e experiment_name -o 'Organism name' -n 5 -f path/to/reads_f.fastq -r path/to/reads_r.fastq"

# make sure we have enough args
if [ -z $5 ]
then
    echo "Not enough arguments"
    echo $USAGE
    exit 1
fi

while getopts "e:o:n:f:r:p:g:" opt; do
    case $opt in
	e)
	    NAME=$OPTARG
	    echo "Experiment name: $OPTARG" >&2
	    ;;
	o)
	    ORGNAME=$OPTARG
	    echo "Organism name: $OPTARG" >&2
	    ;;
	n)
	    NSTRAINS=$OPTARG
	    echo "number of strains: $NSTRAINS" >&2
	    ;;
	f)
	    FREADS=$OPTARG
	    echo "file path F reads: $FREADS" >&2
	    ;;
	r)
	    RREADS=$OPTARG
	    echo "file path R reads: $RREADS" >&2
	    ;;
	p)
	    PROKFILE=$OPTARG
	    echo "path to prokaryotes.txt: $PROKFILE" >&2
	    ;;
	g)
	    GENOMESDIR=$OPTARG
	    echo "Genomes dir: $GENOMESDIR" >&2
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    echo $USAGE
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    exit 1

  esac
done
OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_ANI_${NAME}/"
MINIDIR=${OUTDIR}/mini_assembly/
LOGFILE="${OUTDIR}/log.txt"
mkdir ${OUTDIR}
mkdir ${MINIDIR}

#############################   Run Mini Assembly  #########################
echo "running mini assembly"  2>> "${LOGFILE}"
${SCRIPTPATH}/select_ref_by_ANI/mini_assembly.sh $FREADS $RREADS $MINIDIR  &> "${LOGFILE}"



if [ ! -d "$GENOMESDIR" ] && [ -z "$(ls -A $GENOMESDIR)" ]
then
    #GENOMESDIR="${OUTDIR}/genomes_for_run_ANI"
    mkdir $GENOMESDIR
    ###################   Get potenetial close genomes  #########################
    echo "Get potenetial close genomes"  2>> "${LOGFILE}"
    if [ ! -f "$PROKFILE" ]
    then
        PROKFILE="./prokaryotes.txt"  # default, so we dont pass empty args below
    fi
    ${SCRIPTPATH}/select_ref_by_ANI/get_n_random_complete_genomes.sh -o "$ORGNAME" -n $NSTRAINS -f $PROKFILE > ${OUTDIR}/accessions 2>> "${LOGFILE}"


    #######################   Download close genomes  ###########################
    echo "Downloading genomes"  2>> "${LOGFILE}"

    while read accession
    do
	get_genomes.py -q $accession -o $GENOMESDIR   2>> "${LOGFILE}"
    done < $OUTDIR/accessions
else
    echo "using existing genomes directory"   2>> "${LOGFILE}"
    # delete any existing "conigs.fasta" file from the dir, as those would be from old runs.
    if [ -f "${GENOMESDIR}/contigs.fasta" ]
    then
	rm "${GENOMESDIR}/contigs.fasta"
    fi
fi

echo "copy the mini_assembly result to the potential genomes dir"  2>> "${LOGFILE}"
cp $MINIDIR/spades/contigs.fasta $GENOMESDIR

##########################   Run ANI analysis  ###############################
echo "Running pyani"  >&2
average_nucleotide_identity.py -v -i $GENOMESDIR -g -o $OUTDIR/pyani  2>> "${LOGFILE}"


############# remove contigs from genomes dir if we plan on reusing   ########
rm $GENOMESDIR/contigs.fasta

# extract best hit
echo "extract best hit"  2>> "${LOGFILE}"

python ${SCRIPTPATH}/select_ref_by_ANI/parse_closest_ANI.py $OUTDIR/pyani/ANIm_percentage_identity.tab > ${OUTDIR}/best_reference  2>> "${LOGFILE}"

cat ${OUTDIR}/best_reference

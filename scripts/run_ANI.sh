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

# make sure we have enough args (10, because we have 5 flagged arguments)

if [ "$#" -lt 10 ]
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
PYANILOGFILE="${OUTDIR}/pyanilog.txt"
mkdir ${OUTDIR}
mkdir ${MINIDIR}

#############################   Run Mini Assembly  #########################
echo "running mini assembly"  &>> "${LOGFILE}"
${SCRIPTPATH}/select_ref_by_ANI/mini_assembly.sh $FREADS $RREADS $MINIDIR  &>> "${LOGFILE}"



if [ ! -d "$GENOMESDIR" ] && [ -z "$(ls -A $GENOMESDIR)" ]
then
    #GENOMESDIR="${OUTDIR}/genomes_for_run_ANI"
    mkdir $GENOMESDIR
    ###################   Get potenetial close genomes  #########################
    echo "Get potenetial close genomes"  &>> "${LOGFILE}"
    if [ ! -f "$PROKFILE" ]
    then
        PROKFILE="./prokaryotes.txt"  # default, so we dont pass empty args below
    fi
    ${SCRIPTPATH}/select_ref_by_ANI/get_n_random_complete_genomes.sh -o "$ORGNAME" -n $NSTRAINS -f $PROKFILE > ${OUTDIR}/accessions 2>> "${LOGFILE}"


    #######################   Download close genomes  ###########################
    echo "Downloading genomes"  &>> "${LOGFILE}"

    while read accession
    do
	get_genomes.py -q $accession -o $GENOMESDIR   &>> "${LOGFILE}"
    done < $OUTDIR/accessions
else
    echo "using existing genomes directory"   &>> "${LOGFILE}"
    # delete any existing "conigs.fasta" file from the dir, as those would be from old runs.
    if [ -f "${GENOMESDIR}/contigs.fasta" ]
    then
	rm "${GENOMESDIR}/contigs.fasta"
    fi
fi

echo "copy the mini_assembly result to the potential genomes dir"  &>> "${LOGFILE}"
cp $MINIDIR/spades/contigs.fasta $GENOMESDIR

##########################   Run ANI analysis  ###############################
echo "Running pyani"   &>> "${LOGFILE}"
pyani index $GENOMESDIR -v -l $PYANILOGFILE
# then fix the columns; we dont care about the name
paste <(cut $GENOMESDIR/classes.txt -f 1,2) <(cut $GENOMESDIR/classes.txt -f 2) > tmp_classes
mv tmp_classes $GENOMESDIR/classes.txt
paste <(cut $GENOMESDIR/labels.txt -f 1,2) <(cut $GENOMESDIR/labels.txt -f 2) > tmp_labels
mv tmp_labels $GENOMESDIR/labels.txt


if [ -d .pyani ]
then
    echo "pyani db already exists" &>> $PYANILOGFILE
else
    pyani createdb
fi

pyani anim $GENOMESDIR $OUTDIR/pyani --workers 4 -v -l $PYANILOGFILE --labels $GENOMESDIR/labels.txt --classes $GENOMESDIR/classes.txt

pyani report -v -l $PYANILOGFILE --runs $OUTDIR/pyani

# get the most recent run
this_run=$(tail -n 1 $OUTDIR/pyani/runs.tab | cut -f 1)

pyani report -v -l $PYANILOGFILE --run_matrices $this_run $OUTDIR/pyani


# average_nucleotide_identity.py -v -i $GENOMESDIR -g -o $OUTDIR/pyani  &>> "${LOGFILE}"


############# remove contigs from genomes dir if we plan on reusing   ########
rm $GENOMESDIR/contigs.fasta

# extract best hit
echo "extract best hit"  &>> "${LOGFILE}"

python ${SCRIPTPATH}/select_ref_by_ANI/parse_closest_ANI.py $OUTDIR/pyani/matrix_identity_${this_run}.tab > ${OUTDIR}/best_reference  &>> "${LOGFILE}"

bestref=$(cat ${OUTDIR}/best_reference)
echo -e "${NAME}\t${bestref}"

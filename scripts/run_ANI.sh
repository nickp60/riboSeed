#!/bin/bash
set -e
#  Takes 5 mandatory args -
# -e experiement name
# -f Forward Reads
# -r Reverse Reads
# -n number of strains
# -o organism name

# And 3 optional
# -p prokaryotes.txt file
# -g directory of genomes of interest
# -a existing assembly; skip the miniasembly
# hacky way to get the path tho the scripts dir
for prog in skesa pyani nucmer
do
    if hash $prog 2>/dev/null; then
        echo "$prog installed"
    else
        echo "$prog not installed! exiting"
	exit 1
    fi
done

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"


PROKFILE=""
NSTRAINS=""
ASSEMBLY=""
ASSEMBLER="skesa"
RREADS=""
GENOMESDIR="./genomes_for_run_ANI/"
USAGE="USAGE:run_ani.sh  -e experiment_name -o 'Organism name' -n 5 -f path/to/reads_f.fastq -r path/to/reads_r.fastq"

# make sure we have enough args (10, because we have 5 flagged arguments)

if [ "$#" -lt 8 ]
then
    echo "Not enough arguments"
    echo $USAGE
    exit 1
fi

while getopts "e:o:n:f:r:p:g:a:s:" opt; do
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
	a)
	    ASSEMBLY=$OPTARG
	    echo "assembly: $ASSEMBLY" >&2
	    ;;
	s)
	    ASSEMBLER=$OPTARG
	    case "$ASSEMBLER" in
		skesa)
		    # Do stuff
		    ;;
		spades)

		    # Do stuff
		    ;;
		*)
		    echo "Error: assembler must be either spades or skesa"
		    exit 1
		    ;;
	    esac
	    echo "assembler: $ASSEMBLER" >&2
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
if [ -s $ASSEMBLY ]
then
    echo "Downsampling reads"  >> "${LOGFILE}"

    mkdir ${OUTDIR}/downsampled_reads/

    {
	echo "seqtk sample -s100 $FREADS 100000 > ${OUTDIR}/downsampled_reads/reads1.fq" >&2
	seqtk sample -s100 $FREADS 100000 > ${OUTDIR}/downsampled_reads/reads1.fq
    } || {
	echo "seqtk sample -s100 $FREADS .1    > ${OUTDIR}/downsampled_reads/reads1.fq" >&2
	seqtk sample -s100 $FREADS .1    > ${OUTDIR}/downsampled_reads/reads1.fq
    }
    if [ ! -s $RREADS ]
    then
	{
	    echo "seqtk sample -s100 $RREADS 100000   > ${OUTDIR}/downsampled_reads/reads1.fq" >&2
	    seqtk sample -s100 $RREADS 100000 > ${OUTDIR}/downsampled_reads/reads2.fq
	} || {
	    echo "seqtk sample -s100 $RREADS .1   > ${OUTDIR}/downsampled_reads/reads1.fq" >&2
	    seqtk sample -s100 $RREADS .1    > ${OUTDIR}/downsampled_reads/reads2.fq
	}
	SPADES_STRING="-1 ${OUTDIR}/downsampled_reads/reads1.fq -2 ${OUTDIR}/downsampled_reads/reads2.fq"
        SKESA_STRING="${OUTDIR}/downsampled_reads/reads1.fq ${OUTDIR}/downsampled_reads/reads2.fq"
    else
	SPADES_STRING="-s ${OUTDIR}/downsampled_reads/reads1.fq"
        SKESA_STRING="${OUTDIR}/downsampled_reads/reads1.fq"
    fi
    echo "running mini assembly"  >> "${LOGFILE}"
    echo "$SPADES_string"
    mkdir ${MINIDIR}/assembly/
    case "$ASSEMBLER" in
	skesa)
	    skesa --fastq $SKESA_STRING --contigs_out ${MINIDIR}/assembly/contigs.fasta
	    ;;
	spades)
	    spades.py $SPADES_STRING -o ${MINIDIR}/assembly/ -t 4 -m 64
	    ;;
    esac

    ASSEMBLY=${MINIDIR}/assembly/contigs.fasta
fi

if [ ! -d "$GENOMESDIR" ]
then
    #GENOMESDIR="${OUTDIR}/genomes_for_run_ANI"
    mkdir $GENOMESDIR
fi

if  [ -z "$(ls -A $GENOMESDIR)" ]
then
    ###################   Get potenetial close genomes  #########################
    echo "Get potenetial close genomes"  >> "${LOGFILE}"
    if [ ! -s "$PROKFILE" ]
    then
        PROKFILE="./prokaryotes.txt"  # default, so we dont pass empty args below
    fi
    echo "${SCRIPTPATH}/select_ref_by_ANI/get_n_random_complete_genomes.sh -o \"$ORGNAME\" -n $NSTRAINS -f $PROKFILE > ${OUTDIR}/accessions "
    ${SCRIPTPATH}/select_ref_by_ANI/get_n_random_complete_genomes.sh -o "$ORGNAME" -n $NSTRAINS -f $PROKFILE > ${OUTDIR}/accessions


    #######################   Download close genomes  ###########################
    echo "Downloading genomes"  >> "${LOGFILE}"

    while read accession
    do
	get_genomes.py -q $accession -o $GENOMESDIR   >> "${LOGFILE}"2>&1
    done < $OUTDIR/accessions
else
    echo "using existing genomes directory"   >> "${LOGFILE}"
    # delete any existing "conigs.fasta" file from the dir, as those would be from old runs.
    # probaly shoudl make this a try/catch thingy
    # {
    # 	rm "${GENOMESDIR}/contigs_*.*"
    # } || {
    # 	echo "no need to clean genomes dir" &>> "${LOGFILE}"
    # }
fi

echo "copy the mini_assembly result to the potential genomes dir"  >> "${LOGFILE}"
cp $ASSEMBLY $GENOMESDIR/contigs_${NAME}.fasta

##########################   Run ANI analysis  ###############################
echo "Running pyani index"   >> "${LOGFILE}"
if [ -f "$GENOMESDIR/labels.txt" ]
then
    rm $GENOMESDIR/labels.txt
fi
if [ -f "$GENOMESDIR/classes.txt" ]
then
    rm $GENOMESDIR/classes.txt
fi
echo "pyani index $GENOMESDIR -v -l $PYANILOGFILE"  >&2
pyani index $GENOMESDIR -v -l $PYANILOGFILE >> "${LOGFILE}" 2>&1
echo "recreating pyani labels"   >> "${LOGFILE}"
# then fix the columns; we dont care about the name
paste <(cut -f 1,2 $GENOMESDIR/classes.txt) <(cut -f 2 $GENOMESDIR/classes.txt) > tmp_classes
mv tmp_classes $GENOMESDIR/classes.txt
paste <(cut -f 1,2 $GENOMESDIR/labels.txt) <(cut -f 2 $GENOMESDIR/labels.txt) > tmp_labels
mv tmp_labels $GENOMESDIR/labels.txt


if [ -d .pyani ]
then
    echo "pyani db already exists" >> $PYANILOGFILE
else
    echo "pyani createdb"  >&2
    pyani createdb >> "${LOGFILE}" 2>&1
fi

echo "Running pyani anim"   >> "${LOGFILE}"
echo "pyani anim $GENOMESDIR $OUTDIR/pyani --workers 64  -v -l $PYANILOGFILE --labels $GENOMESDIR/labels.txt --classes $GENOMESDIR/classes.txt"  >&2
pyani anim $GENOMESDIR $OUTDIR/pyani --workers 64  -v -l $PYANILOGFILE --labels $GENOMESDIR/labels.txt --classes $GENOMESDIR/classes.txt >> "${LOGFILE}" 2>&1

echo "Generating report of pyani runs"   >> "${LOGFILE}"
echo "pyani report -v -l $PYANILOGFILE --runs $OUTDIR/pyani "  >&2
pyani report -v -l $PYANILOGFILE --runs $OUTDIR/pyani  >> "${LOGFILE}" 2>&1


echo "Getting name of most recent pyani run"   >> "${LOGFILE}"
this_run=$(tail -n 1 $OUTDIR/pyani/runs.tab | cut -f 1)

echo "Generating genome comparison matrix"   >> "${LOGFILE}"
echo "pyani report -v -l $PYANILOGFILE --run_matrices $this_run $OUTDIR/pyani "  >&2
pyani report -v -l $PYANILOGFILE --run_matrices $this_run $OUTDIR/pyani >> "${LOGFILE}" 2>&1


# average_nucleotide_identity.py -v -i $GENOMESDIR -g -o $OUTDIR/pyani  &>> "${LOGFILE}"


############# remove contigs from genomes dir if we plan on reusing   ########
# rm $GENOMESDIR/contigs.fasta
rm $GENOMESDIR/contigs_${NAME}.fasta
rm $GENOMESDIR/contigs_${NAME}.md5
# rm $GENOMESDIR/labels.txt
# rm $GENOMESDIR/classes.txt
# extract best hit
echo "extract best hit"  >> "${LOGFILE}"
echo "python ${SCRIPTPATH}/select_ref_by_ANI/parse_closest_ANI.py $OUTDIR/pyani/matrix_identity_${this_run}.tab > ${OUTDIR}/best_reference"  >&2
python ${SCRIPTPATH}/select_ref_by_ANI/parse_closest_ANI.py $OUTDIR/pyani/matrix_identity_${this_run}.tab > ${OUTDIR}/best_reference

bestref=$(cat ${OUTDIR}/best_reference)
# this should be the only line going to stdout
echo -e "${NAME}\t${bestref}"

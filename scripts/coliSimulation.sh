#!/bin/bash
# version 0.0.4
#  If submitting jobs to a cluster with SGE, set the variables
#  below and the environment variables starting with #$.  Otherwise, you
#  can submit arguments to commandline
#$ -cwd
#$ -j yes
#$ -V
#$ -N simKleb
#$ -pe mpi 8

if [ -z "$1" ]
then
    echo "One required argument: replicate number. \n USAGE: coliSimulation.sh <replicate#>"
    exit 1
else
    REP="$1"
fi

OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_simulation_${REP}/"
mkdir ${OUTDIR}


#source ${HOME}/Virtualenvs/riboSeed/bin/activate
set -euo pipefail
#IFS=$'\n\t'

# get_genomes.py -q BA000007.2 -o ./
if [ -d "${OUTDIR}BA000007.2.fasta" ]
then
    echo "using local copy of genome" 
else
    echo "fetching the genome"
    get_genomes.py -q BA000007.2 -o ${OUTDIR}
fi
LOG=${OUTDIR}log.txt
echo "Running scan and select on reference"

ribo scan ${OUTDIR}BA000007.2.fasta -o ${OUTDIR}scan/ -n Sakai &>> ${LOG}
ribo select ${OUTDIR}scan/scannedScaffolds.gb -o ${OUTDIR}select/ &>> ${LOG}

runs=(
"GA2 50"
 "GA2 75"
 "HS20 100"
 # "HS25 150"
 # "MSv3 250"
)

##############
echo "fetching MG1655 to make reads from"
get_genomes.py -q NC_000913.3 -o ${OUTDIR}

##############




SEED=27
mkdir ${OUTDIR}mauve
for FLANK in 1000
do
for COV in 20 50
do
for ((i = 0; i < ${#runs[@]}; i++))
do
    echo "Running simulation for ${runs[$i]} with ${COV}x coverage, rep ${REP}"    
    runarray=(${runs[$i]})
    tech=${runarray[0]}
    len=${runarray[1]}

    name="sim_${tech}_${len}_COV_${COV}_${FLANK}_rep${REP}"
    thisoutdir="${OUTDIR}${name}/"
    mkdir ${thisoutdir}

    echo " simulating reads from the mg1655 genome with seed ${REP}"

    ~/bin/pirs-2.0.2/pirs simulate -m 300 --threads 1 -l ${len} -x ${COV} -v 10 -B ~/bin/pirs-2.0.2/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz -I ~/bin/pirs-2.0.2/Profiles/InDel_Profiles/phixv2.InDel.matrix -G ~/bin/pirs-2.0.2/Profiles/GC-depth_Profiles/humNew.gcdep_200.dat -o ${thisoutdir}MG1655_${tech} --compress -S ${REP} ${OUTDIR}NC_000913.3.fasta &>> ${LOG}

    echo " running riboSeed with sakai as reference"
    ribo seed -r ${OUTDIR}scan/scannedScaffolds.gb -o ${thisoutdir}seed/ -F ${thisoutdir}MG1655_${tech}_${len}_300_1.fq.gz -R ${thisoutdir}MG1655_${tech}_${len}_300_2.fq.gz  -v 1 -c 4 -t 2 -z ${OUTDIR}select/riboSelect_grouped_loci.txt -s $((len/2)) --memory 16 --clean_temps -l ${FLANK} &>> ${LOG}

    ribo score ${thisoutdir}seed/mauve &>> ${LOG}
    rm ${thisoutdir}seed -rf

done
done
done    



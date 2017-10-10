#!/bin/bash
# version 0.0.5
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
    echo "One required argument: replicate number. \n USAGE: klebSimulation.sh <replicate#>"
    exit 1
else
    REP="$1"
fi

OUTDIRBASE="./"
OUTDIR="${OUTDIRBASE}`date +%F`_simulation_${REP}/"
mkdir ${OUTDIR}


set -euo pipefail
#IFS=$'\n\t'

# download HS11286
if [ -d "${OUTDIR}CP003200.1.fasta" ]
then
    echo "using local copy of genome" 
else
    echo "fetching the genome"
    get_genomes.py -q CP003200.1 -o ${OUTDIR}
fi
LOG=${OUTDIR}log.txt
echo "Running scan and select on reference"
ribo scan ${OUTDIR}CP003200.1.fasta -o ${OUTDIR}scan/ -n HS11286 &>> ${LOG}
ribo select ${OUTDIR}scan/scannedScaffolds.gb -o ${OUTDIR}select/ &>> ${LOG}

runs=(
  "GA2 50"
  "GA2 75"
 "HS20 100"
 # "HS25 150"
 # "MSv3 250"
)

##############
echo "fetching NC_012731.1 (NTUH-K2044) to generate reads from"
get_genomes.py -q NC_012731.1 -o ${OUTDIR}

##############




SEED=27
mkdir ${OUTDIR}mauve
for FLANK in  1000
do
for COV in 20 50
do
for ((i = 0; i < ${#runs[@]}; i++))
do
    echo "replicate ${REP}"
    echo "Running simulation for ${runs[$i]}"    
    runarray=(${runs[$i]})
    tech=${runarray[0]}
    len=${runarray[1]}
    
    name="sim_${tech}_${len}_COV_${COV}_${FLANK}_rep${REP}"
    echo "output directory"
    thisoutdir="${OUTDIR}${name}/"
    mkdir ${thisoutdir}

    echo "simulating reads from the kleb genome"

     ~/bin/pirs-2.0.2/pirs simulate -m 300 --threads 1 -l ${len} -x ${COV} -v 10 -B ~/bin/pirs-2.0.2/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz -I ~/bin/pirs-2.0.2/Profiles/InDel_Profiles/phixv2.InDel.matrix -G ~/bin/pirs-2.0.2/Profiles/GC-depth_Profiles/humNew.gcdep_200.dat -o ${thisoutdir}NTUH_${tech} --compress -S ${REP} ${OUTDIR}NC_012731.1.fasta &>> ${LOG}

    echo "running riboSeed with sakai as reference"
    ribo seed -r ${OUTDIR}scan/scannedScaffolds.gb -o ${thisoutdir}seed/ -F ${thisoutdir}NTUH_${tech}_${len}_300_1.fq.gz -R ${thisoutdir}NTUH_${tech}_${len}_300_2.fq.gz  -v 1 -c 4 -t 1 -z ${OUTDIR}select/riboSelect_grouped_loci.txt --memory 32 -s $((len/2)) -l ${FLANK} &>> ${LOG}

    ribo score ${thisoutdir}seed/mauve &>> ${LOG}
    rm ${thisoutdir}seed -rf
done
done
done    




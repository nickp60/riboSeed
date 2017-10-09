#!/bin/bash
set -eo pipefail
# version 0.2.0
IFS=$'\n\t'


# process original
# only need to do this once;
# we arent doing indels so the coords are the same

if [ -z "$1" ]
then
    SEED=17
else
    SEED="$1"
fi

# get the reference Sakai genome
if [ -e BA000007.2.fasta ]
then
echo "using local copy of reference"
else
get_genomes.py -q BA000007.2 -o ./
fi

# make output dirs here
OUTDIRBASE="./"
mkdir ${OUTDIRBASE}/toyGenome_${SEED}/
# lets log all the outputs to a single log file, shall we?
LOGFILE=${OUTDIRBASE}/degenerate_log_${SEED}.txt

echo "Running degenerate pipeline with the following settings"
echo "SEED: ${SEED}"
echo " Annotate regions in full genome"
ribo scan ./BA000007.2.fasta -o ${OUTDIRBASE}/toyGenome_${SEED}/scan/ &>> ${LOGFILE}
echo " Cluster regions"
ribo select ${OUTDIRBASE}/toyGenome_${SEED}/scan/scannedScaffolds.gb -o ${OUTDIRBASE}/toyGenome_${SEED}/select/ &>> ${LOGFILE}
echo " Extract regions with 5kb flanking"
ribo snag ${OUTDIRBASE}/toyGenome_${SEED}/scan/scannedScaffolds.gb  ${OUTDIRBASE}/toyGenome_${SEED}/select/riboSelect_grouped_loci.txt -o ${OUTDIRBASE}/toyGenome_${SEED}/snag/ -l 5000 -v 1 --just_extract &>> ${LOGFILE}
echo " Combine the extracted regions into a toy genome"
python3.6 ~/GitHub/riboSeed/scripts/concatToyGenome.py ${OUTDIRBASE}/toyGenome_${SEED}/snag/ \*_riboSnag.fasta -o ${OUTDIRBASE}/toyGenome_${SEED}/coli_genome/ &>> ${LOGFILE}
echo " Generate reads from the toy genome"

~/bin/pirs-2.0.2/pirs simulate -m 300 -l 100 -x 30 -v 10 -o ${OUTDIRBASE}/toyGenome_${SEED}/reads -B ~/bin/pirs-2.0.2/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz --compress -I ~/bin/pirs-2.0.2/Profiles/InDel_Profiles/phixv2.InDel.matrix -G ~/bin/pirs-2.0.2/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat ${OUTDIRBASE}/toyGenome_${SEED}/coli_genome/concatenated_seq.fasta &>> ${LOGFILE}

# annotate the toy genome for mauve visualization
ribo scan ${OUTDIRBASE}/toyGenome_${SEED}/coli_genome/ -o  ${OUTDIRBASE}/toyGenome_${SEED}/coli_genome/ref_scan/ -v 1 &>> ${LOGFILE}


FREQS=( 0.0 0.0025 0.0050 0.0075 0.0100 0.0150 0.0200 0.0250 0.0500 0.0750 0.1000 0.1250 0.1500 0.1750 0.2000 0.2250 0.2500 0.2750 0.3000 )

# we will mutate either just the flanking regions (FLANK) or the whole reference sequece (ALL)
for WHERE in FLANK ALL
do
    # set the arg for riboSim
    if [ "$WHERE" == "ALL" ]
    then
	ARG=""
    else
	# make sure 5000 matches the ribo sang call above
	ARG="-e 5000" 
    fi

    OUTDIR="${OUTDIRBASE}`date +%F`_degenerate_output_${WHERE}_${SEED}"
    echo " OUTDIR: ${OUTDIR}"

    mkdir ${OUTDIR}
    mkdir ${OUTDIR}/genomes/
    mkdir ${OUTDIR}/mauve/
    mkdir ${OUTDIR}/score_reports/

    # move a copy of the annotated genome to the "mauve" directory for convenience
    cp ${OUTDIRBASE}/toyGenome_${SEED}/coli_genome/ref_scan/scannedScaffolds.gb ${OUTDIR}/mauve/reference.gb

    # cluster regions in toy genome
    ribo select ${OUTDIRBASE}/toyGenome_${SEED}/coli_genome/ref_scan/scannedScaffolds.gb  -o ${OUTDIR}/ref_select/ -v 1  &>> ${LOGFILE}


    LENGTH=$((${#FREQS[@]}))
    echo " Number of frequencies to run: $LENGTH"
    SUBSEEDS=(`seedRand.py $SEED $LENGTH`)
    for ((i=0; i<${#SUBSEEDS[@]}; ++i));
    do

	SUBSEED=${SUBSEEDS[i]}
	FREQ=${FREQS[i]}
	echo "  Processing mutation frequency of ${FREQ}, seeded with ${SUBSEEDS[i]}"
	###################################################################
	# make a new reference
	mkdir ${OUTDIR}/mutated_regions_${FREQ}/

	for region in ${OUTDIRBASE}/toyGenome_${SEED}/snag/*_riboSnag.fasta
	do
	    echo "   Mutating ${region}"
	    # intrduce the mutations
	    fname_ext=${region##*/}
	    fname=${fname_ext%.*}
	    ribo sim ${region} ${ARG}  -f ${FREQ} -o ${OUTDIR}/degenerate_${FREQ}_${fname}/ -v 1 --seed ${SUBSEED} &>> ${LOGFILE}
	    # copy results to a directory
	    # echo "  copying ${FREQ} results to mutated regions folder"
	    cp ${OUTDIR}/degenerate_${FREQ}_${fname}/${fname_ext} ${OUTDIR}/mutated_regions_${FREQ}/
	done

	echo "   combine the mutated extracted regions into a toy genome"
	python3.6 ~/GitHub/riboSeed/scripts/concatToyGenome.py ${OUTDIR}/mutated_regions_${FREQ}/ \*_riboSnag.fasta -o ${OUTDIRBASE}/toyGenome_${SEED}/mutated_genome_${FREQ}_${WHERE}/ &>> ${LOGFILE}
	###################################################################
	# copy the file
	cp ${OUTDIRBASE}/toyGenome_${SEED}/mutated_genome_${FREQ}_${WHERE}/concatenated_seq.fasta ${OUTDIR}/genomes/${FREQ}_reference.fasta
	# generate the gb file for the sequence
	echo "   running riboScan"
	ribo scan ${OUTDIRBASE}/toyGenome_${SEED}/mutated_genome_${FREQ}_${WHERE}/ -o ${OUTDIR}/scan_${FREQ}/ -v 1  &>> ${LOGFILE}
	{
	    # run riboSeed
	    echo "   running riboSeed"
	    ribo seed -r ${OUTDIR}/scan_${FREQ}/scannedScaffolds.gb  -o ${OUTDIR}/seed_${FREQ}/ ${OUTDIR}/ref_select/riboSelect_grouped_loci.txt -F ${OUTDIRBASE}/toyGenome_${SEED}/reads_100_300_1.fq.gz -R ${OUTDIRBASE}/toyGenome_${SEED}/reads_100_300_2.fq.gz -i 3  -v 1 -l 1000 --clean_temps -c 4 --memory 16 --threads 2 --skip_control &>> ${LOGFILE}
	# move a copy of the annotated genome to the "mauve" directory for convenience
	} || {
	    echo "   Error running riboSeed for ${FREQ} ${WHERE}"
	}
	if [ -e ${OUTDIR}/seed_${FREQ}/final_de_fere_novo_assembly/contigs.fasta ]
	then
	    cp ${OUTDIR}/seed_${FREQ}/final_de_fere_novo_assembly/contigs.fasta ${OUTDIR}/mauve/de_fere_novo_${FREQ}.fasta
	    echo "   running riboScore"

	    ribo score ${OUTDIR}/seed_${FREQ}/mauve/ &>> ${LOGFILE}

	    cp ${OUTDIR}/seed_${FREQ}/mauve_riboScored/riboScore_report.txt ${OUTDIR}/score_reports/${SUBSEED}_${FREQ}_riboScore_report.txt

	fi
	rm ${OUTDIR}/seed_${FREQ}/ -rf
    done

    cat ${OUTDIR}/score_reports/*_riboScore_report.txt > ${OUTDIR}/score_reports/combined_reports.txt
done

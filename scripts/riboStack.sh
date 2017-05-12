#!/bin/bash
GFF="$1"
BAM="$2"
FASTA="$3"
set -e
if [ "$#" -ne 3 ]; then
    echo "Requires three arguments: GFF, BAM file, and a FASTA"
    exit 1
fi
base=${FASTA##*/}
echo $base
base=${base%.*}
echo $base
N=`awk '{if (NR!=1) {print}}' ${FASTA} | wc -c`
echo $N
echo -e "$base\t$N" > temp_geno

if [ -x `command -v bedtools` ] 
 then  
#bedtools
 awk '{if (NR!=1) {print $1,$4,$5}}' OFS="\t" $GFF >temp_region
for i in {1..10}
do
bedtools shuffle -i temp_region -g temp_geno > test_region_$i
done
  
 else 
echo "bedtools is not installed.  Exiting"
exit 1 
 fi

if [ -x `command -v samtools` ] 
then
REFERENCE=`samtools view -b -F 256 $BAM | samtools depth -b temp_region -  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'`
results=()
for j in {1..10}
do
res=`samtools view -b -F 256 $BAM | samtools depth -b test_region_$j -  |  awk '{sum+=$3} END { print sum/NR}'`
results+=("$res")
echo "Reference regions: $REFERENCE"
echo "Sampled regions:"
echo "${results[@]}"
done

 else 
echo "samtools is not installed.  Exiting"
exit 1 
 fi


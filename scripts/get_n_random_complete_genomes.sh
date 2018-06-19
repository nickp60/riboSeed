#!/bin/bash
set -oue
#  Takes three args -
# -o the quoted nameof the organism
# -n the number of strains
# -f the path to an existing prokaryotes.txt file
while getopts "o:n:f:" opt; do
    case $opt in
	o)
	    ORGNAME=$OPTARG
	    echo "Organism name: $OPTARG" >&2

	    ;;
	n)
	    NSTRAINS=$OPTARG
	    echo "number of strains: $NSTRAINS" >&2
	    ;;
	f)
	    PROKFILE=$OPTARG
	    echo "file path to ncbi flatfile: $PROKFILE" >&2
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    echo "USAGE: -o 'Organism name' -n 5 -f ./path/to/prokaryotes"
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    exit 1

  esac
done

if [ ! -f "$PROKFILE" ]
then
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
    PROKFILE=./prokaryotes.txt
fi

# column 9 has the nucc accession if it is a complete genome, or a "-" if empt
# it starts with chromasome
# here we select the lines for all the complete genomes with awk,
# find the lines matching out query
# and save a temp file with the results
cat  $PROKFILE | \
    awk -F "\t" '$9 ~ /chrom*/ { print $0 }' |  \
    grep "$ORGNAME" > \
	 ./tmp_raw_outfile

# if file is not empty
if [ ! -s ./tmp_raw_outfile ]
then
    echo "grepping for '$ORGNAME' returned no results"
    exit 1
fi

# now we shuffle the file, get the top n lines, and use some sed to split apart the
# chromosome:NZ_CP013218.1/CP013218.1; plasmid paadD:NZ_CP014695.1/CP014695.1; plasmid plinE154:NZ_CP014694.1/CP014694.1
# to
# NZ_CP013218.1
# Note that we only get the first chromasome for a given entry. Sorry vibrioists

# shuf ./tmp_raw_outfile | head -n $NSTRAINS | cut -d "\t" -f 9
shuf ./tmp_raw_outfile | \
    head -n $NSTRAINS | \
    cut -f 9 | \
    sed "s/chro.*://" | \
    sed "s/\/.*//" > tmp_accessions

mkdir

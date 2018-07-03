#!/bin/bash
set -e
#  Takes three args -
# -o the quoted nameof the organism
# -n the number of strains
# -f the path to an existing prokaryotes.txt file
PROKFILE=""
NSTRAINS=""
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
    if [ ! -f "prokaryotes.txt" ]
    then
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt >&2
    fi

    PROKFILE=./prokaryotes.txt
    if [ ! -s "$PROKFILE" ]
    then
	echo "The file '$PROKFILE' is empty/corrupted; please delete and try again"
	exit 1
    fi
fi

# column 9 has the nucc accession if it is a complete genome, or a "-" if empt
# it starts with chromasome
# here we select the lines for all the complete genomes with awk,
# find the lines matching out query
# and save a temp file with the results

if [ -f  /tmp/prok_subset_raw_outfile ]
then
    rm /tmp/prok_subset_raw_outfile
fi


cat  $PROKFILE | \
    awk -F "\t" '$9 ~ /chrom*/ { print $0 }' |  \
    grep "$ORGNAME" > \
	 /tmp/prok_subset_raw_outfile

# if file is empty, raise an error

if [ ! -s /tmp/prok_subset_raw_outfile ]
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

if [ $(command -v shuf) ]
then
    echo "using shuf" >&2
    SHUF=shuf
else
    echo "using gshuf" >&2
    SHUF=gshuf
fi

echo "selecting $NSTRAINS random strains" >&2

if [ $NSTRAINS != "" ]
then
    $SHUF /tmp/prok_subset_raw_outfile | \
	head -n $NSTRAINS | \
	cut -f 9 | \
	# sed "s/chro.*://" | \
	sed "s/^chro[^:]*://" | \
	# sed "s/\/.*//"
        sed "s/[;\/].*//"	
else
        $SHUF /tmp/prok_subset_raw_outfile | \
	cut -f 9 | \
	# sed "s/chro.*://" | \
	sed "s/^chro[^:]*://" | \
	# sed "s/\/.*//"  # handle instances lacking both genbank and refseq accs
        sed "s/[;\/].*//"	
fi

echo "done" >&2

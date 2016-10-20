# RiboSeed Pipeline

## Description

RiboSeed is an supplemental assembly method to try to address the issue of multiple ribosomal regions in a genome.  It takes advantage of the fact that while each region is identical, the regions flanking are unique, and therfore should be able to be used to seed an alignment.

The pipeline (currently) consists of optional preprocessing and two main stages:

### 0: Preprocessing with scanScaffolds.sh

This is an (ever-growing) shell script to preprocess seqeuencees straight from DNA fasta.  The issue with many legacy annotations, assemblies, and scaffold collections is they are often poorly annotated at best, and unannoated at worst.  This is a quick way to happiness without using the full Prokka annotation path. It requires barrnap [cite] and seqret [cite] to be availible in your path.

#### Usage

Move your data to its own directory, whether it is a single fasta genome or multiple fasta scaffolds.  For each file it finds in this dir with a given extension, the script renames complex headers (sketchy), scans with barrnap and captures the output gff.  It then edits the gff to add fake incrementing locus_tags, and uses the original sequence file through seqret to make a genbank file that contains just annotated rRNA features. the last step is a concatenation which, whether or not there are multiple files, makes a single (possiby multi-entry) genbank file ripe for riboseeding.
* argument 1 is directory containing contig fastas with the prefix given by arg 2
*    ( must end with sep "/")
* argument 2 is suffix, such as .fa or .fna, etc
* argument 3 is the desired output directory
* argument 4 is kingdom: bac  euk mito arc
* argument 5 is threshold [float 0-1], where 1 is 100% identity
* output scanScaffolds_combined.gb in current directory

## 1: Selection and Extraction


For extracting other conserved regions for use in seeded assembly:

* `otherSelect` [under development]


For extracting ribosomal regions for use in seeded assembly

* `riboSelect.py` searches the genome for rRNA annotations, clusters them into likely ribosomal groups, and outputs a colon-separated list of clustered rRNA locus tags by record id.

You will probably want to preview your file to figure out the syntax used. (ie, 16s vs 16S, rRNA vs RRNA, etc...)

For fungal or other Eukaryotic genomes, reannotation is frequently required.  For `riboSelect.py`, you will need to change `--specific_features` appropriatly to 5_8S:18S:28S.

NOTE: the format is very simple, and due to the relatively small number of such coding sequences in bacterial genomes, this can be constructed by hand if the clusters do not look appropiate. The format is "genome_sequence_id locus_tag1:locus_tag2", where each line represents a cluster. See example below, where 14 rRNA's are clustered into 6 groups:

```
CM000577.1 FGSG_20052:FGSG_20051:FGSG_20053
CM000577.1 FGSG_20048:FGSG_20047
CM000577.1 FGSG_20049:FGSG_20050
CM000577.1 FGSG_20054:FGSG_20056:FGSG_20055
CM000577.1 FGSG_20058:FGSG_20057
CM000577.1 FGSG_20075:FGSG_20074
```

#### Usage

```
usage: riboSelect.py [-h] [-o OUTPUT] [-f FEATURE] [-s SPECIFIC_FEATURES]
                     [--keep_temps] [--clobber] [-c CLUSTERS] [-v VERBOSITY]
                     [--debug]
                     genbank_genome

This is used to extract rRNA regions from a gb file, returnsa text file with
the clusters

positional arguments:
  genbank_genome        Genbank file (WITH SEQUENCE)

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  -o OUTPUT, --output OUTPUT
                        output directory; default: cwd

optional arguments:
  -f FEATURE, --feature FEATURE
                        Feature, rRNA or RRNA; default: rRNA
  -s SPECIFIC_FEATURES, --specific_features SPECIFIC_FEATURES
                        colon:separated -- specific features; default:
                        16S:23S:5S
  --keep_temps          view intermediate clustering files; default: False
  --clobber             overwrite previous output files; default: False
  -c CLUSTERS, --clusters CLUSTERS
                        number of rDNA clusters;if submitting multiple
                        records, must be a colon:separated list that matches
                        number of genbank records. Default is inferred from
                        specific feature with fewest hits
  -v VERBOSITY, --verbosity VERBOSITY
                        1 = debug(), 2 = info(), 3 = warning(), 4 = error()
                        and 5 = critical(); default: 2
  --debug               Enable debug messages

```


* `riboSnag.py` takes the list of clustered locus tags and extracts their sequences with flanking regions, optionally turning the coding sequences to N's to minimize bias towards reference. Is used to pull out regions of interest from a Genbank file. Outputs a directory with a fasta file for each clustered region (and a log file).

#### Usage:

```
usage: riboSnag.py [-o OUTPUT] [-f FEATURE] [-w WITHIN]
                   [-m MINIMUM] [-l FLANKING] [-r] [-c] [-p PADDING]
                   [-v VERBOSITY] [--clobber] [-h]
                   genbank_genome clustered_loci

Use to extract regions of interest based on supplied locus tags.

positional arguments:
  genbank_genome        Genbank file (WITH SEQUENCE)
  clustered_loci        output from riboSelect

optional arguments:
  -f FEATURE, --feature FEATURE
                        Feature, such as CDS,tRNA, rRNA; default: rRNA

required named arguments:
  -o OUTPUT, --output OUTPUT
                        output directory; default: cwd

optional arguments:
  -w WITHIN, --within_feature_length WITHIN
                        bp's to include within the region; default: 0
  -m MINIMUM, --minimum_feature_length MINIMUM
                        if --replace, and sequence is shorter than 2x
                        --within_feature_length, --within will be modified so
                        that only -m bp of sequnece areturned to N's default:
                        100
  -l FLANKING, --flanking_length FLANKING
                        length of flanking regions, can be colon-separated to
                        give separate upstream and downstream flanking
                        regions; default: 1000
  -r, --replace         replace sequence with N's; default: False
  -c, --circular        if the genome is known to be circular, and an region
                        of interest (including flanking bits) extends past
                        chromosome end, this extends the seqence past
                        chromosome origin forward by 5kb; default: False
  -p PADDING, --padding PADDING
                        if treating as circular, this controls the length of
                        sequence added to the 5' and 3' ends to allow for
                        selecting regions that cross the chromosom's origin;
                        default: 5000
  -v VERBOSITY, --verbosity VERBOSITY
                        1 = debug(), 2 = info(), 3 = warning(), 4 = error()
                        and 5 = critical(); default: 2
  --clobber             overwrite previous output filesdefault: False
  -h, --help            Displays this help message

```

## 2: Seeded Assebly

* `riboSeed.py` is used to map reads to the extracted regions in an iterative manner, assembling the extracted reads, and then running `SPAdes` assembly to hopefully resolve the contig junctions.

#### Output

This outputs two main directories: `map` and `results`.  If `--temps` is true, temporary files from the mapping scheme will be retained, and is useful for assessing problems.

The results directory will contain a 'mauve' directory with all the extended fragments, the mapped fastq files, and a `de_novo` and `de_fere_novo` folder, containing the results with the *de novo* mapping and supplemented mapping, respectively.

#### Usage:

```
usage: riboSeed.py -F FASTQ1 -R FASTQ2 -r REFERENCE_GENOME -o OUTPUT
                   [-S FASTQS] [-n EXP_NAME] [-m METHOD] [-c CORES] [-k KMERS]
                   [-p PRE_KMERS] [-g MIN_GROWTH] [-s MIN_SCORE_SMALT]
                   [-a MIN_ASSEMBLY_LEN] [--paired_inference] [--subtract]
                   [--keep_unmapped]
                   [--ref_as_contig {None,trusted,untrusted}] [--no_temps]
                   [--skip_control] [-i ITERATIONS] [-v {1,2,3,4,5}]
                   [--target_len TARGET_LEN] [--DEBUG] [--DEBUG_multi]
                   [--smalt_scoring SMALT_SCORING] [-h]
                   [--spades_exe SPADES_EXE] [--samtools_exe SAMTOOLS_EXE]
                   [--smalt_exe SMALT_EXE] [--quast_exe QUAST_EXE]
                   seed_dir

Given regions from riboSnag, assembles the mapped reads

positional arguments:
  seed_dir              path to roboSnag results directory

required named arguments:
  -F FASTQ1, --fastq1 FASTQ1
                        forward fastq reads, can be compressed
  -R FASTQ2, --fastq2 FASTQ2
                        reverse fastq reads, can be compressed
  -r REFERENCE_GENOME, --reference_genome REFERENCE_GENOME
                        fasta reference, used to estimate insert sizes, and
                        compare with QUAST
  -o OUTPUT, --output OUTPUT
                        output directory; default: cwd

optional arguments:
  -S FASTQS, --fastq_single FASTQS
                        single fastq reads
  -n EXP_NAME, --experiment_name EXP_NAME
                        prefix for results files; default: riboSeed
  -m METHOD, --method_for_map METHOD
                        available mappers: smalt; default: smalt
  -c CORES, --cores CORES
                        cores for multiprocessing workers; default: 1
  -k KMERS, --kmers KMERS
                        kmers used for final assembly, separated by commas;
                        default: 21,33,55,77,99,127
  -p PRE_KMERS, --pre_kmers PRE_KMERS
                        kmers used during seeding assemblies, separated bt
                        commas; default: 21,33,55
  -g MIN_GROWTH, --min_growth MIN_GROWTH
                        skip remaining iterations if contig doesnt extend by
                        --min_growth. if 0, ignore; default: 0
  -s MIN_SCORE_SMALT, --min_score_SMALT MIN_SCORE_SMALT
                        min score forsmalt mapping; inferred from read length;
                        default: inferred
  -a MIN_ASSEMBLY_LEN, --min_assembly_len MIN_ASSEMBLY_LEN
                        if initial SPAdes assembly largest contig is not at
                        least as long as --min_assembly_len, exit. Set this to
                        the length of the seed sequence; if it is not
                        achieved, seeding across regions will likely fail;
                        default: 4000
  --paired_inference    if --paired_inference, mapped read's pairs are
                        included; default: False
  --subtract            if --subtract, reads aligned to each reference will
                        not be aligned to future iterations. Probably you
                        shouldnt do thisunless you really happen to want to
  --keep_unmapped       if --keep_unmapped, fastqs are generated containing
                        unmapped reads; default: False
  --ref_as_contig {None,trusted,untrusted}
                        if 'trusted', SPAdes will use the seed sequences as a
                        --trusted-contig; if 'untrusted', SPAdes will treat as
                        --untrusted-contig. if '', seeds will not be used
                        during assembly. See SPAdes docs; default: untrusted
  --no_temps            if --no_temps, mapping files will be removed after all
                        iterations completed; default: False
  --skip_control        if --skip_control, no SPAdes-only de novo assembly
                        will be done; default: False
  -i ITERATIONS, --iterations ITERATIONS
                        if iterations>1, multiple seedings will occur after
                        assembly of seed regions; if setting --target_len,
                        seedings will continue until --iterations are
                        completed or target_len is matched or exceeded;
                        default: 3
  -v {1,2,3,4,5}, --verbosity {1,2,3,4,5}
                        Logger writes debug to file in output dir; this sets
                        verbosity level sent to stderr. 1 = debug(), 2 =
                        info(), 3 = warning(), 4 = error() and 5 = critical();
                        default: 2
  --target_len TARGET_LEN
                        if set, iterations will continue until contigs reach
                        this length, or max iterations (set by --iterations)
                        have been completed. Set as fraction of original seed
                        length by giving a decimal between 0 and 5, or set as
                        an absolute number of base pairs by giving an integer
                        greater than 50. Not used by default
  --DEBUG               if --DEBUG, test data will be used; default: False
  --DEBUG_multi         if --DEBUG_multiprocessing, runs seeding in single
                        loop instead of a multiprocessing pool: False
  --smalt_scoring SMALT_SCORING
                        submit custom smalt scoring via smalt -S scorespec
                        option; default: match=1,subst=-4,gapopen=-4,gapext=-3
  -h, --help            Displays this help message
  --spades_exe SPADES_EXE
                        Path to SPAdes executable; default: spades.py
  --samtools_exe SAMTOOLS_EXE
                        Path to samtools executable; default: samtools
  --smalt_exe SMALT_EXE
                        Path to smalt executable; default: smalt
  --quast_exe QUAST_EXE
                        Path to quast executable; default: quast.py

```

## Known Bugs

* It looks like the `--paired-inference` option will cause an error with `SPAdes` if you are submitting a fastq of singleton/unpaired reads as part of the assembly.

* Submitting `--smalt_scoring` with vastly different scoring schemes usually causes an error.


## Installation

The trickiest part of this whole business is propperly installing SMALT.
Installing with pip3.5 will be the easiest way, but prior to release, clone this repository, and run setup.py.

### Python Requirements:

* Python v3.5 or higher
* Biopython
* pyutilsrnw

### External Requirements

* R (don't ask...)
* SPAdes v3.8 or higher
* SMALT (tested with 0.7.6)
* SAMTools (tested with 1.3.1)
* QUAST (tested with 4.1)

## Note on installation of SMALT

Must have bambamc installed! If you get an error as follows,

```
Error running test to check bambamc lib is installed!

```

it means that SMALT was installed without the bambamc library.  To install properly, see the install link below.  in short,

```
# download libtool if needed
wget ftp://ftp.gnu.org/gnu/libtool/libtool-2.4.tar.xz
tar xf libtool-2.4
cd libtool-2.4
./configure --prefix=$HOME
make
make install
export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH

# get bambamc
git clone https://github.com/gt1/bambamc.git
autoreconf -i -f
./configure --prefix=$HOME
make
make install

# download and install SMALT
wget https://sourceforge.net/projects/smalt/files/smalt-0.7.6-static.tar.gz
tar xf smalt-0.7.6-static
cd smalt-0.7.6-smalt
# this is the key step
./configure --with-bambamc=yes BAMBAMC_CFLAGS="-I$HOME/include" BAMBAMC_LIBS="-L$HOME/lib -lbambamc" --prefix=$HOME
$
make
make install

```
https://sourceforge.net/projects/smalt/files/


## Suggested Running
There are a lot of commandline options for these, so it can help to run concurrently as  script. The included `example_batch.sh` script is run as follows:

```
USAGE: /path/to/genome.gb path/to/genome.fasta path/to/read1 path/to/read2 /path/to/outdir/ n_iterations n_flanking n_cores
All mandatory arguments: genome.gb, genome.fasta, read1, read2, output_dir, iterations, flanking_width, and n_cores

example:

./example_batch.sh ./sample_data/NC_011751.1.gb ./sample_data/NC_011751.1.fasta ./sample_data/toy_set/toy_reads1.fq  ./sample_data/toy_set/toy_reads2.fq ./results_dir/ 3 1000 4


```

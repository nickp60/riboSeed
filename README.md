[![Build Status](https://travis-ci.org/nickp60/riboSeed.svg?branch=master)](https://travis-ci.org/nickp60/riboSeed)
[![Coverage Status](https://coveralls.io/repos/github/nickp60/riboSeed/badge.svg?branch=master)](https://coveralls.io/github/nickp60/riboSeed?branch=master)

![riboSeed: its whats for dinner](https://github.com/nickp60/riboSeed/blob/master/icon/logo_1.png)

# RiboSeed Pipeline
Impatient? See our [Quickstart Guide](./quickstart.md)

## Warning
This pipeline is stil very much in testing. Please back up any and all data used, and work within a virtualenv just in case.

## Description

RiboSeed is an supplemental assembly method to try to address the issue of multiple ribosomal regions in a genome.  It takes advantage of the fact that while each region is identical, the regions flanking are unique, and therefore can potentially be used to seed an assembly in such a way that rDNA regions are bridged.

The pipeline (currently) consists of optional preprocessing and two main stages:

### 0: Preprocessing with riboScan.py

This preprocesses seqeuences straight from DNA fasta if a genbank file with annotated rRNA products is not available.  The issue with many legacy annotations, assemblies, and scaffold collections is they are often poorly annotated at best, and unannotated at worst.  This is shortcut to happiness without using the full Prokka annotation scheme. It requires [`barrnap`](http://www.vicbioinformatics.com/software.barrnap.shtml) and seqret (from [`emboss`](http://www.ebi.ac.uk/Tools/emboss/))  to be availible in your path.

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

NOTE 2: In order to stremline things, as of version 0.0.3 there will be a commented header line with the feature type in the format "#$ FEATURE <featuretype>", such as "#S FEATURE rRNA".

```
#$ FEATURE rRNA
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

THIS SCRIPT IS NOT LONGER A REQUIRED PART OF THE PIPELINE! It is still included as the plots it generates can be useful for troubleshooting.

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
* `riboSeed.py` is used to map reads to the extracted regions in an iterative manner, assembling the extracted reads into long reads, and then running `SPAdes` assembly to hopefully resolve the contig junctions.  RiboSeed2 differs from the legacy version of riboSeed as riboSeed maps to all the regions at once, which minimizes the complications asscociated with multiple mappings. Instead of mapping to all the regions individually, it concatenates them into a faux genome with 10kb spacers of N's in between.  Because of this, it has the added benefit of running much faster.

#### Output

The results directory will contain a 'final_long_reads' directory with all the extended fragments, the mapped fastq files, and a `de_novo` and `de_fere_novo` folder, containing the results with the *de novo* mapping and supplemented mapping, respectively.

#### Usage:
##### NOTE:
If using a comsumer-grade computer, it may be advantagous to run with -z (--serialize).  SPAdes (or any other assembler) tends to require a lot of RAM. So unless you have about 3gb RAM per core,

minimal usage: riboSeed.py clustered\_accession\_list.txt -F FASTQ1 -R FASTQ2 -r REFERENCE_GENOME -o OUTPUT

```
usage: riboSeed.py -F FASTQ1 -R FASTQ2 -r REFERENCE_GENBANK -o OUTPUT
                    [-S FASTQS] [-n EXP_NAME] [-l FLANKING] [-m {smalt,bwa}]
                    [-c CORES] [-k KMERS] [-p PRE_KMERS] [-I] [-s SCORE_MIN]
                    [--include_shorts] [-a MIN_ASSEMBLY_LEN]
                    [--paired_inference] [--linear] [--padding PADDING]
                    [--keep_unmapped]
                    [--ref_as_contig {None,trusted,untrusted}] [--keep_temps]
                    [--skip_control] [-i ITERATIONS] [-v {1,2,3,4,5}]
                    [--target_len TARGET_LEN] [-t {1,2,4}] [-z]
                    [--smalt_scoring SMALT_SCORING] [-h]
                    [--spades_exe SPADES_EXE] [--samtools_exe SAMTOOLS_EXE]
                    [--smalt_exe SMALT_EXE] [--bwa_exe BWA_EXE]
                    [--quast_exe QUAST_EXE] [--python2_7_exe PYTHON2_7_EXE]
                    clustered_loci_txt

Given regions from riboSnag, assembles the mapped reads

positional arguments:
  clustered_loci_txt    output from riboSelect

required named arguments:
  -F FASTQ1, --fastq1 FASTQ1
                        forward fastq reads, can be compressed
  -R FASTQ2, --fastq2 FASTQ2
                        reverse fastq reads, can be compressed
  -r REFERENCE_GENBANK, --reference_genbank REFERENCE_GENBANK
                        fasta reference, used to estimate insert sizes, and
                        compare with QUAST
  -o OUTPUT, --output OUTPUT
                        output directory; default:
                        /home/nicholas/GitHub/riboSeed

optional arguments:
  -S FASTQS, --fastq_single FASTQS
                        single fastq reads
  -n EXP_NAME, --experiment_name EXP_NAME
                        prefix for results files; default: riboSeed
  -l FLANKING, --flanking_length FLANKING
                        length of flanking regions, can be colon-separated to
                        give separate upstream and downstream flanking
                        regions; default: 1000
  -m {smalt,bwa}, --method_for_map {smalt,bwa}
                        available mappers: smalt and bwa; default: bwa
  -c CORES, --cores CORES
                        cores for multiprocessing workers; default: None
  -k KMERS, --kmers KMERS
                        kmers used for final assembly, separated by commas;
                        default: 21,33,55,77,99,127
  -p PRE_KMERS, --pre_kmers PRE_KMERS
                        kmers used during seeding assemblies, separated bt
                        commas; default: 21,33,55,77,99
  -I, --ignoreS         If true, singletons from previous mappingswill be
                        ignored. try this if you seesamtools merge errors in
                        tracebacks; default: False
  -s SCORE_MIN, --score_min SCORE_MIN
                        min score for smalt mapping; inferred from read
                        length; default: inferred
  --include_shorts      if assembled contig is smaller than
                        --min_assembly_len, contig will still be included in
                        assembly; default: inferred
  -a MIN_ASSEMBLY_LEN, --min_assembly_len MIN_ASSEMBLY_LEN
                        if initial SPAdes assembly largest contig is not at
                        least as long as --min_assembly_len, exit. Set this to
                        the length of the seed sequence; if it is not
                        achieved, seeding across regions will likely fail;
                        default: 6000
  --paired_inference    if --paired_inference, mapped read's pairs are
                        included; default: False
  --linear              if genome is known to not be circular and a region of
                        interest (including flanking bits) extends past
                        chromosome end, this extends the seqence past
                        chromosome origin forward by --padding; default: False
  --padding PADDING     if treating as circular, this controls the length of
                        sequence added to the 5' and 3' ends to allow for
                        selecting regions that cross the chromosome's origin;
                        default: 5000
  --keep_unmapped       if --keep_unmapped, fastqs are generated containing
                        unmapped reads; default: False
  --ref_as_contig {None,trusted,untrusted}
                        if 'trusted', SPAdes will use the seed sequences as a
                        --trusted-contig; if 'untrusted', SPAdes will treat as
                        --untrusted-contig. if '', seeds will not be used
                        during assembly. See SPAdes docs; default: untrusted
  --keep_temps          if not --keep_temps, mapping files will be removed
                        once they are no no longer needed during the
                        iterations; default: False
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
  -t {1,2,4}, --threads {1,2,4}
                        if your cores are hyperthreaded, set number threads to
                        the number of threads per processer.If unsure, see
                        'cat /proc/cpuinfo' under 'cpu cores', or 'lscpu'
                        under 'Thread(s) per core'.: 1
  -z, --serialize       if --serialize, runs seeding in single loop instead of
                        a multiprocessing pool: False
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
  --bwa_exe BWA_EXE     Path to BWA executable; default: bwa
  --quast_exe QUAST_EXE
                        Path to quast executable; default: quast.py
  --python2_7_exe PYTHON2_7_EXE
                        Path to python2.7 executable, cause; QUAST won't run
                        on python3. default: python2.7
```

## Key Parameters

Results can be tuned by changing several of the default parameters.

* `--score_min`: With either SMALT or BWA, this can be used to set the minimum mapping score. If using BWA, the default is not to supply a minimum and to rely on the BWA default.  If submitting a `--score_min` to BWA, double check that it is appropriate.  It appears to be extremely sensitive to read length, and having a too-low threshold for minimum mapping can seriously ruin ones day.  Check out IGB or similar to view your mappings if greater than, say, 5% or the reads are mapping in subsequent iterations.  If using SMALT, the default minimum is chosen using this formula:
1.0 - (1.0 / (2.0 + *i*)), where *i* is the 0-based iteration.  This makes it progressivly more stringent with each iteration, starting with a minimum score of half the read length. Again, visualize your mappings if anything looks amiss.

* `--flanking_length`: Default is 1000.  That seems to be a good compramise between gaining unique sequence and not relying too much on the reference.

* `--kmers` and `--pre_kmers`: Adjust these as you otherwise would for a *de novo* assembly.

* `--min_assembly_len`:  For many bacteria, this is about 7000bp, as the rDNA regions for a typical operon of 16S 23S and 5S coding sequences combined usually are about that long.  If you are using non-standard rDNA regions, this should be adjusted to prevent spurious assemblies.

* `--ref_as_contig`:  This can be used to guide how SPAdes treats the long read sequences during the assembly.

* `--iterations`:  Each iteration typically increases the length of the long read by approximatly 5%.

* `--smalt_scoring`: You can adjust the SMALT scoring matrix to fine-tune the mapping stringency.




## Known Bugs

* It looks like the `--paired-inference` option will cause an error with `SPAdes` if you are submitting a fastq of singleton/unpaired reads as part of the assembly.

* Submitting `--smalt_scoring` with vastly different scoring schemes usually causes an error.

## Running Tests

The tests for the module can be found under the `tests` directory. I run them with the unittests module.  The tests assume the installation of all the reccommended tools.


## Installation
riboSeed is on Pypi, so you can install with pip or within a virtualenv (probably recommended):

`pip3.5 install riboSeed`

or

`mkvirtualenv riboSeedVE -p python3.5 -i riboSeed `

The trickiest part of this whole business is properly installing SMALT. BWA is definitly the easiest option.

You can also clone this repository, and run setup.py.

### Python Requirements:

* Python v3.5 or higher
* Biopython
* pyutilsrnw

### External Requirements

riboScan.py

* Barrnap
* Seqret

riboSelect.py

* R (don't ask...)

riboSnag.py

* PRANK or Mafft
* BLAST+ suite
* Barrnap (must be 0.7 or above)
** note that barrnap has certain Perl requirements that may not be included on your machine.  Ensure barrnap runs fine before trying riboSnag.py

riboSeed.py

* SPAdes v3.8 or higher
* BWA (tested with 0.7.12-r1039)
* (or SMALT (tested with 0.7.6), see note below)
* SAMTools (must be 1.3.1 or above)
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
### `example_batch.sh`
There are a lot of commandline options for riboSeed, so it can help to run the pipeline as script. The included `example_batch.sh` script is run as follows:


```
USAGE: /path/to/genome.gb path/to/genome.fasta path/to/read1 path/to/read2 /path/to/outdir/ n_iterations n_flanking n_cores
All mandatory arguments: genome.gb, genome.fasta, read1, read2, output_dir, iterations, flanking_width, and n_cores

example:

./example_batch.sh ./sample_data/NC_011751.1.gb ./sample_data/toy_set/toy_reads1.fq  ./sample_data/toy_set/toy_reads2.fq ./results_dir/ 3 1000 4

```
We recommend copying this file to your project directory, and customizing it as needed.


### `sge_batch.sh`
If you have access to a hpc, this script makes it easier to submit riboSeed jobs.  Just set up a python virtualenv, edit the 7 fields in this script, make any other modifications needed to fit your job, and submit with qsub.

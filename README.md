[![Build Status](https://travis-ci.org/nickp60/riboSeed.svg?branch=master)](https://travis-ci.org/nickp60/riboSeed)
[![Coverage Status](https://coveralls.io/repos/github/nickp60/riboSeed/badge.svg?branch=master)](https://coveralls.io/github/nickp60/riboSeed?branch=master)

![riboSeed: its whats for dinner](https://github.com/nickp60/riboSeed/blob/master/icon/logo_1.png)

# riboSeed Pipeline
Impatient? See our [Quickstart Guide](./quickstart.md)

A brief overview of the theory can be found [here](https://nickp60.github.io/riboSeed.html).

## Table of Contents

* [`Description`](./README.md#description)
* [`riboScan: Preprocessing`](./README.md#0-preprocessing)
* [`riboSelect: Determine rDNA Operon Structure`](./README.md#1-selection-and-extraction)
* [`riboSeed: Seeded Assembly`](./README.md#2-seeded-assembly)
* [`riboSeed Key Parameters`](./README.md#key-parameters)
* [`riboSwap: Assembly Refinement`](./README.md#3-assembly-refinement)
* [`Installation`](./README.md#installation)
* [`External Requirements`](./README.md#external-requirements)

## Before We Start
This pipeline is still very much in testing. Please back up any and all data used, and work within a virtualenv.

Genome assembly gobbles RAM. If you, like me, are working on a 4gb RAM lappy, don't run riboSeed in parallel and instead run in series by using the `--serialize` option.  That should prevent you from running out of RAM during the final SPAdes calls.

## Description

riboSeed is an supplemental assembly refinement method to try to address the issue of multiple ribosomal regions in a genome, as these create repeates unresolvable by short read sequencing.  It takes advantage of the fact that while each region is identical, the regions flanking are unique, and therefore can potentially be used to seed an assembly in such a way that rDNA regions are bridged.

The pipeline (currently) consists of optional preprocessing and two main stages:

## 0: Preprocessing

### `splitMultifasta.sh`
<!-- `splitMultifasta.sh` takes a single argument for the fasta file containing multiple contigs, and output them as individual contigs to the current working directory.  This is not sophisticated, it is essentially one awk call that is already shorter than this description.  `riboScan.py` will only work with single entry fastas. -->
<!-- #### usage -->
<!-- make a directory for the new files, and `cd` into it.  `splitMultifasta ../contigs.fa` -->

### `riboScan.py`
`riboScan.py` preprocesses sequences straight from a multifasta or  one or more fasta. The issue with many legacy annotations, assemblies, and scaffold collections is rDNAs are often poorly annotated at best, and unannotated at worst.  This is shortcut to happiness without using the full Prokka annotation scheme. It requires [`barrnap`](http://www.vicbioinformatics.com/software.barrnap.shtml) and seqret (from [`emboss`](http://www.ebi.ac.uk/Tools/emboss/))  to be available in your path.
#### Usage

riboScan can either use a directory of fastas or one (multi)fasta file.  If using a directory of fastas, provide the appropriate extension using the `-e` option. If using a (multi)fasta as input, it write out each entry to its own fasta in the `contigs` subdirectory that it makes in the output. For each of the fastas, the script renames complex headers (sketchy), scans with barrnap and captures the output gff.  It then edits the gff to add fake incrementing locus_tags, and uses the original sequence file through seqret to make a GenBank file that contains just annotated rRNA features. The last step is a concatenation which, whether or not there are multiple files, makes a single (possibly multi-entry) genbank file perfect for riboSeed-ing.

```
usage: riboScan.py -o OUTPUT [-e EXT] [-k {bac,euk,arc,mito}] [-t ID_THRESH]
                   [-b BARRNAP_EXE] [-n NAME] [-c {1,2,4,8,16}]
                   [-s SEQRET_EXE] [-v {1,2,3,4,5}] [-h]
                   contigs

Given a directory of one or more chromosomes as fasta files, this facilitates
reannotation of rDNA regions with Barrnap and outputs all sequences as a
single, annotated genbank file

positional arguments:
  contigs_dir           directory containing one or more chromosomal sequences
                        in fasta format

required named arguments:
  -o OUTPUT, --output OUTPUT
                        output directory; default:
                        /home/nicholas/GitHub/riboSeed


optional arguments:
  -k {bac,euk,arc,mito}, --kingdom {bac,euk,arc,mito}
                        whether to look for eukaryotic, archaeal, or bacterial
                        rDNA; default: bac
  -e extension          extension of the chromosomal sequences, usually
                        '.fa', '.fasta' or similar; default: .fa
  -t ID_THRESH, --id_thresh ID_THRESH
                        partial rRNA hits below this threshold will be
                        ignored. default: 0.5
  -b BARRNAP_EXE, --barrnap_exe BARRNAP_EXE
                        path to barrnap executable; default: barrnap
  -n NAME, --name NAME  name to give the contig files; default: infered from
                        file
  -s SEQRET_EXE, --seqret_exe SEQRET_EXE
                        path to seqret executable, usually installed with
                        emboss; default: seqret
  -v {1,2,3,4,5}, --verbosity {1,2,3,4,5}
                        Logger writes debug to file in output dir; this sets
                        verbosity level sent to stderr. 1 = debug(), 2 =
                        info(), 3 = warning(), 4 = error() and 5 = critical();
                        default: 2
  -h, --help            Displays this help message
```
NOTE: If using a reference with long names or containing special characters, use the --name argument to rename the contigs to something a bit more convenient and less prone to errors when piping results.

## 1: Selection and Extraction

### `riboSelect.py`
`riboSelect.py` searches the genome for rRNA annotations, clusters them into likely ribosomal groups, and outputs a colon-separated list of clustered rRNA locus tags by record id.

You will probably want to preview your file to figure out the syntax used. (ie, 16s vs 16S, rRNA vs RRNA, etc...)

If not using `riboScan.py` or if not working with a bacterial genome, you will need to change `--specific_features` appropriately to reflect the annotations in your reference (ie, for a fungal genome, use `--specific_features 5_8S:18S:28S`).

NOTE: the format of the output text file is very simple, and due to the relatively small number of such coding sequences in bacterial genomes, this can be constructed by hand if the clusters do not look appropriate. The format is `genome_sequence_id locus_tag1:locus_tag2`, where each line represents a cluster. See example below, where 14 rRNAs are clustered into 6 groups:

NOTE 2: In order to streamline things, as of version 0.0.3 there will be a commented header line with the feature type in the format "#$ FEATURE <featuretype>", such as "#S FEATURE rRNA".

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

This is used to identify and cluster rRNA regions from a gb file, returns a
text file with the clusters

positional arguments:
  genbank_genome        Genbank file (WITH SEQUENCE)

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  -o OUTPUT, --output OUTPUT
                        output directory;default:
                        /home/nicholas/GitHub/riboSeed

optional arguments:
  -f FEATURE, --feature FEATURE
                        Feature, rRNA or RRNA; default: rRNA
  -s SPECIFIC_FEATURES, --specific_features SPECIFIC_FEATURES
                        colon:separated -- specific features; default:
                        16S:23S:5S
  --keep_temps          view intermediate clustering filesdefault: False
  --clobber             overwrite previous output files: default: False
  -c CLUSTERS, --clusters CLUSTERS
                        number of rDNA clusters;if submitting multiple
                        records, must be a colon:separated list whose length
                        matches number of genbank records. Default is inferred
                        from specific feature with fewest hits
  -v VERBOSITY, --verbosity VERBOSITY
                        1 = debug(), 2 = info(), 3 = warning(), 4 = error()
                        and 5 = critical(); default: 2
  --debug               Enable debug messages



```

## Assessing rDNA regions
### `riboSnag.py`
`riboSnag.py` takes the list of clustered locus tags and extracts their sequences with flanking regions, optionally turning the coding sequences to N's to minimize bias towards reference. Is used to pull out regions of interest from a Genbank file. Outputs a directory with a fasta file for each clustered region (and a log file).

Additionally, it does a lot of plotting to visualize the Shannon entropy, coverage, occurrences, and other useful metrics.

THIS SCRIPT IS NOT LONGER A REQUIRED PART OF THE PIPELINE! It is still included as the plots it generates can be useful for troubleshooting.

#### Usage:

```
usage: riboSnag.py [-o OUTPUT] [-n NAME] [-l FLANKING] [--msa_kmers] [-c]
                   [-p PADDING] [-v VERBOSITY] [--clobber] [--no_revcomp]
                   [--skip_check] [--msa_tool {mafft,prank}]
                   [--prank_exe PRANK_EXE] [--mafft_exe MAFFT_EXE]
                   [--barrnap_exe BARRNAP_EXE]
                   [--makeblastdb_exe MAKEBLASTDB_EXE]
                   [--kingdom {mito,euk,arc,bac}] [-h]
                   genbank_genome clustered_loci

Use to extract regions of interest based on supplied Locus tags and evaluate
the extracted regions

positional arguments:
  genbank_genome        Genbank file (WITH SEQUENCE)
  clustered_loci        output from riboSelect

required named arguments:
  -o OUTPUT, --output OUTPUT
                        output directory; default:
                        /home/nicholas/GitHub/riboSeed

optional arguments:
  -n NAME, --name NAME  rename the contigs with this prefixdefault: date
                        (YYYYMMDD)
  -l FLANKING, --flanking_length FLANKING
                        length of flanking regions, in bp; default: 1000
  --msa_kmers           calculate kmer similarity based on aligned sequences
                        instead of raw sequences;default: False
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
  --no_revcomp          default returns reverse complimented seq if majority
                        of regions on reverse strand. if --no_revcomp, this is
                        overwriddendefault: False
  --skip_check          Dont bother calculating Shannon Entropy; default:
                        False
  --msa_tool {mafft,prank}
                        Path to PRANK executable; default: mafft
  --prank_exe PRANK_EXE
                        Path to PRANK executable; default: prank
  --mafft_exe MAFFT_EXE
                        Path to MAFFT executable; default: mafft
  --barrnap_exe BARRNAP_EXE
                        Path to barrnap executable; default: barrnap
  --makeblastdb_exe MAKEBLASTDB_EXE
                        Path to makeblastdb executable; default: makeblastdb
  --kingdom {mito,euk,arc,bac}
                        kingdom for barrnap; default: bac
  -h, --help            Displays this help message
```

## 2: Seeded Assembly
### `riboSeed.py`
`riboSeed.py` maps reads to a genome and (1) extracts reads mapping to rDNA regions, (2) perfoms subassemblies on each pool of extracted reads to recover the rDNA complete with flanking regions (resulting in a pseudocontig) (3) concatenates a;; pseudocontigs into them into a pseudogenome with 5kb spacers of N's in between, (5) map remaining reads to the pseudogenome, and (6) repeat steps 1-5 for a given number of iterations (default 3 iterations). Finally, riboSeed runs SPAdes assemblied with and without the pseudocontigs and the resulting assemblies are assessed with QUAST.

#### Output

The results directory will contain a 'final_long_reads' directory with all the pseudocontigs, the mapped fastq files, and `final_de_novo_assembly` and `final_de_fere_novo_assembly` folders, containing the SPAdes results.

##### NOTE:
If using a consumer-grade computer, it will be advantagous to run with `-z/--serialize` enabled to run asseblies in serial rather than parallel.

#### Usage:
minimal usage:
```riboSeed.py clustered_accession\list.txt -F FASTQ1 -R FASTQ2 -r REFERENCE_GENOME -o OUTPUT```

```
usage: riboSeed.py -r REFERENCE_GENBANK -o OUTPUT [-F FASTQ1] [-R FASTQ2]
                   [-S1 FASTQS1] [-n EXP_NAME] [-l FLANKING] [-m {smalt,bwa}]
                   [-c CORES] [-k KMERS] [-p PRE_KMERS] [-s SCORE_MIN]
                   [-a MIN_ASSEMBLY_LEN] [--include_shorts] [--linear]
                   [--ref_as_contig {None,trusted,untrusted}] [--keep_temps]
                   [--skip_control] [-i ITERATIONS] [-v {1,2,3,4,5}]
                   [--target_len TARGET_LEN] [-t {1,2,4}] [-z]
                   [--smalt_scoring SMALT_SCORING] [--mapper_args MAPPER_ARGS]
                   [-h] [--spades_exe SPADES_EXE]
                   [--samtools_exe SAMTOOLS_EXE] [--smalt_exe SMALT_EXE]
                   [--bwa_exe BWA_EXE] [--quast_exe QUAST_EXE]
                   [--python2_7_exe PYTHON2_7_EXE]
                   clustered_loci_txt

Given cluster file of rDNA regions from riboSelect and either paired-end or
single-end reads, assembles the mapped reads into pseduocontig 'seeds', and
uses those with the reads to runde fere novo and de novo assembly with SPAdes

positional arguments:
  clustered_loci_txt    output from riboSelect

required named arguments:
  -r REFERENCE_GENBANK, --reference_genbank REFERENCE_GENBANK
                        genbank reference, used to estimate insert sizes, and
                        compare with QUAST
  -o OUTPUT, --output OUTPUT
                        output directory; default:
                        /home/nicholas/GitHub/riboSeed

optional arguments:
  -F FASTQ1, --fastq1 FASTQ1
                        forward fastq reads, can be compressed
  -R FASTQ2, --fastq2 FASTQ2
                        reverse fastq reads, can be compressed
  -S1 FASTQS1, --fastq_single1 FASTQS1
                        single fastq reads
  -n EXP_NAME, --experiment_name EXP_NAME
                        prefix for results files; default: riboSeed
  -l FLANKING, --flanking_length FLANKING
                        length of flanking regions, in bp; default: 1000
  -m {smalt,bwa}, --method_for_map {smalt,bwa}
                        available mappers: smalt and bwa; default: bwa
  -c CORES, --cores CORES
                        cores for multiprocessing; default: None
  -k KMERS, --kmers KMERS
                        kmers used for final assembly, separated by commas;
                        default: 21,33,55,77,99,127
  -p PRE_KMERS, --pre_kmers PRE_KMERS
                        kmers used during seeding assemblies, separated bt
                        commas; default: 21,33,55,77,99
 -s SCORE_MIN, --score_min SCORE_MIN
                        If using smalt, this sets the '-m' param; default with
                        smalt is inferred from read length. If using BWA,
                        reads mapping with ASscore lower than this will be
                        rejected; default with SWA is half of read length
   -a MIN_ASSEMBLY_LEN, --min_assembly_len MIN_ASSEMBLY_LEN
                        if initial SPAdes assembly largest contig is not at
                        least as long as --min_assembly_len, exit. Set this to
                        the length of the seed sequence; if it is not
                        achieved, seeding across regions will likely fail;
                        default: 6000
  --include_shorts      if assembled contig is smaller than
                        --min_assembly_len, contig will still be included in
                        assembly; default: inferred
  --linear              if genome is known to not be circular and a region of
                        interest (including flanking bits) extends past
                        chromosome end, this extends the seqence past
                        chromosome origin forward by --padding; default: False
  --ref_as_contig {None,trusted,untrusted}
                        if 'trusted', SPAdes will use the seed sequences as a
                        --trusted-contig; if 'untrusted', SPAdes will treat as
                        --untrusted-contig. if '', seeds will not be used
                        during assembly. See SPAdes docs; default: untrusted
  --keep_temps          if not --keep_temps, mapping files will be removed
                        once they are no no longer needed during the
                        iterations; default: False
  --skip_control        if --skip_control, no de novo assembly will be done;
                        default: False
  -i ITERATIONS, --iterations ITERATIONS
                        if iterations>1, multiple seedings will occur after
                        subassembly of seed regions; if setting --target_len,
                        seedings will continue until --iterations are
                        completed or --target_len is matched or exceeded;
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
  -z, --serialize       if --serialize, runs seeding and assembly without
                        multiprocessing. This is recommended for machines with
                        less than 8GB RAM: False
  --smalt_scoring SMALT_SCORING
                        if mapping with SMALT, submit custom smalt scoring via
                        smalt -S scorespec option; default:
                        match=1,subst=-4,gapopen=-4,gapext=-3
  --mapper_args MAPPER_ARGS
                        submit custom parameters to mapper. And by mapper, I
                        mean bwa, cause we dont support this option for SMALT,
                        sorry. This requires knowledge of your chosen mapper's
                        optional arguments. Proceed with caution! default: -L
                        0,0 -U 0
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
                        Path to python2.7 executable, cause QUAST won't run on
                        python3. default: python2.7
```

## Key Parameters

Results can be tuned by changing several of the default parameters.

* `--score_min`: With either SMALT or BWA, this can be used to set the minimum mapping score. If using BWA, the default is not to supply a minimum and to rely on the BWA default.  If submitting a `--score_min` to BWA, double check that it is appropriate.  It appears to be extremely sensitive to read length, and having a too-low threshold for minimum mapping can seriously ruin ones day.  Check out IGB or similar to view your mappings if greater than, say, 5% or the reads are mapping in subsequent iterations.  If using SMALT, the default minimum is chosen using this formula:
1.0 - (1.0 / (2.0 + *i*)), where *i* is the 0-based iteration.  This makes it progressively more stringent with each iteration, starting with a minimum score of half the read length. Again, visualize your mappings if anything looks amiss.

* `--flanking_length`: Default is 2000.  That seems to be a good compromise between gaining unique sequence and not relying too much on the reference.

* `--kmers` and `--pre_kmers`: Adjust these as you otherwise would for a *de novo* assembly.

* `--min_assembly_len`:  For bacteria, this is about 7000bp, as the rDNA regions for a typical operon of 16S 23S and 5S coding sequences combined usually are about that long.  If you are using non-standard rDNA regions, this should be adjusted to prevent spurious assemblies.

* `--ref_as_contig`:  This can be used to guide how SPAdes treats the long read sequences during the assembly (`trusted` or `untrusted`).

* `--iterations`:  Each iteration typically increases the length of the long read by approximately 5%.

* `--smalt_scoring`: You can adjust the SMALT scoring matrix to fine-tune the mapping stringency.

## 3: Assembly Refinement
### `riboSwap.py`
Infrequently, `riboSeed` has joined together contigs that appear incorrect according to your reference.  If you are at all unhappy with a bridging, `riboSwap.py` allows swapping of a "bad" contig for one or more syntenic contigs from the *de novo* assembly.
#### USAGE
```
usage: riboSwap.py -o OUTPUT [-v {1,2,3,4,5}] [-h]
                   de_novo_file de_fere_novo_file bad_contig good_contigs

Given de novo and de fere novo contigs files, a misjoined de fere novo contig
name, and a colon:separated list of de novo contig names, replace the
offending contig with the de novo contig(s)

positional arguments:
  de_novo_file          multifasta containing de novo contigs
  de_fere_novo_file     multifasta containing de fere novo contigs
  bad_contig            name of the bad contig
  good_contigs          colon separated good contigs for replacement

required named arguments:
  -o OUTPUT, --output OUTPUT
                        output directory; default:
                        /home/nicholas/GitHub/riboSeed

optional arguments:
  -v {1,2,3,4,5}, --verbosity {1,2,3,4,5}
                        Logger writes debug to file in output dir; this sets
                        verbosity level sent to stderr. 1 = debug(), 2 =
                        info(), 3 = warning(), 4 = error() and 5 = critical();
                        default: 2
  -h, --help            Displays this help message
```

## Known Bugs

* Submitting `--smalt_scoring` with vastly different scoring schemes usually causes an error.

## Running Tests

The tests for the module can be found under the `tests` directory. I run them with the unittests module.  The tests assume the installation of all the recommended tools.


## Installation
#### From Pypi
riboSeed is on Pypi, so you can install with pip, preferably within a virtualenv (recommended):

`virtualenv -p python3.5 venv-riboSeed`
`source venv-riboSeed/bin/activate`
`pip3.5 install riboSeed`

#### From GitHub
You can also clone this repository, and run `python3.5 setup.py install`.

In the scripts directory, there is a script called `runme.py` which run the pipeline on a small sample dataset.  it should output a folder to your current working dir called `integration_tests`.


### Dependencies
The trickiest part of this whole business is properly installing SMALT. BWA is definitely the easiest option, and the current default mapper, so don't bother with SMALT unless you need to.

### Python Requirements:

* Python >= v3.5
* Biopython v1.68
* pysam v0.9.1.4,
* pyutilsnrw >= 0.0.768
* matplotlib v1.5.3
* pandas v0.18.1

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
USAGE: /path/to/genome.gb  path/to/read1 path/to/read2 /path/to/outdir/ n_iterations n_flanking n_cores
All mandatory arguments: genome.gb, genome.fasta, read1, read2, output_dir, iterations, flanking_width, and n_cores

example:

./example_batch.sh ./sample_data/NC_011751.1.gb ./sample_data/toy_set/toy_reads1.fq  ./sample_data/toy_set/toy_reads2.fq ./results_dir/ 3 1000 4

```
We recommend copying this file to your project directory, and customizing it as needed.


### `sge_batch.sh`
If you have access to a hpc, this script makes it easier to submit riboSeed jobs.  Just set up a python virtualenv, edit the 7 fields in this script, make any other modifications needed to fit your job, and submit with `qsub`.

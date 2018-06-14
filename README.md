[![Build Status](https://travis-ci.org/nickp60/riboSeed.svg?branch=master)](https://travis-ci.org/nickp60/riboSeed)
[![PyPI version](https://badge.fury.io/py/riboSeed.svg)](https://badge.fury.io/py/riboSeed)
[![Coverage Status](https://coveralls.io/repos/github/nickp60/riboSeed/badge.svg?branch=master)](https://coveralls.io/github/nickp60/riboSeed?branch=master)
[![DOI](https://zenodo.org/badge/68617544.svg)](https://zenodo.org/badge/latestdoi/68617544)
[![Code Health](https://landscape.io/github/nickp60/riboSeed/master/landscape.svg?style=flat)](https://landscape.io/github/nickp60/riboSeed/master)
[![Documentation Status](https://readthedocs.org/projects/riboseed/badge/?version=latest)](http://riboseed.readthedocs.io/en/latest/?badge=latest)

![riboSeed](https://github.com/nickp60/riboSeed/blob/master/icon/logo_1.png)

# riboSeed Pipeline
Impatient? See our [Quickstart Guide](./quickstart.md)

A brief overview of the theory can be found [here](https://nickp60.github.io/riboSeed.html).

The riboSeed manuscript can be found [here](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky212/4955760).

## Citation
```
Nicholas R Waters, Florence Abram, Fiona Brennan, Ashleigh Holmes, Leighton Pritchard;
riboSeed: leveraging prokaryotic genomic architecture to assemble across ribosomal regions,
Nucleic Acids Research, gky212, https://doi.org/10.1093/nar/gky212
```

Interested in the figures/tables/analyses in the manuscript?  See the [README](https://github.com/nickp60/riboSeed/blob/master/scripts/README.md) in the `scripts` dir.

## Table of Contents

* [`Description`](./README.md#description)
* [`scan: Preprocessing`](./README.md#0-preprocessing)
* [`select: Determine rDNA Operon Structure`](./README.md#1-selection-and-extraction)
* [`seed: Seeded Assembly`](./README.md#2-seeded-assembly)
* [`riboSeed Key Parameters`](./README.md#key-parameters)
* [`swap: Assembly Refinement`](./README.md#3-assembly-refinement)
* [`Installation`](./README.md#installation)
* [`External Requirements`](./README.md#external-requirements)

## Before We Start
Please back up any and all data used, and work within a virtualenv.

Genome assembly gobbles RAM. If you, like me, are working on a 4gb RAM lappy, don't run riboSeed in parallel and instead run in series by using the `--serialize` option.  That should prevent you from running out of RAM during the final SPAdes calls.

## Description

riboSeed is an supplemental assembly refinement method to try to address the issue of multiple ribosomal regions in a genome, as these create repeats unresolvable by short read sequencing.  It takes advantage of the fact that while each region is identical, the regions flanking are unique, and therefore can potentially be used to seed an assembly in such a way that rDNA regions are bridged.

Preprocessing
- `scan`
- `select`

*De fere novo assembly*
- `seed`

Visualizations/assessment
- `snag`
- `stack`
- `sketch`
- `swap`
- `score`

## 0: Preprocessing

### `scan`
`scan` preprocesses sequences straight from a multifasta or  one or more fasta. The issue with many legacy annotations, assemblies, and scaffold collections is rDNAs are often poorly annotated at best, and unannotated at worst.  This is shortcut to happiness without using the full Prokka annotation scheme. It requires [`barrnap`](http://www.vicbioinformatics.com/software.barrnap.shtml) and seqret (from [`emboss`](http://www.ebi.ac.uk/Tools/emboss/))  to be available in your path.
#### Usage

riboScan can either use a directory of fastas or one (multi)fasta file.  If using a directory of fastas, provide the appropriate extension using the `-e` option. If using a (multi)fasta as input, it write out each entry to its own fasta in the `contigs` subdirectory that it makes in the output. For each of the fastas, the script renames complex headers (sketchy), scans with barrnap and captures the output gff.  It then edits the gff to add fake incrementing locus_tags, and uses the original sequence file through seqret to make a GenBank file that contains just annotated rRNA features. The last step is a concatenation which, whether or not there are multiple files, makes a single (possibly multi-entry) genbank file perfect for riboSeed-ing.

NOTE: If using a reference with long names or containing special characters, use the --name argument to rename the contigs to something a bit more convenient and less prone to errors when piping results.


### `select`
`select` searches the genome for rRNA annotations, clusters them into likely ribosomal groups, and outputs a colon-separated list of clustered rRNA locus tags by record id.

You will probably want to preview your file to figure out the syntax used. (ie, 16s vs 16S, rRNA vs RRNA, etc...)

If not using `scan` or if not working with a prokaryotic genome, you will need to change `--specific_features` appropriately to reflect the annotations in your reference (ie, for a fungal genome, use `--specific_features 5_8S:18S:28S`).

NOTE: the format of the output text file is very simple, and due to the relatively small number of such coding sequences in bacterial genomes, this can be constructed by hand if the clusters do not look appropriate. The format is `genome_sequence_id locus_tag1:locus_tag2`, where each line represents a cluster. See example below, where 14 rRNAs are clustered into 6 groups:

NOTE 2: In order to streamline things, as of version 0.0.3 there will be a commented header line with the feature type in the format "#$ FEATURE <featuretype>", such as `#$ FEATURE rRNA`.

```
#$ FEATURE rRNA
CM000577.1 FGSG_20052:FGSG_20051:FGSG_20053
CM000577.1 FGSG_20048:FGSG_20047
CM000577.1 FGSG_20049:FGSG_20050
CM000577.1 FGSG_20054:FGSG_20056:FGSG_20055
CM000577.1 FGSG_20058:FGSG_20057
CM000577.1 FGSG_20075:FGSG_20074
```

## 2: *De fere novo* Assembly
### `seed`
`seed` maps reads to a genome and (1) extracts reads mapping to rDNA regions, (2) performs subassemblies on each pool of extracted reads to recover the rDNA complete with flanking regions (resulting in a pseudocontig) (3) concatenates a;; pseudocontigs into them into a pseudogenome with 5kb spacers of N's in between, (5) map remaining reads to the pseudogenome, and (6) repeat steps 1-5 for a given number of iterations (default 3 iterations). Finally, riboSeed runs SPAdes assemblies with and without the pseudocontigs and the resulting assemblies are assessed with QUAST.

#### Output

The results directory will contain a 'final_long_reads' directory with all the pseudocontigs, the mapped fastq files, and `final_de_novo_assembly` and `final_de_fere_novo_assembly` folders, containing the SPAdes results.

##### NOTE:
If using a consumer-grade computer, it will be advantageous to run with `-z/--serialize` enabled to run assemblies in serial rather than parallel.


## Key Parameters

Results can be tuned by changing several of the default parameters.

* `--score_min`:  This can be used to set the minimum mapping score. If using BWA, the default is not to supply a minimum and to rely on the BWA default.  If submitting a `--score_min` to BWA, double check that it is appropriate.  It appears to be extremely sensitive to read length, and having a too-low threshold for minimum mapping can seriously ruin ones day.  Check out IGB or similar to view your mappings if greater than, say, 5% or the reads are mapping in subsequent iterations.

* `-l, --flanking_length`: Default is 2000.  That seems to be a good compromise between gaining unique sequence and not relying too much on the reference.

* `--kmers` and `--pre_kmers`: Adjust these as you otherwise would for a *de novo* assembly.

* `--min_assembly_len`:  For bacteria, this is about 7000bp, as the rDNA regions for a typical operon of 16S 23S and 5S coding sequences combined usually are about that long.  If you are using non-standard rDNA regions, this should be adjusted to prevent spurious assemblies.

* `--ref_as_contig`:  This can be used to guide how SPAdes treats the long read sequences during the assembly (`trusted` or `untrusted`).  By default, this is inferred from mapping percentage (`trusted` if over 85% of reads map to the reference)

* `--iterations`:  Each iteration typically increases the length of the long read by approximately 5%.

## 3: Visualization/Assessment

### `snag`
`snag` takes the list of clustered locus tags and extracts their sequences with flanking regions, optionally turning the coding sequences to N's to minimize bias towards reference. Is used to pull out regions of interest from a Genbank file. Outputs a directory with a fasta file for each clustered region (and a log file).

Additionally, it does a lot of plotting to visualize the Shannon entropy, coverage, occurrences, and other useful metrics.


#### Usage:

```
usage: ribo snag  [-o OUTPUT] [-n NAME] [-l FLANKING] [--msa_kmers] [-c]
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
  -n NAME, --name NAME  rename the contigs with this prefix; default: date
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
  --clobber             overwrite previous output files default: False
  --no_revcomp          default returns reverse complimented seq if majority
                        of regions on reverse strand. if --no_revcomp, this is
                        overwridden; default: False
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

### `stack`
Because assembly using short reads often collapses rDNA repeats, it is not uncommon to find a reference genome that has less than the actual number of rDNAs.  riboStack uses `bedtools` and `samtools` to determine the coverage across rDNA regions, and compares that coverage depth to 10 sets of randomly selected non-rDNA regions.  If the number of rDNAs in the reference matches the number of rDNAs in your sequenced isolate, the coverage should be pretty similar. However, if the coverage in your rDNA regions is significantly higher, than there are likely more rDNAs in your sequenced isolate that there are in the reference, which is something to be aware of.

It requires a mapping BAM file and the riboScan output directory as input.



### `swap`
Infrequently, `seed` has joined together contigs that appear incorrect according to your reference.  If you are at all unhappy with a bridging, `swap` allows swapping of a "bad" contig for one or more syntenic contigs from the *de novo* assembly.
#### USAGE
```
usage: ribo swap -o OUTPUT [-v {1,2,3,4,5}] [-h]
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

### X server

You may run into issues where you get an error about "Unable to connect to X server: None" or localhost:N. Sorry about that; any tips would be useful;  a quick glance at the commit history will show I have spent much time trying to resolve it, without any luck.  If you do run into this, try the following:
  - connect to the machine with an X session (`ssh -X hostname`)
  - avoid using `gnu screen` if possible, but if you do need to use it, start the `screen` session after ensuring you have a `$DISPLAY` availible through starting the host session with `-X`

### Pysam on MacOS
If you are on MacOS, you may run into an issue with Pysam.
```
ImportError: dlopen(/Users/nicholas/miniconda3/envs/ribo/lib/python3.5/site-packages/pysam/libchtslib.cpython-35m-darwin.so, 2): Library not loaded: @rpath/liblzma.5.dylib
  Referenced from: /Users/nicholas/miniconda3/envs/ribo/lib/libhts.2.dylib
  Reason: Incompatible library version: libhts.2.dylib requires version 8.0.0 or later, but liblzma.5.dylib provides version 6.0.0
```
In this case, try installing by first making a conda env with the `environment.yaml` file, and then installing riboSeed from pip.
```
conda env create -y environment.yaml
source activate ribo
pip install riboSeed
```

If you run into malloc issues similar to https://github.com/ablab/spades/issues/9, we recommend running in a VM.

### smalt scoring
Submitting `--smalt_scoring` with vastly different scoring schemes usually causes an error.

## Running Tests

The tests for the module can be found under the `tests` directory. I run them with the unittests module.  The tests assume the installation of all the recommended tools.


## Installation

### From conda (new and recommended!)
Conda is a cross-platform, cross-language package management system.  If you haven't already installed conda, follow [these instructions here](https://bioconda.github.io/index.html), and install the python3.6 version.  Once you have that done, add the appropriate channels.

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

and then install riboSeed and all of its dependencies with one command:

```
conda install riboseed
```

(Note the lowercase "s")


#### From Pypi
riboSeed is on Pypi, so you can install with pip, preferably within a virtualenv (recommended):

```
virtualenv -p python3.5 venv-riboSeed
source venv-riboSeed/bin/activate
pip3.5 install riboSeed
```


#### From TestPypi
To install the bleeding-edge version, install from testpypi:
```
virtualenv -p python3.5 venv-riboSeed
source venv-riboSeed/bin/activate
pip install --extra-index-url https://testpypi.python.org/pypi riboSeed
```

#### From GitHub
You can also clone this repository, and run `python3.5 setup.py install`.

In the scripts directory, there is a script called `runme.py` which run the pipeline on a small sample dataset.  it should output a folder to your current working dir called `integration_tests`.


### Dependencies
Python requirements can be found in the `requirements.txt` file.

### External Requirements

riboScan.py

* Barrnap (must be 0.7 or above)
* Seqret

riboSelect.py

* None

riboSnag.py

* PRANK or Mafft
* BLAST+ suite
* Barrnap (must be 0.7 or above)

riboSeed.py

* SPAdes v3.8 or higher
* BWA (tested with 0.7.12-r1039)
* SAMTools (must be 1.3.1 or above)
* QUAST (tested with 4.5)

NOTE: barrnap has certain Perl requirements that may not be included on your machine. Ensure barrnap runs fine before trying `snag.py`.  Or try [python barrnap](https://github.com/nickp60/barrnap/).


## Suggested Running

The `ribo run` command orchestrates the most commonly used sequence of calls to `scan`, `select`, `seed`, `sketch`, `score`, and so on.

```
usage: ribo run [-r reference.fasta] -c config_file [-o /output/dir/]
                [-n experiment_name] [-K {bac,euk,arc,mito}] [-S 16S:23S:5S]
                [--clusters str] [-C str] [-F reads_F.fq] [-R reads_R.fq]
                [-S1 reads_S.fq] [-s int]
                [--ref_as_contig {ignore,infer,trusted,untrusted}] [--linear]
                [-j] [--score] [-l int] [-k 21,33,55,77,99,127]
                [--force_kmers] [-p 21,33,55,77,99] [-d int] [--clean_temps]
                [-i int] [-v {1,2,3,4,5}] [--cores int] [--memory int]
                [--damn_the_torpedos] [-t {1,2,4}] [-z] [-h] [--version]

Run the riboSeed pipeline of scan, select, seed, sketch, and score. Uses a
config file to wrangle all the args not available via these commandline args.
This can either be run by providing (as minimum) a reference, some reads, and
an output directory; or, if you have a completed config file, you can run it
with just that.

optional arguments:
  -r reference.fasta, --reference_fasta reference.fasta
                        path to a (multi)fasta or a directory containing one
                        or more chromosomal sequences in fasta format.
                        Required, unless using a config file
  -c config_file, --config config_file
                        config file; if none given, create one; default:
                        /home/nicholas/GitHub/riboSeed
  -o /output/dir/, --output /output/dir/
                        output directory; default: /home/nicholas/GitHub/riboS
                        eed/2018-06-14T1353_riboSeed_pipeline_results/
  -n experiment_name, --experiment_name experiment_name
                        prefix for results files; default: inferred
  -K {bac,euk,arc,mito}, --Kingdom {bac,euk,arc,mito}
                        whether to look for eukaryotic, archaeal, or bacterial
                        rDNA; default: bac
  -S 16S:23S:5S, --specific_features 16S:23S:5S
                        colon:separated -- specific features; default:
                        16S:23S:5S
  --clusters str        number of rDNA clusters;if submitting multiple
                        records, must be a colon:separated list whose length
                        matches number of genbank records. Default is inferred
                        from specific feature with fewest hits
  -C str, --cluster_file str
                        clustered_loci file output from riboSelect;this is
                        created by default from run_riboSeed, but if you don't
                        agree with the operon structure predicted by
                        riboSelect, you can use your alternate clustered_loci
                        file. default: None
  -F reads_F.fq, --fastq1 reads_F.fq
                        path to forward fastq file, can be compressed
  -R reads_R.fq, --fastq2 reads_R.fq
                        path to reverse fastq file, can be compressed
  -S1 reads_S.fq, --fastq_single1 reads_S.fq
                        path to single fastq file
  -s int, --score_min int
                        If using smalt, this sets the '-m' param; default with
                        smalt is inferred from read length. If using BWA,
                        reads mapping with ASscore lower than this will be
                        rejected; default with BWA is half of read length
  --ref_as_contig {ignore,infer,trusted,untrusted}
                        ignore: reference will not be used in subassembly.
                        trusted: SPAdes will use the seed sequences as a
                        --trusted-contig; untrusted: SPAdes will treat as
                        --untrusted-contig. infer: if mapping percentage over
                        80%, 'trusted'; else 'untrusted'. See SPAdes docs for
                        details. default: infer
  --linear              if genome is known to not be circular and a region of
                        interest (including flanking bits) extends past
                        chromosome end, this extends the seqence past
                        chromosome origin forward by --padding; default: False
  -j, --just_seed       Don't do an assembly, just generate the long read
                        'seeds'; default: False
  --score               run riboScore too! default: False
  -l int, --flanking_length int
                        length of flanking regions, in bp; default: 1000
  -k 21,33,55,77,99,127, --kmers 21,33,55,77,99,127
                        kmers used for final assembly, separated by commas
                        such as21,33,55,77,99,127. Can be set to 'auto', where
                        SPAdes chooses. We ensure kmers are not too big or too
                        close to read length; default: 21,33,55,77,99,127
  --force_kmers         skip checking to see if kmerchoice is appropriate to
                        read length. Sometimes kmers longer than reads can
                        help in the final assembly, as the long reads
                        generated by riboSeed contain kmers longer than the
                        read length
  -p 21,33,55,77,99, --pre_kmers 21,33,55,77,99
                        kmers used during seeding assemblies, separated bt
                        commas; default: 21,33,55,77,99
  -d int, --min_flank_depth int
                        a subassembly won't be performed if this minimum depth
                        is not achieved on both the 3' and5' end of the
                        pseudocontig. default: 0
  --clean_temps         if --clean_temps, mapping files will be removed once
                        they are no no longer needed during the mapping
                        iterations to save space; default: False
  -i int, --iterations int
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
  --cores int           cores used; default: None
  --memory int          cores for multiprocessing; default: 8
  --damn_the_torpedos   Ignore certain errors, full speed ahead!
  -t {1,2,4}, --threads {1,2,4}
                        if your cores are hyperthreaded, set number threads to
                        the number of threads per processer.If unsure, see
                        'cat /proc/cpuinfo' under 'cpu cores', or 'lscpu'
                        under 'Thread(s) per core'.: 1
  -z, --serialize       if --serialize, runs seeding and assembly without
                        multiprocessing. We recommend this for machines with
                        less than 8GB RAM: False
  -h, --help            Displays this help message
  --version             show program's version number and exit
```
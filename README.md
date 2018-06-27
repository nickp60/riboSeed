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

* [`Reference Selection`](./README.md#reference-selection)
* [`Description`](./README.md#description)
* [`Installation`](./README.md#installation)
* [`Suggested Running`](./README.md#suggested-running)
* [`Contributing`](./README.md#contributing)

## Reference Selection
`riboSeed` requires an appropriate reference genome for the *de fere novo* assembly.  We outline a few different ways to select the best refference:

- [Protocol for using Kraken.](http://riboseed.readthedocs.io/en/latest/REFERENCE.html#method-1-kraken) This kmer-based method requires having a Kraken database handy, but it provides good resolution while also checking for contamination.
- [Protocol for using Reads2Type.](http://riboseed.readthedocs.io/en/latest/REFERENCE.html#method-2-reads2type-and-cgfind) This method is quick and easy if you are less comfortable with commandline tools, but at the cost of resolution.
- [Using ANI (in development).](./choose_reference_with_ANI.md) This method should provide a appropriate reference based on ANI.  It can be time-consuming depending on the number of references considered, but can identify the best reference in a hands-free way, perfect for automation.


## Before We Start
Please back up any and all data used, and work within a virtualenv.

Genome assembly gobbles RAM. If you, like me, are working on a 4gb RAM lappy, don't run riboSeed in parallel and instead run in series by using the `--serialize` option.  That should prevent you from running out of RAM during the final SPAdes calls.

## Description

riboSeed is an supplemental assembly refinement method to try to address the issue of multiple ribosomal regions in a genome, as these create repeats unresolvable by short read sequencing.  It takes advantage of the fact that while each region is identical, the regions flanking are unique, and therefore can potentially be used to seed an assembly in such a way that rDNA regions are bridged.

For a description of each submodule, follow the links below to the readthedocs manual page.

Preprocessing
- [`scan` | annotate reference genome rRNAs](http://riboseed.readthedocs.io/en/latest/PIPELINE.html#scan)
- [`select` | identify rDNA operons](http://riboseed.readthedocs.io/en/latest/PIPELINE.html#select)

*De fere novo assembly*
- [`seed` | perform interative subassembly](http://riboseed.readthedocs.io/en/latest/PIPELINE.html#seed)

Visualizations/assessment
- [`snag` | extract and visualize rDNA regions](http://riboseed.readthedocs.io/en/latest/PIPELINE.html#snag)
- [`stack` | calculate coverage at rDNAs in final assembly](http://riboseed.readthedocs.io/en/latest/PIPELINE.html#stack)
- [`sketch` | plot the relative rDNA regions in a handful of genomes](http://riboseed.readthedocs.io/en/latest/PIPELINE.html#sketch)
- [`swap` | switch questionable contigs ](http://riboseed.readthedocs.io/en/latest/PIPELINE.html#swap)
- [`score` | automated scoring for rDNA assemblies](http://riboseed.readthedocs.io/en/latest/PIPELINE.html#score)
- [`spec` | speculate the nunber of rDNA operons based on assembly graph](http://riboseed.readthedocs.io/en/latest/PIPELINE.html#spec)


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
External requirements can be found in the `environment.yml`, and can be used to create a conda environment: (`conda env create -f environment.yml`)

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

## Contributing
Pull requests are more than welcome!

### Known Bugs

#### X server

You may run into issues where you get an error about "Unable to connect to X server: None" or localhost:N. Sorry about that; any tips would be useful;  a quick glance at the commit history will show I have spent much time trying to resolve it, without any luck.  If you do run into this, try the following:
  - connect to the machine with an X session (`ssh -X hostname`)
  - avoid using `gnu screen` if possible, but if you do need to use it, start the `screen` session after ensuring you have a `$DISPLAY` availible through starting the host session with `-X`

#### Pysam on MacOS
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

#### smalt scoring
Submitting `--smalt_scoring` with vastly different scoring schemes usually causes an error.

### Running Tests

The tests for the module can be found under the `tests` directory. I run them with the unittests module.  The tests assume the installation of all the recommended tools.

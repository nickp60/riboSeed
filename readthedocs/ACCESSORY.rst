Accessory Scripts
===============

Assessment
----

``riboScore.py``
~~~
Suppose you have a whole bunch of assemblies to assess. The most rigorous way of checking the assemblies would be to use Mauve (or a similar whole genome alignment visualizer tool) for the job, and manually check the quality of each assembly, listening to the ends of the contigs, seeking one-ness with the data.  Thats all well and good if you are (a) independantly wealthy and enjoy doing this sort of thing, (b) seeking a meditative state through mindless clicking, or (c) an undergrad assistant, but for the rest of us, we are willing to sacrifice a bit of accuracy for throughput.  This is, after all, why we aren't sequencing on gels anymore.

``riboScore.py`` outputs two types of score repors as text files: one which is
easy for humans to read, and the other that can be easily combined with
hundreds like it to make various types of graphs etc.

::
   usage: riboScore.py [-h] [-o OUTPUT] [-l FLANKING] [-p MIN_PERCENT]
                    [-f ASSEMBLY_EXT] [-g REF_EXT] [-F] [-v {1,2,3,4,5}]
                    indir

   This does some simple blasting to detect correctness of riboSeed results

   positional arguments:
     indir                 dir containing a genbank file and other file

   optional arguments:
     -h, --help            show this help message and exit
     -o OUTPUT, --output OUTPUT
                           directory in which to place the output files
     -l FLANKING, --flanking_length FLANKING
                        length of flanking regions, in bp; default: 1000
     -p MIN_PERCENT, --min_percent MIN_PERCENT
                        minimum percent identity
     -f ASSEMBLY_EXT, --assembly_ext ASSEMBLY_EXT
                        extenssion of reference, usually fasta
     -g REF_EXT, --ref_ext REF_EXT
                        extension of reference, usually .gb
     -F, --blast_Full      if true, blast full sequences along with just the
                        flanking. Interpretation is not implemented currently
                        as false positives cant be detected this way
     -v {1,2,3,4,5}, --verbosity {1,2,3,4,5}
                        Logger writes debug to file in output dir; this sets
                        verbosity level sent to stderr. 1 = debug(), 2 =
                        info(), 3 = warning(), 4 = error() and 5 = critical();
                        default: 2


Visualization
-------------

``riboSnag.py``
~~~~~~~~~~~~~~~

``riboSnag.py`` takes the list of clustered locus tags and extracts
their sequences with flanking regions, optionally turning the coding
sequences to N's to minimize bias towards reference. Is used to pull out
regions of interest from a Genbank file. Outputs a directory with a
fasta file for each clustered region (and a log file).

Additionally, it does a lot of plotting to visualize the Shannon
entropy, coverage, occurrences, and other useful metrics.

Usage:
^^^^^^

::

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

``riboStack.py``
~~~~~~~~~~~~~~~~

Decause assembly using short reads often collases rDNA repeats, it is
not uncommon to find a reference genome that has less than the actual
number of rDNAs. riboStack uses ``bedtools`` and ``samtools`` to
determine the coverage across rDNA regiosn, adn compares that coverage
depth to 10 sets of randomly selected non-rDNA regions. If the number of
rDNAs in the reference matches the number of rDNAs in your sequecned
isolate, the coverage should be pretty similar. However, if the coverage
in your rDNA regions is significantly higher, than there are likely more
rDNAs in your sequenced isoalte that there are in the reference, which
is something to be aware of.

It requires a mapping BAM file and the riboScan output directory as
input.


Utilities
-------------

``riboSwap.py``
~~~~~~~~~~~~~~~

Infrequently, ``riboSeed`` has joined together contigs that appear
incorrect according to your reference. If you are at all unhappy with a
bridging, ``riboSwap.py`` allows swapping of a "bad" contig for one or
more syntenic contigs from the *de novo* assembly. #### USAGE

::

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


``seedRand.py``
~~~~~

There is no convenient unix command to generate seeded random numbers from the
command line.  This script uses numpy (if availible) or the built-in random
module to generate *n* random numbers given a seed.

Note:  numpy *should* give you the same random numbers given the same seed
across platforms:  this is *not* the case with python's build-in random module.

::

   usage: seedRand.py [-h] seed n

   Given a seed, return a pseudrando integer between 1 and 9999, separated by
   newlines, to stdout. usage : `seedRand.py 27 10` would return 10 random
   numbers seeded with 27

   positional arguments:
     seed        seed
     n           number of random numbers to return, must be > 0

   optional arguments:
     -h, --help  show this help message and exit

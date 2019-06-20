Choosing an Appropriate Reference
===============

Introduction
----

Any sort of analysis that involves a reference strain begs the question: which strain should I use? Luckily, riboSeed's *de fere novo* approach gives you a bit of flexiblility here, in that the method is not affected my major genomic restructuring events, as long as they do not occur in the rDNA region.  Still, using a close reference maximizes your chance for a successful assembly.

There are three ways that we recommend selecting a reference.  The first uses average nucleotide identity to select an optimal reference for a given set of potential genomes -- perfect if you are sequencing a popular bug and are spoiled for reference choice.  The Kraken method gives great results when nothing is known a priori, and results in a high degree of certainty, but requires a bit of legwork. The last method can be done entirely through your web browser, but is much less robust. For popular bugs, we reccomend ANI; otherwise we recommend Kraken, but if you have already pulled out all of your hair from bioinformatics software installation, reads2Type/cgfind is a nice, painless alternative.


Reference Selection via ANI
---------------------------

With any assembly that uses a reference, the choice of that reference is
crucial. Here, we outline a protocol for using Average Nucleotide
Identity via pyANI.  We call this pipeline Plenty of Bugs, as it
will help you find a compatible match for your sequenced isolate.

`PlentyOfBugs found here <https://github.com/nickp60/plentyofbugs>`__


This is easiest to do with Docker (or Singularity)


The following script will identify 25 random E. coli genomes,
download then, build the pyani database, do a quick assembly of your isolate,
and find the closest reference to your isolate.

::

    docker run --rm -t -v  ${PWD}:/data/ nickp60/plentyofbugs:0.87 -f /data/tests/references/toy_reads1.fq -o "Escherichia coli" -n 25 -e sample_ecoli -d /data/results/



Kraken Method
-----


Kraken is a kmer-based phylogeny tool that can be used to idenify the strains present in a metagenomic dataset;  the installation and usage `instructions can be found here <https://ccb.jhu.edu/software/kraken/>`__

- Download and install Kraken, along with the MiniKraken database from their website.
- Run Kraken on your isolate's reads, and generate the Kraken report.

::

    kraken --db MiniKraken reads1.fq reads2.fastq  > sequences.kraken
    kraken-translate --db MiniKraken sequences.kraken > sequences.labels
    kraken-report --db MiniKraken sequence.kraken

Because the MiniKraken database was built from all the complete genomes from RefSeq, it should be easy to identify which strain in the database has the closest match to your sequenced isolate.

PS:  This is a great time to check if you have any contamination in your sample;  thanks, Kraken!


reads2type and  `cgFind <https://nickp60.github.io/cgfind>`__ Method
------


reads2type is also a kmer-based phylogeny tool, but it relies on a lightweight, prebuilt database that allows the analysis to be performed in your web browser, and it doesn't require you to upload your whole read file to a webserver.  It works by taking one read at a time from your file, generating 55-mers, and comparing to its databse. If there is not enough taxanomic information to indentify the isolate off of that read alone, additional reads will be processed until a single taxonomy is achieved.  This method works best on trimmed reads. `Instructions and the webserver can be found here <https://cge.cbs.dtu.dk/services/Reads2Type/>`__

Once you have a genus and species, you can use ``cgfind``, a tool we developed to provide easy access to downloadable genomes based on the complete prokaryotic genomes found in NCBI.  `it can be found here <https://nickp60.github.io/cgfind>`__  Just enter your genus and species name, and select one of the available strains to download.

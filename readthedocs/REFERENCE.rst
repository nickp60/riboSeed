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
Identity via pyANI.

All the subscripts are in the ``scripts/select_ref_by_ANI`` dir, and the
main runner script can be found at ``scripts/run_ANI.sh``

Required Tools
~~~~~~~~~~~~~~

``shuf``
^^^^^^^^

You must have coreutils installed on OSX to provide the ``gshuf``
replacement for ``shuf``.

``pyani``
^^^^^^^^^

``pyani`` makes it easy to perform average nucleotide identity
computations. Their development team has some exciting updates in the
works, so we will be working with the ``classify`` branch from GitHub.

::

    git clone https://github.com/widdowquinn/pyani.git
    cd pyani
    git checkout classify
    python setup.py install

``run_ANI.sh``: Helper Script Usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script helps automate this process. It takes 5 mandatory args:

-  ``-e`` experiement name
-  ``-f`` Forward Reads
-  ``-r`` Reverse Reads
-  ``-n`` number of strains
-  ``-o`` organism name

And 3 optional arguments. This allows you to reuse the same
prokaryotes.txt file from NCBI and the same comparison genomes (rather
than having to download them fresh each time). The last arg ``-a``
allows you to skip the assembly step and use an existing assembly for
ANI anaysis.

-  ``-p`` prokaryotes.txt file
-  ``-g`` directory of genomes of interest
-  ``-a`` existing assembly

This returns the string of the file name of the best reference, as well
as writes that to a file in the output directory called
``best_reference``.

What's it doing under the hood?

Description of Components
~~~~~~~~~~~~~~~~~~~~~~~~~

1. Randomly Selecting Potential References (if not using ``-g``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use the list of prokaryote assemblies from
`NCBI <ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt>`__.
If you have it downloaded already, you can pass it as an arg to the
get\_n\_genomes\_script.

The following script will get 25 random E. coli genomes

::

    ./scripts/select_ref_by_ANI/get_n_random_complete_genomes.sh -o "Escherichia coli" -n 25 > tmp_accessions

use the get\_genomes.py script to download all of those into a directory
called ``potential_references``

::

    mkdir potential_references
    while read x
    do
      echo "fetching $x"
      get_genomes.py -q $x -o potential_references
    done < tmp_accessions

2. Mini Assembly (if not using ``-a``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now, because we have a chicken/egg situation as we need an assembly to
determine ANI to determine the best reference to make an assembly, we
start with a minimal assembly with 1/10th of the reads using the
``mini_assembly.sh`` scsript.

::

    ./scripts/mini_assembly.sh ./tests/references/toy_reads1.fq ./tests/references/toy_reads2.fq TEST

    # copy the resulting assembly scaffold into the same dir as the other genomes we selected
    cp 2018-06-19_mini_assembly_TEST/spades/scaffolds.fasta ./potential_references/

3. Run pyANI
^^^^^^^^^^^^

::

    average_nucleotide_identity.py -i potential_references -g -o ./pyani

Look at the resulting ``pyani/ANIm_percentage_identity.tab`` file; the
best hit will be the one in the column/row for "contigs" with the
closest score to 1.




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

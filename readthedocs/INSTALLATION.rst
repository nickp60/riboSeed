Installation
============

From conda (new and recommended!)
---

Conda is a cross-platform, cross-language package management system.  If you haven't already installed conda, follow `these instructions here <https://bioconda.github.io/index.html>`__, and install the python3.6 version.  Once you have that done, add the appropriate channels.

::

   conda config --add channels defaults
   conda config --add channels conda-forge
   conda config --add channels bioconda


and then install riboSeed and all of its dependencies with one command:

::

   conda install riboseed


(Note the lowercase "s")


From Pypi
-----------

riboSeed is on Pypi, so you can install with pip, preferably within a
virtualenv:

::

    virtualenv -p python3.5 venv-riboSeed
    source venv-riboSeed/bin/activate
    pip3.5 install riboSeed

From TestPypi
-----------

To install the bleeding-edge version, install from testpypi:

::

    virtualenv -p python3.5 venv-riboSeed
    source venv-riboSeed/bin/activate
    pip install --extra-index-url https://testpypi.python.org/pypi riboSeed

From GitHub
-----------

You can also clone this repository, create a virtual environment, and run
``python3.5 setup.py install``.

..
   In the scripts directory, there is a script called ``runme.py`` which
   run the pipeline on a small sample dataset. it should output a folder to
   your current working dir called ``integration_tests``.

Dependencies
-----------

External Requirements
~~~~~~~~~~~~~~~~~~~~~

riboScan.py

-  Barrnap (must be 0.7 or above)
-  EMBOSS's Seqret

riboSelect.py

-  None

riboSnag.py

-  PRANK or Mafft
-  BLAST+ suite
-  Barrnap (must be 0.7 or above)

riboSeed.py

-  SPAdes v3.8 or higher
-  BWA (tested with 0.7.12-r1039)
-  SAMTools (must be 1.3.1 or above)
-  QUAST (tested with 4.1)

NOTE: barrnap has certain Perl requirements that may not be included on
your machine. Ensure barrnap runs fine before trying riboSeed. Or
try `python barrnap <https://github.com/nickp60/barrnap/>`__.

Installation
============

From Pypi (recommended)
-----------

riboSeed is on Pypi, so you can install with pip, preferably within a
virtualenv (recommended):

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

You can also clone this repository, and run
``python3.5 setup.py install``.

In the scripts directory, there is a script called ``runme.py`` which
run the pipeline on a small sample dataset. it should output a folder to
your current working dir called ``integration_tests``.

Dependencies
-----------

Python Requirements:
~~~~~~~~~~~~~~~~~~~~

-  Python >= v3.5
-  Biopython v1.68
-  pysam v0.9.1.4,
-  pyutilsnrw >= 0.0.768
-  matplotlib v1.5.3
-  pandas v0.18.1

External Requirements
~~~~~~~~~~~~~~~~~~~~~

riboScan.py

-  Barrnap (must be 0.7 or above)
-  Seqret

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
your machine. Ensure barrnap runs fine before trying ``riboSnag.py``. Or
try `python barrnap <https://github.com/nickp60/barrnap/>`__.

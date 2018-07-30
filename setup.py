"""
A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
import re
from codecs import open
from os import path
import sys

# see https://stackoverflow.com/questions/25192794
try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

here = path.abspath(path.dirname(__file__))

VERSIONFILE = "riboSeed/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


if sys.version_info <= (3, 0):
    sys.stderr.write("ERROR: riboseed requires Python 3.5 " +
                     "or above...exiting.\n")
    sys.exit(1)

## parse requirements file
install_reqs = parse_requirements("requirements.txt",
                                  session=False)
requirements = [str(ir.req) for ir in install_reqs]
setup(
    name='riboSeed',
    version=verstr,

    description='riboSeed: assemble across rDNA regions',
    # long_description=long_description,
    long_description="""
    its a blast, you gotta try it!  check out the GitHub
    repo for the read README.md file at https://github.com/nickp60/riboSeed
    """,
    url='https://github.com/nickp60/riboSeed',
    author='Nick Waters',
    author_email='nickp60@gmail.com',
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='bioinformatics assembly genomics development',
    packages=['riboSeed', 'tests'],
    install_requires=requirements,
    include_package_data=True,
    package_data={
       '': [path.join(__name__, "riboSeed", "integration_data/*")],
    },
    entry_points={
       'console_scripts': [
           'ribo=riboSeed.ribo:main',
       ],
    },
    scripts=[
        'scripts/seedRand.py',
        "scripts/runme.py"
    ],
)

#!/bin/bash
# script for installing samtools cause otherwise I cant get things 
#   to work with travis ci
wget 'https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2'
tar xf samtools-1.3.1.tar.bz2
cd samtools-1.3.1/
make
make prefix=/home/travis/ install
export samtools=/home/travis/bin/samtools

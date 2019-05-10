FROM continuumio/miniconda3
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install riboseed seqtk
RUN conda remove riboseed
RUN git clone https://github.com/nickp60/riboSeed
RUN git clone https://github.com/widdowquinn/pyani
RUN cd riboSeed && python setup.py install
RUN cd pyani && git checkout development && python setup.py develop
# test
RUN pyani
RUN pyani --help
RUN ribo --help
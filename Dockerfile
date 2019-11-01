FROM python:3.5
#  Future me, using a conda base docker image wasn't working
#  I tried using a miniconda docker container, I tried installing
#  miniconda within a python:3.5 container, etc. For some reason
#  every time things defaulted back to 3.7.  SPAdes did some strange
#  changes handling long reads between 3.9.1 and 3.10, so we are locked
#  in to 3.9, which locks us in to python 3.5, hence Docker
RUN wget  https://github.com/ncbi/SKESA/releases/download/v2.3.0/skesa.centos6.10 && mv skesa.centos6.10 skesa && chmod +x skesa && mv skesa /bin/
RUN apt-get update -y -qq
RUN apt-get install -y \
	ncbi-blast+ \
	bwa \
	mafft \
	bcftools \
	bedtools \
	samtools \
	libboost-all-dev \
	hmmer

RUN git clone https://github.com/tseemann/barrnap.git
ENV PATH="/barrnap/bin:${PATH}"

RUN mkdir /usr/spades
RUN cd /usr/spades/ && \
	wget http://cab.spbu.ru/files/release3.9.1/SPAdes-3.9.1-Linux.tar.gz && \
	tar xzvf SPAdes-3.9.1-Linux.tar.gz && \
	rm SPAdes-3.9.1-Linux.tar.gz
RUN pip install riboSeed
ENV PATH="/usr/spades/SPAdes-3.9.1-Linux/bin:${PATH}"
# RUN ribo --help
ENTRYPOINT [ "/usr/local/bin/ribo" ]

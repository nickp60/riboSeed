#!/bin/bash
echo "making a dir called `pwd`/riboSeedbin where we will install samtools,"
echo "  SPAdes, QUAST, barrnap, bwa, SRAtoolkit"
# make riboSeed bin dir
RB="`pwd`/riboSeedbin"
mkdir $RB
# move there
cd $RB
#Make backup of you bash profile in case of bad stuff
BACKUP_BP=$RB/.bash_profile_riboseed_backup
cp $HOME/.bash_profile $BACKUP_BP
echo "# The following things were added with the riboSeed installer"
# download EMBOSS for seqret
curl -L ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz | tar xz
cd EMBOSS-6.6.0
./configure --prefix=$RB
make
make install
cd ../



## Get SPAdes and QUAST, no instal needed
curl -L http://cab.spbu.ru/files/release3.9.1/SPAdes-3.9.1-Darwin.tar.gz | tar zx
curl -L https://downloads.sourceforge.net/project/quast/quast-4.4.tar.gz | tar zx

#echo "export PATH=$HOME/bin/SPAdes-3.9.1-Darwin:$PATH" >> $HOME/bash_profile
echo "export PATH=$RB/SPAdes-3.9.1-Darwin/bin:\$PATH" >> $HOME/.bash_profile
echo "export PATH=$RB/quast-4.4:\$PATH" >> $HOME/.bash_profile

# install samtools
curl -L https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj
cd samtools-1.3.1/
make
make prefix=$RB install
cd ../
echo "export PATH=$RB/samtools-1.3.1:\$PATH" >> $HOME/.bash_profile

# install barrnap
curl -L https://github.com/tseemann/barrnap/archive/0.7.tar.gz | tar zx
echo "export PATH=$RB/barrnap-0.7/bin:\$PATH" >> $HOME/.bash_profile


curl -L http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2 | tar xj
cd bwa-0.7.12
make
cd ../
echo "export PATH=$RB/bwa-0.7.12:\$PATH" >> $HOME/.bash_profile

curl -L https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1/sratoolkit.2.8.1-mac64.tar.gz | tar zx
echo "export PATH=$RB/sratoolkit.2.8.1-mac64/bin:\$PATH" >> $HOME/.bash_profile

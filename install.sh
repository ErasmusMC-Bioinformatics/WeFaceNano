#!/bin/bash
cd ~
# Install dependencies
apt-get update --fix-missing
apt-get -y install r-base perlbrew cpanminus ncbi-blast+ blast2 libgd-dev libbz2-dev cmake build-essential g++ qtbase5-dev libqt5svg5-dev libxml-libxml-perl
pip3 install django==2.2.8 biopython NanoPlot tabulate cgecore --user
# Bandage
wget https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip
unzip Bandage_Ubuntu_static_v0_8_1.zip
cd ~
# Porechop
git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install
cd ~
# Minimap
git clone https://github.com/lh3/minimap2
cd minimap2 && make
# Miniasm
git clone https://github.com/lh3/miniasm
cd miniasm && make
cd ~
# Flye
git clone https://github.com/fenderglass/Flye
cd Flye
python setup.py build
cd ~
# Racon
git clone --recursive https://github.com/isovic/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ~
# KmerGenie
wget http://kmergenie.bx.psu.edu/kmergenie-1.7051.tar.gz
tar -xvzf kmergenie-1.7051.tar.gz
cd kmergenie-1.7051 && make
cd ~
# Resfinder database
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
# Plasmidfinder
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
cd plasmidfinder
cd ~
# Simple Circularise
git clone https://github.com/Kzra/Simple-Circularise
# cpan packages
cpan -u
cpan install Getopt::Long Bio::SeqIO Bio::SearchIO Try::Tiny::Retry GD
# WeFaceNano
# git clone https://github.com/ErasmusMC-Bioinformatics/WeFaceNano

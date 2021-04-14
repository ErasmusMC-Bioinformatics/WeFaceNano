#!/bin/bash

# install simple-circularise
git clone https://github.com/Kzra/Simple-Circularise

# install kmergenie
wget http://kmergenie.bx.psu.edu/kmergenie-1.7051.tar.gz
tar xvzf kmergenie-1.7051.tar.gz
cd kmergenie-1.7051 && make

# cleanup
rm -rf ../kmergenie-1.7051.tar.gz*
cd ../

# get reference databases
mkdir -p nanopore-drive/blastdb
mv RB_REF/ nanopore-drive/
mv plasmidb/ nanopore-drive/
cd nanopore-drive
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db

cd plasmidfinder_db
./INSTALL.py
cd ../resfinder_db
./INSTALL.py
cd ../
# ln -s plasmidfinder_db plasmidb




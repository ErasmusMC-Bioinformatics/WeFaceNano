cd $HOME
git clone https://github.com/Kzra/Simple-Circularise
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db
wget http://kmergenie.bx.psu.edu/kmergenie-1.7051.tar.gz
tar -xvzf kmergenie-1.7051.tar.gz
cd $HOME/kmergenie-1.7051 && make
rm -rf $HOME/kmergenie-1.7051.tar.gz*
cd $HOME

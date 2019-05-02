# Plasmid and Resistance Identification by MinION ULtrafast (PRIMUL)

[Dependencies](#dependencies)

[installation](#installation)

[Usage](#usage)

[Results](#results)

## <a name="dependencies"></a>Dependencies

* Python 3.6 or higher
* perl5
* Biopython 1.73 or higher
* Django 2.1.2 or higher
* Albacore 2.3.4 or higher
* Porechop 0.2.4 or higher
* Canu 1.8 or higher
* Minimap 2.15-r915 or higher
* Miniasm 0.3-r179 or higher
* Racon 1.3.3 or higher
* Resfinder database
* Plasmidfinder
* NCBI-BLAST+ 2.6.0 or higher
* NCBI-BLAST 2 (version that includes blastall and formatdb needed to run resfinder)
* NanoPlot 1.20.1 or higher
* Simple-Circularise
* KmerGenie 1.7051
* cpan packages: Getopt::Long, Bio::SeqIO, Bio::SearchIO, Try::Tiny::Retry and GD

## <a name="installation"></a>How to install

1. Install dependencies:

```bash
# Python 3, Perl5, NCBI-BLAST+, NCBI-BLAST 2, libgd
sudo apt-get update
sudo apt-get install python3.6 perlbrew cpanminus ncbi-blast+ blast2 libgd-dev cmake

# Django, Biopython, NanoPlot
pip install django biopython NanoPlot --user

# PRIMUL
git clone https://github.com/ErasmusMC-Bioinformatics/PRIMUL
python PRIMUL/manage.py migrate

# Albacore 2.3.4
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.3.4-cp36-cp36m-manylinux1_x86_64.whl
sudo pip install ont_albacore-2.3.4-cp36-cp36m-manylinux1_x86_64.whl
rm ont_albacore-2.3.4-cp36-cp36m-manylinux1_x86_64.whl

# Porechop
git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install

# Minimap
git clone https://github.com/lh3/minimap2
cd minimap2 && make

# Miniasm
git clone https://github.com/lh3/miniasm
cd miniasm && make

# Canu 1.8
git clone https://github.com/marbl/canu.git
cd canu/src
make -j <number of threads>

# Racon
git clone --recursive https://github.com/isovic/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make

# KmerGenie
wget http://kmergenie.bx.psu.edu/kmergenie-1.7051.tar.gz
tar -xvzf kmergenie-1.7051.tar.gz
cd kmergenie-1.7051 && make

# Resfinder database
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git

# Plasmidfinder
cd /path/to/some/dir
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
cd plasmidfinder

# Simple Circularise
git clone https://github.com/Kzra/Simple-Circularise

# cpan packages
cpan install Getopt::Long, Bio::SeqIO, Bio::SearchIO, Try::Tiny::Retry and GD
```

2. Mount a network drive to /mnt/d/Nanopore to start using PRIMUL with the default settings. You can change the default nanopore drive in the settings.py by changing the NANOPORE_DRIVE variable.

3. Create a django superuser and start the PRIMUL server:

```bash
django-admin createsuperuser --username username --password password --database path/to/PRIMUL/db.sqlite3
python path/to/PRIMUL/manage.py runserver 0.0.0.0:8080
```

4. Go to 127.0.0.1:8080 to see if the server is running and the homepage is visible. If the homepage is visible you can now log in and start using the plasmid pipeline.

## <a name="usage"></a>Plasmid Pipeline Usage

To start using the Nanopore pipeline a user account has to be created by a superuser.

If you are logged in as a superuser, a new user can be created by clicking the sign up button and entering a username and password. When the new user is created a new folder will appear on the network drive with the username. This drive can be mounted on the users local machine and new data can be added using that mounted drive.
The default drive will be /mnt/d/Nanopore. In the Nanopore folder, the data structure should be as follows:

```text
[username]
    |
    |---> Run01 (barcoded)
    |       |
    |       |---> BC01
    |       |      |---> reads01.fastq
    |       |      |---> reads02.fastq
    |       |      |---> ...
    |       |---> BC02
    |       |      |---> reads01.fastq
    |       |      |---> reads02.fastq
    |       |      |---> ...
    |       |---> BC..
    |              |---> ...
    |
    |---> Run02 (not barcoded FASTQ)
    |       |
    |       |---> fastq
    |               |---> reads01.fastq
    |               |---> reads02.fastq
    |               |---> ...
    |
    |---> Run03 (not barcoded FAST5)
    |       |
    |       |---> fast5
    |               |---> read01.fast5
    |               |---> read02.fast5
    |               |---> ...
    |
    |---> Run ...
```

Click on the pipelines button to start. Select FAST5 is you want to basecall using Albacore 2. If you want to skip the basecalling because you already have the FASTQ files, please select the FASTQ option. Select the project folder and enter a name for your output folder. This folder will be created in the results folder on the network drive. Using the default settings the results will be stored in this folder: <b>/mnt/d/Nanopore/[username]/results/</b>. When using a different path the results will be stored in that location.

If the FAST5 option is selected, please select the following Albacore settings:

* Barcoding [yes or no]
* Flowcell configuration file.

When selecting the FASTQ option the Albacore basecalling step will be skipped and Canu or Miniasm will use the inputfolder to start the assembly. For the Canu assembly, please enter the genome size. Try to guess if you don't know the exact size. For Miniasm this step can be skipped but instead the KmerGenie tool will run to find the optimal kmer-size for the selected reads. If the run is barcoded the tool will run for all barcodes and will calculate the optimal kmer-size for all barcodes.

After the assembly the Simple-Circularise tool will look for repeats and will circularise the contigs.

The plasmid pipeline can run two additional tools (BLAST and resfinder). Please check the tools you want to run during the pipeline.

When the BLAST option is selected, please select the database you want to use (It is possible to use the remote nt/nr database, a local nt/nr database or a local plasmid database). To use a local database please make sure the database is on the mounted network drive in the following path: <b>/mnt/d/Nanopore/plasmidb/</b>

<b>*If the default drive path is changed make sure the database is in the correct location.</b>

If the resfinder option is selected, please select the resfinder settings you want to use to find the antibiotic resistance genes.
The options are the threshold % (80%, 85%, 90%, 95%, 100%) and the minimum length of the overlap (60%, 70%, 80%, 90%, 100%).
There is also an option to search for a specific gene or search for all available genes within resfinder.

## <a name="results"></a>View Pipeline Results
When the pipeline is finished running you will be send back to the homepage. On this page you will see a dropdown menu with all generated resultfolders for the user that is logged in. To see any of the generated results select the result you want to view and click the VIEW button. After a few minutes/seconds a resultpage will be shown with the following information:

1. Plots created with NanoPlot (readlength and yield histogram)
2. A table with the assembly information containing the contigs and lengths of the contigs.
3. A table with the antimicrobial resistance genes found by resfinder for each of the contigs.
4. The BLAST results with the top hit for all contigs, percentage identity and alignment length
5. Interactive circle diagrams based on the BLAST results.

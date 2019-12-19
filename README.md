# WeFaceNano

[Dependencies](#dependencies)

[installation](#installation)

[Usage](#usage)

[Results](#results)

## <a name="dependencies"></a>Dependencies

* Python 3.6 or higher
* perl5
* R 3.6.1 or higher
* Biopython 1.73 or higher
* Django 2.2.8
* Bandage 0.8.1 or higher
* Porechop 0.2.4 or higher
* Minimap 2.15-r915 or higher
* Miniasm 0.3-r179 or higher
* Flye 2.4.2 or higher
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

### Install from source

1. Install dependencies by running the install.sh script.

2. Clone the latest version of WeFaceNano from GitHub [```git clone https://github.com/ErasmusMC-Bioinformatics/WeFaceNano```]

3. Add the static folder to PATH.

4. Mount a network drive to the Nanopore path (/mnt/d/Nanopore) to start using WeFaceNano with the default settings. You can change the default nanopore drive in the settings.py by changing the NANOPORE_DRIVE variable.

5. Start the WeFaceNano server by running the wefacenano_start.sh script.

6. Log in using admin.

7. Go to 127.0.0.1:8008 to see if the server is running and the homepage is visible. If the homepage is visible you can now log in and start using the plasmid pipeline.

### Build docker

1. Clone the latest version of WeFaceNano from GitHub [```git clone https://github.com/ErasmusMC-Bioinformatics/WeFaceNano```].

2. Open the command line (Windows) or terminal (Linux).

3. Build the docker: ```docker build -f wefacenano.Dockerfile -t wefacenano .```

4. Create a Nanopore working folder and add the blast database and plasmid database to this folder. The BLAST database folder should have the name blastdb with the database name called nt and nr. The plasmid folder and database should have the name plasmidb.

5. Run the docker: ```docker run -p 8008:8008 -v /path/to/Nanopore_Folder:/mnt/d/Nanopore -it wefacenano```.

6. Log in using admin.

7. Go to 127.0.0.1:8008 to see if the server is running and the homepage is visible. If the homepage is visible you can now log in and start using the plasmid pipeline.

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
    |---> Run ...
```

Click on the pipelines button to start. Select FAST5 is you want to basecall using Albacore 2. If you want to skip the basecalling because you already have the FASTQ files, please select the FASTQ option. Select the project folder and enter a name for your output folder. This folder will be created in the results folder on the network drive. Using the default settings the results will be stored in this folder: <b>/mnt/d/Nanopore/[username]/results/</b>. When using a different path the results will be stored in that location.

When selecting the FASTQ option the Albacore basecalling step will be skipped and Flye or Miniasm will use the inputfolder to start the assembly. For the Flye assembly, please enter the genome size. Try to guess if you don't know the exact size. For Miniasm this step can be skipped. If the run is barcoded the tool will run for all barcodes and will calculate the optimal kmer-size for all barcodes. After the Miniasm assembly the Simple-Circularise tool will look for repeats and will try to circularise the contigs.

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

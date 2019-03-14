# Plasmid and Resistance Identification by MinION ULtrafast (PRIMUL)

[Dependencies](#dependencies)

[installation](#installation)

[Usage](#usage)

[Results](#results)

## <a name="dependencies"></a>Dependencies

* Python 3.5 or higher
* perl5
* Biopython 1.70 or higher
* Django 2.0 or higher
* Albacore 2.3.4 or higher
* Canu 1.8 or higher
* Minimap 2.15-r915 or higher
* Miniasm 0.3-r178 or higher
* Resfinder database
* BLAST+ 2.6.0 or higher
* BLAST 2 (version that includes blastall and formatdb needed to run resfinder)
* NanoPlot 1.20.1 or higher
* Simple-Circularise
* libgd 2.2.5
* cpan packages: Getopt::Long, Bio::SeqIO, Bio::SearchIO, Try::Tiny::Retry and GD

## <a name="installation"></a>How to install

1. Get the newest version of PRIMUL by cloning of downloading the code from [github](https://github.com/ErasmusMC-Bioinformatics/PRIMUL) and
install all dependencies.

2. Download the Resfinder database

```bash
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
cd resfinder_db
```

3. Mount a network drive to /mnt/d/Nanopore to start using PRIMUL with the default settings. You can change the default nanopore drive in the settings.py by changing the NANOPORE_DRIVE variable.

4. Create a django superuser and start the PRIMUL server:

```bash
django-admin createsuperuser --username username --password password --database path/to/PRIMUL/db.sqlite3
python path/to/PRIMUL/manage.py runserver 0.0.0.0:8080
```

5. Go to 127.0.0.1:8080 to see if the server is running and the homepage is visible. If the homepage is visible you can now log in and start using the plasmid pipeline.

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

When selecting the FASTQ option the Albacore basecalling step will be skipped and Canu or Miniasm will use the inputfolder to start the assembly. For the Canu assembly, please enter the genome size. Try to guess if you don't know the exact size. For Miniasm this step can be skipped. After the assembly the Simple-Circularise tool will look for repeats and will circularise the contigs.

The plasmid pipeline can run two additional tools (BLAST and resfinder). Please check the tools you want to run during the pipeline.

When the BLAST option is selected, please select the database you want to use (It is possible to use the remote nt/nr database, a local nt/nr database or a local plasmid database). To use a local database please make sure the database is on the mounted network drive in the following path: <b>/mnt/d/Nanopore/plasmidb/</b>

<b>*If the default drive path is changed make sure the database is in the correct location.</b>

If the resfinder option is selected, please select the resfinder settings you want to use to find the antibiotic resistance genes.
The options are the threshold % (80%, 85%, 90%, 95%, 100%) and the minimum length of the overlap (60%, 70%, 80%, 90%, 100%).
There is also an option to search for a specific gene or search for all available genes within resfinder.

## <a name="results"></a>View Pipeline Results
When the pipeline is finished running you will be send back to the homepage. On this page you will see a dropdown menu with all generated resultfolders for the user that is logged in. To see any of the generated results select the result you want to view and click the VIEW button. After a few minutes/seconds a resultpage will be shown with the following information:

1. Plots created with NanoPlot (readlength and yield histogram)
2. The assembly information containing the contigs and lengths of the contigs.
3. Antimicrobial resistance genes found by resfinder for each of the contigs.
4. The BLAST results with the top hit for all contigs, percentage identity and alignment length
5. Interactive circle diagrams based on the BLAST results.

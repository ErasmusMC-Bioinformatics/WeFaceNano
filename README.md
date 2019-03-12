# PRIMUL
Nanopore plasmid pipeline.

[Dependencies](#dependencies)

[installation](#installation)

[Usage](#usage)

[Results](#results)

## <a name="dependencies"></a>Dependencies

* Python 3.5 or higher
* perl5
* libgd 2.2.5
* cpan packages: Getopt::Long, Bio::SeqIO, Bio::SearchIO, Try::Tiny::Retry and GD
* Biopython 1.70 or higher
* Django 2.0 or higher
* Albacore 2.1.0 or higher
* Canu 1.6 or higher
* Minimap 2.2-r519-dirty or higher
* Miniasm 0.2-r168-dirty or higher
* Resfinder 2.1 (edited version of resfinder.pl is available in the static folder)
* BLAST+ 2.7.1 or higher
* MinionQC
* Simple-Circularise

## <a name="installation"></a>How to install

Get the newest version of PRIMUL by cloning of downloading the code from [github](https://github.com/ErasmusMC-Bioinformatics/PRIMUL) and
install all dependencies as seen in the [dependencies list](#dependencies).

Please mount a network drive to /mnt/d/Nanopore in order for PRIMUL to work without having to make any changes to the code.
If a different storage option will be used changes will have to be made in the PRIMUL code. To do this, open views.py and change the NANOPORE_DRIVE variable.

To run the PRIMUL server enter:

```bash
python path/to/PRIMUL/manage.py runserver 0.0.0.0:8080
```

Go to 127.0.0.1:8080 to see if the server is running and the homepage is visible. If the homepage is visible you can now log in (or create a user) and start using the plasmid pipeline.

To add a superuser please use the following command in the terminal: 

```bash
django-admin createsuperuser --username username --password password --database /path/to/database
```

## <a name="usage"></a>Plasmid Pipeline Usage

To start using the Nanopore pipeline a user account has to be created by a superuser.

If you are logged in as a superuser, a new user can be created by clicking the sign up button and entering a username and password. When the new user is created a new folder will appear on the network drive with the username. This drive can be mounted on the users local machine and new data can be added using that mounted drive.

When using FASTQ data all FASTQ files must be added to a folder called fastq, this folder should be added within the project folder. If FAST5 files are being used they have to be added to a folder called fast5 within the project folder.
The FASTQ structure should look like this:

<b>Run Folder ---> Barcode Folder ---> FASTQ, FAST5 or FASTA</b>

The FAST5 structure should look like this:

<b>Run01 ---> BC01 ---> read01.fastq, read02.fastq</b>

Click on the pipelines button to start selecting a pipeline to run. For now the only available pipeline is the plasmid pipeline.

Select the input filetype (FAST5 or FASTQ). Select FAST5 is you want to basecall using Albacore 2.1.10. 
To skip the basecalling because you already have the FASTQ files, please select FASTQ. Select the project folder with the FAST5 or FASTQ files and enter a name for your output folder. This folder will be created in the results folder on the network drive.

If the FAST5 option is selected, please select the following Albacore settings:

* Barcoding [yes or no]
* Flowcell configuration file.

When selecting the FASTQ option the Albacore basecalling step will be skipped and Canu or Miniasm will use the inputfolder to start the assembly. For the Canu assembly, please enter the genome size. Try to guess if you don't know the exact size. For Miniasm this step can be skipped.

The plasmid pipeline can run two additional tools (BLAST and resfinder). Please check the tools you want to run during the pipeline.

When the BLAST option is selected, please select the database you want to use and which BLAST task (megablast or blastn) you want to use. There are options for a local blast database and remote databases.

If the resfinder option is selected, please select the resfinder settings you want to use to find the antibiotic resistance genes.
The options are the threshold % (80%, 85%, 90%, 95%, 100%) and the minimum length of the overlap (60%, 70%, 80%, 90%, 100%).
There is also an option to search for a specific gene or search for all available genes within resfinder.

The results will be stored in the network drive. A results folder will be created within the user folder to store all future results created by PRIMUL.

## <a name="results"></a>View Pipeline Results
When the pipeline is finished running you will be send back to the homepage. On this page you will see a dropdown menu with all generated resultfolders for the user that is logged in. To see any of the generated results select the result you want to view and click the VIEW button. After a few minutes/seconds a resultpage will be shown with the following information:

1. The assembly information containing the contigs and lengths of the contigs.
2. Antimicrobial resistance genes found by resfinder for each of the contigs.
3. The BLAST results with the top hit for all contigs, percentage identity and alignment length
4. Plasmid circle diagrams based on the BLAST results.

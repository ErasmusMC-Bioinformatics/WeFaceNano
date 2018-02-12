# RISPIC
Nanopore plasmid pipeline.

[Dependencies](#dependencies)

[installation](#installation)

[Usage](#usage)

[Results](#results)

# <a name="dependencies"></a>Dependencies
* Python 2.7.12
* perl5
* libgd 2.2.5
* cpan packages: Getopt::Long, Bio::SeqIO, Bio::SearchIO, Try::Tiny::Retry and GD
* Biopython 1.70
* Django 1.11.3
* Albacore 2.1.10
* Canu 1.6
* Resfinder 2.1 (edited version of resfinder.pl is available in the static folder)
* BLAST+ 2.7.1
* blastall and formatdb (both are available in blast2)
* MinionQC

# <a name="usage"></a>How to install
Get the newest version of RISPIC by cloning of downloading the code from [github](https://github.com/ErasmusMC-Bioinformatics/RISPIC) and
install all dependencies as seen in the [dependencies list](#dependencies).

Please mount a network drive to /media/Nanopore in order for RISPIC to work without having to make any changes to the code.
If a different storage option will be used changes will have to be made in the RISPIC code.
results/results/
To run the RISPIC server enter path/to/RISPIC/manage.py runserver 0.0.0.0:8080

Go to 127.0.0.1:8080 to see if the server is running and the homepage is visible. If the homepage is visible 
you can now log in (or create a user) and start using the plasmid pipeline.

To add a superuser please use the following command in the terminal: django-admin createsuperuser --username username --password password --database /path/to/database 

# <a name="usage"></a>Plasmid Pipeline Usage
To start using the Nanopore pipeline a user account has to be created by a superuser.

A new user can bew created by clicking the sign up button and entering a username and password. 
A new folder will be created on the network drive with the username. This drive can be synced and accessed with the winscp software.

Click on the pipelines button to start selecting a pipeline to run. For now the only available 
pipeline is the plasmid pipeline.

Select the input filetype (FAST5 or FASTQ). Select FAST5 is you want to basecall using Albacore 2.1.10. 
To skip the basecalling because you already have FASTQ files created by Metrichor, please select FASTQ.
Enter a name for your output folder. This folder will be created in the results folder on the network drive.

If the FAST5 option is selected, please select the following Albacore settings:
* Barcoding [yes or no]
* Flowcell configuration file.

When selecting the FASTQ option the Albacore basecalling step will be skipped and Canu will use the inputfolder to start the assembly.

For the Canu assembly, please enter the genome size. Try to guess if you don't know the exact size.

The plasmid pipeline can run two additional tools (BLAST and resfinder). Please check the tools you want to run during the pipeline.

When the BLAST option is selected, please select the database you want to use 
and which BLAST task (megablast or blastn) you want to use. There are options for a local blast database and remote databases.

If the resfinder option is selected, please select the resfinder settings you want to use to find the antibiotic resistance genes.
The options are the threshold % (80%, 85%, 90%, 95%, 100%) and the minimum length of the overlap (60%, 70%, 80%, 90%, 100%).
There is also an option to search for a specific gene or search for all available genes within resfinder.

The results will be stored in the network drive. A results folder will be created within the user folder to store all future results created by RISPIC.

# <a name="results"></a>View Pipeline Results
When the pipeline is finished running a results page will be generated. This page contains the BLAST results from all barcodes and contigs.
The results are sorted by barcodes and each barcode has a table with three columns:
1. The plasmid / accession number and name found by BLAST and length in bp from the genbank information. 
2. The contig name and length.
3. Antimicrobial resistance genes found by resfinder
 

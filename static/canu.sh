#!/usr/bin/env bash

echo
VERSION=0
SEED=42
# Get options
while getopts p:d:g:i:hs: option
do
 case "${option}"
 in
 p) PREFIX=${OPTARG};;
 i) FASTQ=${OPTARG};;
 d) OUTFOLDER=${OPTARG};;
 g) GENOMESIZE=${OPTARG};;
 h) echo "Usage: bash canu.sh -p <prefix filename> -d <output folder> -d <FASTQ folder> --- Pipeline that runs canu for long read assembly."
    echo "Run by selecting"
    echo
    exit 1;;
 s) SEED=${OPTARG};;
 :) echo "missing argument for -%s\n" "$OPTARG" >&2
       exit 1
       ;;
 esac
done

~/canu/Linux-amd64/bin/canu -p ${PREFIX} -d ${OUTFOLDER} genomeSize=${GENOMESIZE} -nanopore-raw ${FASTQ}/*.fastq stopOnReadQuality=false gnuplotTested=true

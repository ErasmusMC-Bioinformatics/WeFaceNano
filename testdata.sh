#!/bin/bash

mkdir -p nanopore-drive/admin/testrun/fasta/
mkdir -p nanopore-drive/admin/testrun-small/RB01/
mkdir -p nanopore-drive/admin/testrun-small/RB02/
cd nanopore-drive/admin/testrun/fasta/

wget ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100387/MinION_reads-12_samples.rar
unrar e MinION_reads-12_samples.rar .

# create a folde with a subset of the data for speed
cp RB01.fasta ../../testrun-small/RB01/
cp RB02.fasta ../../testrun-small/RB02/

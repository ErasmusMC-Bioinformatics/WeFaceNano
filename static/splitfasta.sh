#!/usr/bin/env bash

i=0;
FOLDER=$2;

while read line ; do
  if [ ${line:0:1} == ">" ] ; then
    ((i++))
    echo "$line" >> ${FOLDER}/unitig"${i}".fasta
  else
    echo "$line" >> ${FOLDER}/unitig"${i}".fasta
  fi
done < $1
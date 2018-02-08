#!/bin/bash

echo
VERSION=0
SEED=42
# Get options
while getopts v:i:o:hs: option
do
 case "${option}"
 in
 v) VERSION=${OPTARG};;
 i) FASTQFile=${OPTARG};;
 o) FASTAOUT="${OPTARG}.fasta";;
 h) echo "Usage: miniasm.bash -v <MiniMap Version number> -i <FASTQ input> -o <FASTA output> --- Pipeline that runs MiniMap and Miniasm and creates a FASTA file."
    echo
    exit 1;;
 s) SEED=${OPTARG};;
 :) echo "missing argument for -%s\n" "$OPTARG" >&2
       exit 1
       ;;
 esac
done

GAFOUTPUT="$FASTQFile.gfa"
TEMPPAF="$FASTQFile.paf.gz"
#FASTAOUT="$FASTQFile.out.fasta"
#MAPLOG="$FASTQFile.minimap.log"
#ASMLOG="$FASTQFile.miniasm.log"
echo -n "-Starting MiniMap / Miniasm pipeline at: "
date
# Cleanup from previous runs:
if [ -f "$TEMPPAF" ]; then
	echo "Trying to removing any existing files - if they don't exist expect an 'rm' error"
	rm $TEMPPAF
	rm $GAFOUTPUT
	rm $FASTAOUT
	#rm $MAPLOG
	#rm $ASMLOG
fi
echo
### 1) MiniMap
case "$FASTQFile"
in
'') echo "No or invalid FASTQ file"
	exit 1;;
esac
echo "-Running MiniMap with input file '$FASTQFile'"
echo -n "MiniMap Start: "
date
echo
case "$VERSION"
in
0) echo "Select a MiniMap version (-v1 or -v2)"
	exit 1;;
1) time ~/minimap/minimap -x ava10k "$FASTQFile" "$FASTQFile" | gzip -1 > "$TEMPPAF";;
2) time ~/minimap2/minimap2 -x ava-ont "$FASTQFile" "$FASTQFile" | gzip -1 > "$TEMPPAF";;
esac
echo
echo -n "MiniMap2 End: "
date
echo
### 2) Miniasm
echo -n "Miniasm Start: "
date
time ~/miniasm/miniasm -f "$FASTQFile" "$TEMPPAF" > "$GAFOUTPUT" #>(tee $ASMLOG >&2) 
echo -n "Miniasm End: "
date
echo
### 3) Convert to FASTA
awk '/^S/{print ">"$2"\n"$3}' "$GAFOUTPUT" | fold > "$FASTAOUT"
echo
### 4) Summary:
echo "---------------Summary---------------"
echo " "
echo " MINIMAP Version	: $VERSION"
echo " FASTQ Input   		: $FASTQFile"
echo " FASTA Output  		: $FASTAOUT"
echo -n " No. of utigs		: "
grep ">" "$FASTAOUT" | wc -l  | cut -f1
echo " "
echo "---------------Summary---------------"
echo 
echo -n "Removing $GAFOUTPUT ..... "
rm $GAFOUTPUT
echo "$GAFOUTPUT removed"
echo -n "Removing $TEMPPAF ..... "
rm $TEMPPAF
echo "$TEMPPAF removed"
echo
echo -n "-Finishing MiniMap / Miniasm pipeline at: "
date 
echo

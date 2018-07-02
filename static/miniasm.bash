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
GFAOUTPUT="$FASTQFile.gfa"
TEMPPAF="$FASTQFile.paf.gz"
echo -n "-Starting MiniMap / Miniasm pipeline at: "
date
# Cleanup from previous runs:
if [ -f "$TEMPPAF" ]; then
	echo "Trying to removing any existing files - if they don't exist expect an 'rm' error"
	rm $TEMPPAF
	rm $GFAOUTPUT
	rm "$FASTQFile.racon.fasta"
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
1) time ~/minimap/minimap -Sw5 -L100 -m0 "$FASTQFile" "$FASTQFile" | gzip -1 > "$TEMPPAF";;
2) time ~/minimap2/minimap2 -x ava-ont "$FASTQFile" "$FASTQFile" | gzip -1 > "$TEMPPAF";;
esac
echo
echo -n "MiniMap End: "
date
echo
### 2) Miniasm
echo -n "Miniasm Start: "
date
time ~/miniasm/miniasm -f "$FASTQFile" "$TEMPPAF" > "$GFAOUTPUT"
echo -n "Miniasm End: "
date
echo
echo -n "Create concensus"
awk '$1 == "S" {print ">"$2;print $3}' $GFAOUTPUT > "$FASTQFile.raw.fasta"
echo
echo "Concensus build"
echo
echo "Run Minimap against the concensus"
echo
time ~/minimap2/minimap2 "$FASTQFile.raw.fasta" "$FASTQFile" > "$FASTQFile.mapped.paf"
echo
echo "Run racon"
echo
time racon -t8 "$FASTQFile" "$FASTQFile.mapped.paf" "$FASTQFile.raw.fasta" > "$FASTQFile.racon.fasta"
echo
echo "---------------Summary---------------"
echo " "
echo " MINIMAP Version	: $VERSION"
echo " FASTQ Input   		: $FASTQFile"
echo " FASTA Output  		: $FASTQFile.racon.fasta"
echo -n " No. of utigs		: "
grep ">" "$FASTQFile.racon.fasta" | wc -l  | cut -f1
echo " "
echo "---------------Summary---------------"
echo
echo -n "Filter FASTA sequences by length"
echo
awk '!/^>/ { next } { getline seq } length(seq) >= 15000 { print $0 "\n" seq }' "$FASTQFile.racon.fasta" > "$FASTAOUT"
echo
echo "No. of utigs after filtering:"
grep ">" "$FASTAOUT" | wc -l  | cut -f1
echo
echo -n "Removing $GAFOUTPUT ..... "
rm "$GFAOUTPUT"
echo "$GFAOUTPUT removed"
echo -n "Removing $TEMPPAF ..... "
rm "$TEMPPAF"
echo "$TEMPPAF removed"
echo -n "Removing $FASTQFile.raw.fasta ..... "
rm "$FASTQFile.raw.fasta"
echo -n "Removing $FASTQFile.mapped.paf ..... "
rm "$FASTQFile.mapped.paf"
echo -n "Removing temp racon file ..... "
rm "$FASTQFile.racon.fasta"
echo -n "-Finishing MiniMap / Miniasm pipeline at: "
date
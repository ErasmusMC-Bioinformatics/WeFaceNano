#!/usr/bin/env bash

BARCODING=0
MLST=0
RES=0
SPECIES=0

while getopts bmrsh:-c:-i:-o:-f: option
do
#key="$1"
case "${option}" in
    b|barcoding)
    BARCODING=1
#     shift # past argument
#     shift # past value
    ;;
    c|config)
    CONFIG=${OPTARG}
#     shift # past argument
#     shift # past value
    ;;
    i|input)
    FAST5Files=${OPTARG}
#     shift # past argument
#     shift # past value
    ;;
    o|output)
    FILEOUT=${OPTARG}
#     shift # past argument
#     shift # past value
    ;;
    f|format)
    FORMAT=${OPTARG}
#     shift # past argument
#     shift # past value
    ;;
    m|mlst)
    MLST=1
#     shift # past argument
#     shift # past value
    ;;
    r|res)
    RES=1
#     shift # past argument
#     shift # past value
    ;;
    s|species)
    SPECIES=1
#     shift # past argument
#     shift # past value
    ;;
    h|help) 
    echo "bash npanalysis (-b -m -s -r) -c -i -o -f"
    echo 
    echo "-c|--config:"
    echo "Select the config file that corresponds with the used flowcell and kit combination"
    echo "i.e. -c r94_450bps_linear.cfg"
    echo 
    echo "flowcell      kit         barcoding   config file"
    echo "FLO-MIN106    SQK-DCS108  -           r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-LSK108  -           r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-LWB001  included    r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-PCS108  -           r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-RAB201  included    r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-RAD002  -           r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-RAD003  -           r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-RAS201  -           r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-RBK001  included    r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-RLB001  included    r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-RLI001  -           r94_450bps_linear.cfg"
    echo "FLO-MIN106    SQK-RNA001  -           r94_70bps_rna_linear.cfg"
    echo "FLO-MIN106    VSK-VBK001  -           r94_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-DCS108  -           r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-LSK108  -           r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-LWB001  included    r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-LWP001  -           r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-PCS108  -           r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-RAB201  included    r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-RAD002  -           r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-RAD003  -           r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-RAS201  -           r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-RBK001  included    r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-RLB001  included    r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-RLI001  -           r95_450bps_linear.cfg"
    echo "FLO-MIN107    SQK-RNA001  -           r94_70bps_rna_linear.cfg"
    echo "FLO-MIN107    VSK-VBK001  -           r95_450bps_linear.cfg"
    echo 
    echo "-i|--input: "
    echo "input folder with fast5 reads - please enter the folder with the fast5 files"
    echo 
    echo "-o|--output: "
    echo "Albacore basecaller output folder"
    echo 
    echo "-b|--barcoding: "
    echo "Add -b or --barcoding to run Albacore in barcoding mode."
    exit 1
    ;;
esac
done

# Albacore
if [ "${BARCODING}" = 1 ]; then
      gnome-terminal -e 'bash -c "read_fast5_basecaller.py --barcoding -r -i '"${FAST5Files}"' -c '"${CONFIG}"' -t $(nproc) -o '"${FORMAT}"' -s '"${FILEOUT}"' 2> /dev/null"' &
      sleep 5
      while (ps -a | grep -c read_fast5); do
        echo "Albacore is running..."
      done
      FINISHED=${FILEOUT}/finished.dat
      touch ${FINISHED}
else
      gnome-terminal -e 'bash -c "read_fast5_basecaller.py -r -i '"${FAST5Files}"' -c '"${CONFIG}"' -t $(nproc) -o '"${FORMAT}"' -s '"${FILEOUT}"'"' &
      sleep 5
      while (ps -a | grep -c read_fast5); do
        echo "Albacore is running..."
      done
      FINISHED=${FILEOUT}/finished.dat
      touch ${FINISHED}
fi
# speciesTyping
#if [ "${SPECIES}" = 1 ]; then
#      gnome-terminal -e 'bash -c "jsa.util.streamServer -port 3456 | bwa mem -t 10 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 ~/japsa/SpeciesTyping/Bacteria/genomeDB.fasta  - 2> /dev/null | jsa.np.rtSpeciesTyping -bam - -index ~/japsa/SpeciesTyping/Bacteria/speciesIndex --read 5 -time 10 -out '"${FILEOUT}"'/species/speciesTypingResults.out"' &
#      sleep 5
#      gnome-terminal -e 'bash -c "jsa.np.npreader -realtime -streams 0.0.0.0:3456 -folder '"${FILEOUT}"'/workspace/pass --output - -format fastq --exhaustive"'
#fi
# MLST
#if [ "${MLST}" = 1 ]; then
#      gnome-terminal -e 'bash -c "jsa.util.streamServer --port 3457 | bwa mem -t 8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y ~/japsa/MLST/Klebsiella_pneumoniae/bwaIndex/genes.fasta - | jsa.np.rtMLST -bam - -mlst ~/japsa/MLST/Klebsiella_pneumoniae/ -read 5 -time 10  --out '"${FILEOUT}"'/MLST.dat"' &
#      sleep 5
#      gnome-terminal -e 'bash -c "jsa.np.npreader -realtime -streams 0.0.0.0:3457 -folder '"${FILEOUT}"'/workspace/pass --output - -format fastq --exhaustive"'
#fi
# resFinder
#if [ "${RES}" = 1 ]; then
#      gnome-terminal -e 'bash -c "jsa.util.streamServer --port 3458 | bwa mem -t 2 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y ~/japsa/ResGene/resFinder/DB.fasta - 2> /dev/null | jsa.np.rtResistGenes -bam - -score=0.00001 -read 5 -time 10 -resDB ~/japsa/ResGene/resFinder -tmp /tmp/_tmp_ -o '"${FILEOUT}"'/resfinder/resFinder.out"' &
#      sleep 60
#      gnome-terminal -e 'bash -c "jsa.np.npreader -realtime -streams 0.0.0.0:3458 -folder '"${FILEOUT}"'/workspace/pass --output - -format fastq --exhaustive"'
#      gnome-terminal -e 'bash -c "jsa.np.npreader --realtime --folder='"${FILEOUT}"'/workspace/pass --output=- --exhaustive | bwa mem -t12 -x ont2d -K10000 /home/myfair/japsa/ResGene/resFinder/DB.fasta - 2> /dev/null | jsa.np.rtResistGenes --output='"${FILEOUT}"'/resfinder/resfinder.out --bamFile=- --resDB=/home/myfair/japsa/ResGene/resFinder --time=60 --tmp=tmp/resTest 2> '"${FILEOUT}"'/res.log"'
#fi

# npReader
#gnome-terminal -e 'bash -c "jsa.np.npreader -realtime -streams 0.0.0.0:3456,0.0.0.0:3457,0.0.0.0:3458 -folder '"${FILEOUT}"'/workspace/pass --output output.fastq -format fastq --exhaustive"'
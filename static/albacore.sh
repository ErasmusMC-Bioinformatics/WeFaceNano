#!/usr/bin/env bash

BARCODING=0
MLST=0
RES=0
SPECIES=0

while getopts bmrsh:-c:-i:-o:-f: option
do
case "${option}" in
    b|barcoding)
    BARCODING=1
    ;;
    c|config)
    CONFIG=${OPTARG}
    ;;
    i|input)
    FAST5Files=${OPTARG}
    ;;
    o|output)
    FILEOUT=${OPTARG}
    ;;
    f|format)
    FORMAT=${OPTARG}
    ;;
    m|mlst)
    MLST=1
    ;;
    r|res)
    RES=1
    ;;
    s|species)
    SPECIES=1
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
echo "$BARCODING"
if [ "${BARCODING}" = 1 ]; then
      gnome-terminal -e 'bash -c "read_fast5_basecaller.py --barcoding -r -i '"${FAST5Files}"' -c '"${CONFIG}"' -t $(nproc) -o '"${FORMAT}"' -s '"${FILEOUT}"' 2> /dev/null"' &
      sleep 5
      while (ps -a | grep -c read_fast5); do
        echo "Albacore is running..."
      done
      FINISHED=${FILEOUT}/finished.dat
      touch ${FINISHED}
fi
if [ "${BARCODING}" = 0 ]; then
      gnome-terminal -e 'bash -c "read_fast5_basecaller.py -r -i '"${FAST5Files}"' -c '"${CONFIG}"' -t $(nproc) -o '"${FORMAT}"' -s '"${FILEOUT}"' 2> /dev/null"' &
      sleep 5
      while (ps -a | grep -c read_fast5); do
        echo "Albacore is running..."
      done
      FINISHED=${FILEOUT}/finished.dat
      touch ${FINISHED}
fi
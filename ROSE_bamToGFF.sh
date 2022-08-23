#!/usr/bash

CLEAR='\033[0m'
RED='\033[0;31m'

function usage(){
    if [ -n "$1" ]; then
        echo -e "${RED}→ $1${CLEAR}"
    fi
    echo "Usage: $0 [-h help] [-b bam] [-i input] [-o output] [-s sense] [-f floor] [-e extension] [-r rpm] [-m matrix] [-v verbose]"
    echo " #Required arguments"
    echo " -b, --bam        .bam file to process"
    echo " -i, --input      Enriched region .gff3 file"
    echo " -o, --output     Output file name"
    echo ""
    echo " #Additional arguments"
    echo " -s, --sense      Strand to map to (default='both')"
    echo " -f, --floor      Read floor threshold necessary to count towards density (default=1)"
    echo " -e, --extension  Extends reads by n bp (default=200)"
    echo " -r, --rpm        Normalizes density to reads per million (rpm)"
    echo " -m, --matrix     Variable bin sized matrix (default=1)"
    echo " -v, --verbose    Print verbose messages"
    echo ""
    echo "Example: $0 -b ./data/MM1S_MED1.hg18.bwt.sorted.bam -i example/gff/HG18_MM1S_MED1_12.5kb_stitched_TSS_distal.gff3 -o example/mappedGFF/HG18_MM1S_MED1_12.5kb_stitched_TSS_distal_MM1S_MED1.hg18.bwt.sorted.bam_mapped.gff3 -s both -e 200 -r True -m 1 -v true"
    exit 1
}

#Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) usage ;;
        -v|--verbose) VERBOSE="$2"; shift ;;
        -b|--bam) BAM="$2"; shift ;;
        -i|--input) INPUT="$2"; shift ;;
        -o|--output) OUTPUT="$2"; shift ;;
        -s|--sense) SENSE="$2"; shift ;;
        -f|--floor) FLOOR="$2"; shift ;;
        -e|--extension) EXTENSION="$2"; shift ;;
        -r|--rpm) RPM="$2"; shift ;;
        -m|--matrix) MATRIX="$2"; shift ;;
    esac
    shift
done

#Verifying arguments
if [ -z "$VERBOSE" ]; then VALUE_V=false; else VALUE_V=true; fi;
if [ -z "$BAM" ]; then usage ".bam file is not specified"; else VALUE_B=$BAM; fi;
if [ -z "$INPUT" ]; then usage "Input enriched region file is not specified"; else VALUE_I=$INPUT; fi;
if [ -z "$OUTPUT" ]; then usage "Output directory name is not specified"; else VALUE_O=$OUTPUT; fi;
if [ -z "$SENSE" ]; then VALUE_S="both"; else VALUE_S=$SENSE; fi;
if [ -z "$FLOOR" ]; then VALUE_F=1; else VALUE_F=$FLOOR; fi;
if [ -z "$EXTENSION" ]; then VALUE_E=200; else VALUE_E=$EXTENSION; fi;
if [ -z "$RPM" ]; then VALUE_R=true; else VALUE_R=$RPM; fi;
if [ -z "$MATRIX" ]; then VALUE_M=1; else VALUE_M=$MATRIX; fi;

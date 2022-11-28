#!/usr/bin/env bash

CLEAR='\033[0m'
RED='\033[0;31m'

function usage(){
    if [ -n "$1" ]; then
        echo -e "${RED}→ $1${CLEAR}"
    fi
    echo "Usage: $0 [-h help] [-g genome] [-i input] [-o output] [-r rankby] [optional flags]"
    echo " #Required arguments"
    echo " -g, --genome     Genome build (MM8, MM9, MM10, HG18, HG19, HG38)"
    echo " -i, --input      File (.bed, .gff or .gtf) containing enhancer binding sites"
    echo " -o, --output     Name of output directory where data will be stored"
    echo " -r, --rankby     .bam file to rank enhancers by"
    echo " -a, --annot      UCSC table track annotation file"
    echo ""
    echo " #Additional arguments"
    echo " -v, --verbose    Print verbose messages (default=true)"
    echo ""
    echo " #Additional arguments for ROSE_main.py"
    echo " -c, --control    .bam file to rank enhancers by"
    echo " -s, --stitch     Max linking distance for stitching (default=12500)"
    echo " -t, --tss        Distance from TSS to exclude (0 = no TSS exclusion) (default=0)"
    echo " -d, --debug      Enhancer stitching debugging output (default=False)"
    echo ""
    echo " #Additional arguments for ROSE_bamToGFF.py"
    echo " -n, --sense      Strand to map to (default='both')"
    echo " -f, --floor      Read floor threshold necessary to count towards density (default=1)"
    echo " -x, --extension  Extends reads by n bp (default=200)"
    echo " -p, --rpm        Normalizes density to reads per million (rpm) (default=true)"
    echo " -m, --matrix     Variable bin sized matrix (default=1)"
    echo ""
    echo "Example: $0 -g hg18 -i ./data/HG18_MM1S_MED1.gff -o output -r ./data/MM1S_MED1.hg18.bwt.sorted.bam -a ./data/annotation/hg18_refseq.ucsc -c ./data/MM1S_WCE.hg18.bwt.sorted.bam -s 12500 -t 2500 -n both -x 200 -p true -m 1 -v true"
    exit 1
}

#Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) usage ;;
        -v|--verbose) VERBOSE="$2"; shift ;;
        -g|--genome) GENOME="$2"; shift ;;
        -i|--input) INPUT="$2"; shift ;;
        -o|--output) OUTPUT="$2"; shift ;;
        -r|--rankby) RANKBY="$2"; shift ;;
        -a|--annot) ANNOT="$2"; shift ;;
        -c|--control) VALUE_C="$2"; shift ;;
        -s|--stitch) STITCH="$2"; shift ;;
        -t|--tss) TSS="$2"; shift ;;
        -d|--debug) DEBUG="$2"; shift ;;
        -n|--sense) SENSE="$2"; shift ;;
        -f|--floor) FLOOR="$2"; shift ;;
        -x|--extension) EXTENSION="$2"; shift ;;
        -p|--rpm) RPM="$2"; shift ;;
        -m|--matrix) MATRIX="$2"; shift ;;
    esac
    shift
done

#Verifying arguments
if [ -z "$VERBOSE" ]; then VALUE_V=false; else VALUE_V=true; fi;
if [ -z "$GENOME" ]; then usage "Genome build is not specified"; else VALUE_G=$GENOME; fi;
if [ -z "$INPUT" ]; then usage "Input file is not specified"; else VALUE_I=$INPUT; fi;
if [ -z "$OUTPUT" ]; then usage "Output directory name is not specified"; else VALUE_O=$OUTPUT; fi;
if [ -z "$RANKBY" ]; then usage ".bam file is not specified"; else VALUE_R=$RANKBY; fi;
if [ -z "$ANNOT" ]; then usage "UCSC annotation file is not specified"; else VALUE_A=$ANNOT; fi;
if [ -z "$STITCH" ]; then VALUE_S=12500; else VALUE_S=$STITCH; fi;
if [ -z "$TSS" ]; then VALUE_T=0; else VALUE_T=$TSS; fi;
if [ -z "$DEBUG" ]; then VALUE_D=false; else VALUE_D=true; fi;
if [ -z "$SENSE" ]; then VALUE_N="both"; else VALUE_N=$SENSE; fi;
if [ -z "$FLOOR" ]; then VALUE_F=1; else VALUE_F=$FLOOR; fi;
if [ -z "$EXTENSION" ]; then VALUE_X=200; else VALUE_X=$EXTENSION; fi;
if [ -z "$RPM" ]; then VALUE_P=true; else VALUE_P=$RPM; fi;
if [ -z "$MATRIX" ]; then VALUE_M=1; else VALUE_M=$MATRIX; fi;


#Initialising variables
BAM_FILES=($VALUE_R)
if [ "$VALUE_C" ]; then
    BAM_FILES+=($VALUE_C)
fi

#Create stitched enhancers .gff3 file
python3 src/ROSE_main.py -g $VALUE_G -i $VALUE_I -r $VALUE_R -a $VALUE_A -o $VALUE_O -c $VALUE_C -s $VALUE_S -t $VALUE_T -d $VALUE_D -v $VALUE_V

#Creating variable names
ORIGINAL=$(find ${PWD}/${VALUE_O}/gff/ -name "$(basename ${VALUE_I} | cut -d "." -f 1).gff3")
STITCHED=$(find ${PWD}/${VALUE_O}/gff/ -name "$(basename ${VALUE_I} | cut -d "." -f 1)*_distal.gff3" | head -n 1)

#Mapping reads to stitched enhancers gff
for BAM in ${BAM_FILES[@]}; do
    python3 src/ROSE_bamToGFF.py -b $BAM -i $STITCHED -s $VALUE_N -f $VALUE_F -e $VALUE_X -r $VALUE_P -m $VALUE_M -v $VALUE_V &
    python3 src/ROSE_bamToGFF.py -b $BAM -i $ORIGINAL -s $VALUE_N -f $VALUE_F -e $VALUE_X -r $VALUE_P -m $VALUE_M -v $VALUE_V &
    wait
done

#Calculate read density for each stitched enhancer locus
python3 src/ROSE_mapCollection.py -s $STITCHED -g $VALUE_I -b "${BAM_FILES[@]}" -d ${PWD}/${VALUE_O}/mappedGFF

#Creating more variable names
DENSITY=$(find ${PWD}/${VALUE_O}/ -name "$(basename ${VALUE_I} | cut -d "." -f 1)*_enhancer_region_map.txt")

#Identifing and visualising superenhancers
if [ "$VALUE_C" ]; then
    Rscript src/ROSE_callSuper.R -o ${PWD}/${VALUE_O} -d $DENSITY -g $VALUE_I -c $VALUE_C
else
    Rscript src/ROSE_callSuper.R -o ${PWD}/${VALUE_O} -d $DENSITY -g $VALUE_I
fi
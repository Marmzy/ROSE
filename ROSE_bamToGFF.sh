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

#Verifying arguments' content
if [ -z "$(find $(dirname ${VALUE_B}) -name $(basename ${VALUE_B}).bai)" ]; then usage "No associated .bai file found with .bam file"; fi;
if [ -z $(echo "+", "-", ".", "both" | grep -wo $VALUE_S) ]; then usage "Sense flag argument must be '+', '-', '.' or 'both'"; fi;


#Get total number of mapped alignments
if $VALUE_R; then
    MMR=$(samtools flagstat $VALUE_B | grep "mapped (" | cut -d" " -f 1)
    # MMR=17414095
else
    MMR=1
fi

if $VALUE_V; then
    echo "Using a MMR value of $MMR"
fi

#Add genome file option for parse args???
#Preparing output file names
EXTENDED_GFF=$(dirname $(dirname ${VALUE_I}))/mappedGFF/$(basename ${VALUE_I::-5})_extended.gff3
EXTENDED_BED=$(dirname $(dirname ${VALUE_I}))/mappedGFF/$(basename ${VALUE_I::-5})_extended.bed
# EXTENDED_SAM=$(dirname $(dirname ${VALUE_I}))/mappedGFF/$(basename ${VALUE_I::-5})_extended.sam

#Extend stitched enhancers regions gff file and convert it to .bed format
bedtools slop -i $VALUE_I -g "/home/calvin/Projects/ROSE/data/human.hg18.genome" -b $VALUE_E > $EXTENDED_GFF
gff2bed < $EXTENDED_GFF > $EXTENDED_BED

#Find reads mapped in stitched enhancer loci
mkdir -p $(dirname $(dirname ${VALUE_I}))/mappedGFF/regions
while read line; do
    REGION=$(echo $line | awk '{print $1":"$2+1"-"$3}')
    samtools view $VALUE_B $REGION | awk '$6 !~ /N/ {print}' > $(dirname $(dirname ${VALUE_I}))/mappedGFF/regions/${REGION}.sam
done < $EXTENDED_BED

#????
python3 ROSE_bamToGFF.py -b $VALUE_B -i $VALUE_I -r $(dirname $(dirname ${VALUE_I}))/mappedGFF/regions -s $VALUE_S -f $VALUE_F -e $VALUE_E -m $MMR -x $VALUE_M -v







#Find reads mapped in stitched enhancer loci
# samtools view -L $EXTENDED_BED $VALUE_B | awk '$6 !~ /N/ {print}' > $EXTENDED_SAM                            <- replace with loop and output each match to new file?
#Split output file per chromosome
# mkdir -p $(dirname $(dirname ${VALUE_I}))/mappedGFF/chromosomes
# CHROMS=$(samtools idxstats ${VALUE_B} | cut -f1 | grep -v '*')
# for chr in ${CHROMS[@]}; do
#     CHR_FILE=$(dirname ${EXTENDED_SAM})/chromosomes/$(basename ${EXTENDED_SAM::-4}_${chr}.sam)
#     awk -v chr="$chr" -v out="$CHR_FILE" '$3==chr {print > out}' $EXTENDED_SAM
# done
# python3 ROSE_bamToGFF.py -i $VALUE_I -c $(dirname $(dirname ${VALUE_I}))/mappedGFF/chromosomes -s $VALUE_S -f $VALUE_F -e $VALUE_E -m $MMR -x $VALUE_M -v
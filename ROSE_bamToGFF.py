#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd

from collections import defaultdict
from pathlib import Path
from src.utils.annotation import get_strand
from src.utils.file_helper import check_file, check_path
from src.utils.locus import Locus


def str2bool(
    v: str
) -> bool:
    """Convert string to boolean

    Args:
        v (str): boolean string

    Raises:
        argparse.ArgumentTypeError: String is not named "true" or "false"

    Returns:
        bool: Booleanised string
    """
    
    if v.lower() == "true":
        return True
    elif v.lower() == "false":
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def parseArgs() -> argparse.Namespace:
    """Parse arguments from CLI

    Returns:
        argparse.Namespace: Argparse space containing parsed arguments
    """

    parser = argparse.ArgumentParser(description="Make locus objects from reads mapped to stitched enhancer loci")

    #Required arguments
    parser.add_argument("-b", "--bam", type=str, help=".bam file")
    parser.add_argument("-i", "--input", type=str, help="Stitched enhancer loci .gff3 file")
    parser.add_argument("-r", "--region",  type=str, help="Regions directory path")
    parser.add_argument("-m", "--mmr", type=int, help="Total number of mapped reads")

    #Optional arguments
    parser.add_argument("-e", "--extension", type=int, nargs="?", default=200, help="Extend reads by n bp (default: 200)")
    parser.add_argument("-f", "--floor", type=int, nargs="?", default=1, help="Read floor threshold necessary to count towards density (default: 1)")
    parser.add_argument("-s", "--sense", type=str, nargs="?", default="both", help="Strand to map to (default: 'both')")
    parser.add_argument("-x", "--matrix", type=int, nargs="?", default=1, help="Variable bin sized matrix (default: 1)")
    parser.add_argument("-v", "--verbose", type=str2bool, nargs="?", const=True, default=False, help="Print verbose messages")


    #Printing arguments to the command line
    args = parser.parse_args()

    print("Called with args:")
    print(f"{args}\n")

    #Ensuring that argument files exist
    check_file(args.input)

    return args
        

def main() -> None:

    #Parse arguments from the command line
    args = parseArgs()

    #Initialising variables
    extendedReads = []
    newGFF = []

    if args.mmr != 1:
        mmr = round(args.mmr/1000000, 4)
    else:
        mmr = args.mmr

    #Reading the gff3 file as a dataframe
    gff_df = pd.read_csv(args.input, sep="\t", header=None, comment="#")
    
    #Loop over stitched enhancer loci
    for row in zip(*gff_df.to_dict("list").values()):

        #Create locus from stitched enhancer locus
        gffLocus = Locus(row[0], row[3], row[4], row[6], row[8])

        #Get reads that mapped to the stitched anhancer locus
        sam_file = str(Path(args.region, f"{row[0]}:{row[3]-args.extension}-{row[4]+args.extension}.sam"))
        mini_df = pd.read_csv(sam_file, sep="\t", header=None)

        #Loop over mapped reads
        for read in zip(*mini_df.to_dict("list").values()):
            if get_strand(read[1]) == "+":
                #Read length does not account for gaps -> add in later update
                locus = Locus(read[2], read[3], read[3]+len(read[9])+args.extension, get_strand(read[1]))
            else:
                locus = Locus(read[2], read[3]-args.extension, read[3]+len(read[9]), get_strand(read[1]))
            extendedReads.append(locus)

        #Define sense and antisense reads
        if gffLocus._sense == "+" or gffLocus._sense == ".":
            senseReads = [read for read in extendedReads if read._sense == "+"]
            antiReads = [read for read in extendedReads if read._sense == "-"]
        else:
            senseReads = [read for read in extendedReads if read._sense == "-"]
            antiReads = [read for read in extendedReads if read._sense == "+"]

        #Create dictionary of number of reads mapped at each genomic position 
        if args.sense == "+" or args.sense == "both" or args.sense == ".":
            sense_unique, sense_counts = np.unique(np.concatenate([np.arange(read._start, read._end+1) for read in senseReads]), return_counts=True)
            senseHash = defaultdict(int, zip(sense_unique, sense_counts))
        if args.sense == "-" or args.sense == "both" or args.sense == ".":
            anti_unique, anti_counts = np.unique(np.concatenate([np.arange(read._start, read._end+1) for read in antiReads]), return_counts=True)
            antiHash = defaultdict(int, zip(anti_unique, anti_counts))

        #Remove positions in hash with less than or equal 'floor' reads mapped
        #and positions outside the stitched enhancer locus
        keys = [k for k in set(list(senseHash.keys()) + list(antiHash.keys())) if senseHash[k]+antiHash[k] > args.floor if gffLocus._start < k < gffLocus._end]

        #Creating bin sizes for calculating read density in stitched enhancer locus
        binSize = (len(gffLocus)-1) / int(args.matrix)
        nBins = int(args.matrix)

        # #Unnecessary check?
        # if not binSize:
        #     clusterLine += ["NA"]*int(args.matrix)
        #     newGFF.append(clusterLine)

        clusterLine = [gffLocus._ID, str(gffLocus)]

        #Calculate stitched enhancer locus mapped read density per bin
        n = 0
        if gffLocus._sense == "+" or gffLocus._sense == "." or gffLocus._sense == "both":
            i = gffLocus._start
            while n < nBins:
                n += 1
                binKeys = [k for k in keys if i < k < i+binSize]
                binDen = float(sum([senseHash[x]+antiHash[x] for x in binKeys])) / binSize
                clusterLine += [round(binDen / mmr, 4)]
                i = i + binSize
        else:
            i = gffLocus._end
            while n < nBins:
                n += 1
                binKeys = [k for k in keys if i-binSize < k < i]
                binDen = float(sum([senseHash[x]+antiHash[x] for x in binKeys])) / binSize
                clusterLine += [round(binDen / mmr, 4)]
                i = i - binSize
        newGFF.append(clusterLine)

    #
    out_df = pd.DataFrame(newGFF, columns=["GENE_ID", "locusLine"] + [f"bin_{n}_{str(Path(args.bam).name)}" for n in range(1, int(args.matrix)+1)])
    gff_name = check_path(Path(Path(args.input).parents[1], "mappedGFF", f"{Path(args.input).stem}_{Path(args.bam).name}_mapped.txt"))
    out_df.to_csv(gff_name, sep="\t", index=False)

                
if __name__ == "__main__":
    main()
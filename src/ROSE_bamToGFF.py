#!/usr/bin/env python

import argparse
import glob
import multiprocessing as mp
import pandas as pd
import re

from classes.bam import Bam
from classes.locus import Locus
from collections import Counter
from itertools import product
from pathlib import Path
from typing import List, Union
from utils.file_helper import check_file, check_path, str2bool


def parseArgs() -> argparse.Namespace:
    """Parse arguments from CLI

    Returns:
        argparse.Namespace: Argparse space containing parsed arguments
    """

    parser = argparse.ArgumentParser(description="Make locus objects from reads mapped to stitched enhancer loci")

    #Required arguments
    parser.add_argument("-b", "--bam", type=str, help="Directory containing .bam files")
    parser.add_argument("-c", "--control", type=str, help="Control .bam file")
    parser.add_argument("-i", "--input", type=str, help="Stitched enhancer loci .gff3 file")
    parser.add_argument("-o", "--original", type=str, help="Enhancer loci .gff3 file")

    #Optional arguments
    parser.add_argument("-s", "--sense", type=str, nargs="?", default="both", help="Strand to map to (default: 'both')")
    parser.add_argument("-f", "--floor", type=int, nargs="?", default=1, help="Read floor threshold necessary to count towards density (default: 1)")
    parser.add_argument("-e", "--extension", type=int, nargs="?", default=200, help="Extend reads by n bp (default: 200)")
    parser.add_argument("-r", "--rpm", type=str2bool, nargs="?", help="Normalize density to reads per million")
    parser.add_argument("-m", "--matrix", type=int, nargs="?", default=1, help="Variable bin sized matrix (default: 1)")
    parser.add_argument("-v", "--verbose", type=str2bool, nargs="?", const=True, default=False, help="Print verbose messages")


    #Printing arguments to the command line
    args = parser.parse_args()

    print(f"Called with args:\n{args}\n")

    #Ensuring that files exist
    check_file(args.control)
    check_file(args.input)
    check_file(args.original)

    #Ensuring sense argument makes sense
    if args.sense not in ["+", "-", ".", "both"]:
        raise ValueError("Argument -s/--sense value must be '+', '-', '.' or 'both'")

    return args


def map_reads(
    bam_file: str,
    gff_df: pd.core.frame.DataFrame,
    rpm: bool,
    extension: int,
    sense: str,
    floor: int,
    matrix: int,
    verbose: bool
) -> List[List[Union[str, float]]]:
    """Map reads to stitched enhancer loci and calculate read density

    Args:
        bam_file (str): .bam file whose reads are to be mapped
        gff_df (pd.core.frame.DataFrame): Stitched enhancer loci dataframe
        rpm (bool): Boolean whether to normalize read density
        extension (int): Number of bp to extend reads by
        sense (str): Strand to map to
        floor (int): Threshold of mapped reads to count towards read density
        matrix (int): Number of bins per stitched enhancer locus to count density for
        verbose (bool): Detailed output boolean

    Returns:
        List[List[Union[str, float]]]: List of lists containing density of reads for each stitched enhancer locus
    """

    #Initialising variables
    newGFF = []

    #Create Bam class
    bam = Bam(bam_file)
    mmr = round(bam.getTotalReads()/1000000, 4) if rpm else 1

    if verbose:
        print(f"MMR value: {mmr}")

    #Check chromosome naming convention
    bam.checkChrStatus()

    #Loop over stitched enhancer loci
    for row in zip(*gff_df.to_dict("list").values()):

        #Create locus object from stitched enhancer locus
        if not bam._chr:
            row[0] = re.sub("chr", "", row[0])
        gffLocus = Locus(row[0], row[3], row[4], row[6], row[8])

        #Get reads that lie in the extended stitched enhancer locus region
        searchLocus = Locus(gffLocus._chr, gffLocus._start-extension, gffLocus._end+extension, gffLocus._sense, gffLocus._ID)
        reads = bam.getReadsLocus(searchLocus)

        #Extend reads
        extendedReads = []
        for locus in reads:
            if locus._sense == "+":
                locus = Locus(locus._chr, locus._start, locus._end+extension, locus._sense, locus._ID)
            if locus._sense == "-":
                locus = Locus(locus._chr, locus._start-extension, locus._end, locus._sense, locus._ID)
            extendedReads.append(locus)

        #Define sense and antisense reads
        if gffLocus._sense == "+" or gffLocus._sense == ".":
            senseReads = [read for read in extendedReads if read._sense == "+"]
            antiReads = [read for read in extendedReads if read._sense == "-"]
        else:
            senseReads = [read for read in extendedReads if read._sense == "-"]
            antiReads = [read for read in extendedReads if read._sense == "+"]

        #Create dictionary of number of reads mapped at each genomic position 
        if sense == "+" or sense == "both" or sense == ".":
            senseHash = Counter([i for read in senseReads for i in range(read._start, read._end+1)])
        if sense == "-" or sense == "both" or sense == ".":
            antiHash = Counter([i for read in antiReads for i in range(read._start, read._end+1)])

        #Remove positions in hash with less than or equal 'floor' reads mapped
        #and positions outside the stitched enhancer locus
        keys = [k for k in set(list(senseHash.keys()) + list(antiHash.keys())) if senseHash[k]+antiHash[k] > floor if gffLocus._start < k < gffLocus._end]

        #Creating bin sizes for calculating read density in stitched enhancer loci
        binSize = (len(gffLocus)-1) / int(matrix)
        nBins = int(matrix)

        clusterLine = [gffLocus._ID, str(gffLocus)]

        #Calculate density of mapped reads per bin in the stitched enhancer locus
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

    return newGFF
        

def main() -> None:
    """Calculate read density per stitched enhancer locus bin
    """

    #Parse arguments from the command line
    args = parseArgs()

    #Reading the gff3 file as a dataframe
    gff_df = pd.read_csv(check_file(args.input), sep="\t", header=None, comment="#")
    original_df = pd.read_csv(check_file(args.original), sep="\t", header=None, comment="#")

    #Loading the .bam files from the directory and making sure they are indexed
    bam_files = glob.glob(str(Path(args.bam, "*.bam")))
    bam_files.append(args.control)

    for bam in bam_files:
        check_file(f"{bam}.bai")

    #Map reads to stitched enhancer loci and calculate read density
    with mp.Pool(mp.cpu_count()) as p:
        results = p.starmap(map_reads, [(bam, gff, args.rpm, args.extension, args.sense, args.floor, args.matrix, args.verbose) for bam, gff in product(bam_files, [gff_df, original_df])])
        
    #Outputting per-bin read density
    for newGFF, (bam, gff) in zip(results, product(bam_files, [args.input, args.original])):
        out_df = pd.DataFrame(newGFF, columns=["GENE_ID", "locusLine"] + [f"bin_{n}_{str(Path(bam).name)}" for n in range(1, int(args.matrix)+1)])
        out_name = check_path(Path(Path(gff).parents[1], "mappedGFF", f"{Path(gff).stem}_{Path(bam).stem}_mapped.txt"))
        out_df.to_csv(out_name, sep="\t", index=False)

                
if __name__ == "__main__":
    main()
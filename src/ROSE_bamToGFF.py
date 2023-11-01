#!/usr/bin/env python

import multiprocessing as mp
import pandas as pd
import re

from collections import Counter
from itertools import product
from pathlib import Path
from src.classes.bam import Bam
from src.classes.locus import Locus
from src.utils.file_helper import check_file, check_path
from typing import List, Union


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
        List[List[Union[str, float]]]: List of lists containing density of
                                       reads for each stitched enhancer locus
    """

    # Initialising variables
    newGFF = []

    # Create Bam class
    bam = Bam(bam_file)
    mmr = round(bam.getTotalReads()/1000000, 4) if rpm else 1

    if verbose:
        print(f"MMR value: {mmr}")

    # Check chromosome naming convention
    bam.checkChrStatus()

    # Loop over stitched enhancer loci
    for row in zip(*gff_df.to_dict("list").values()):

        # Create locus object from stitched enhancer locus
        if not bam._chr:
            row[0] = re.sub("chr", "", row[0])
        gffLocus = Locus(row[0], row[3], row[4], row[6], row[8])

        # Get reads that lie in the extended stitched enhancer locus region
        searchLocus = Locus(gffLocus._chr, gffLocus._start-extension, gffLocus._end+extension,
                            gffLocus._sense, gffLocus._ID)
        reads = bam.getReadsLocus(searchLocus)

        # Extend reads
        extendedReads = []
        for locus in reads:
            if locus._sense == "+":
                locus = Locus(locus._chr, locus._start, locus._end+extension,
                              locus._sense, locus._ID)
            if locus._sense == "-":
                locus = Locus(locus._chr, locus._start-extension, locus._end,
                              locus._sense, locus._ID)
            extendedReads.append(locus)

        # Define sense and antisense reads
        if gffLocus._sense == "+" or gffLocus._sense == ".":
            senseReads = [read for read in extendedReads if read._sense == "+"]
            antiReads = [read for read in extendedReads if read._sense == "-"]
        else:
            senseReads = [read for read in extendedReads if read._sense == "-"]
            antiReads = [read for read in extendedReads if read._sense == "+"]

        # Create dictionary of number of reads mapped at each genomic position
        if sense == "+" or sense == "both" or sense == ".":
            senseHash = Counter([i for read in senseReads for i in range(read._start, read._end+1)])
        if sense == "-" or sense == "both" or sense == ".":
            antiHash = Counter([i for read in antiReads for i in range(read._start, read._end+1)])

        # Remove positions in hash with less than or equal 'floor' reads mapped
        # and positions outside the stitched enhancer locus
        keys = [k for k in set(list(senseHash.keys()) + list(antiHash.keys()))
                if senseHash[k]+antiHash[k] > floor if gffLocus._start < k < gffLocus._end]

        # Creating bin sizes for calculating read density in
        # stitched enhancer loci
        binSize = (len(gffLocus)-1) / int(matrix)
        nBins = int(matrix)

        clusterLine = [gffLocus._ID, str(gffLocus)]

        # Calculate density of mapped reads per bin in
        # the stitched enhancer locus
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


def calc_read_density(
    bam_files: str,
    input: str,
    original: str,
    control: str = "",
    extension: int = 200,
    floor: int = 1,
    matrix: int = 1,
    rpm: bool = True,
    sense: str = "both",
    verbose: bool = True
) -> None:
    """Calculate read density per stitched enhancer locus bin

    Args:
        bam_files (str): List of target .bam file(s)
        control (str): Control .bam file
        input (str): Stiched enhancers .gff3 file
        original (str): Original .gff3 enhancers file
        extension (int): Number of bp to extend reads by
        floor (int): Threshold of mapped reads to count towards read density
        matrix (int): Number of bins per stitched enhancer locus to count density for
        rpm (bool): Boolean whether to normalize read density
        sense (str): Strand to map to
        verbose (bool): Detailed output boolean
    """

    # Reading the gff3 file as a dataframe
    gff_df = pd.read_csv(input, sep="\t", header=None, comment="#")
    original_df = pd.read_csv(original, sep="\t", header=None, comment="#")

    # Loading the .bam files from the directory and
    # making sure they are indexed
    bam_files = bam_files.copy()
    if control:
        bam_files.append(control)

    for bam in bam_files:
        check_file(f"{bam}.bai")

    # Map reads to stitched enhancer loci and calculate read density
    with mp.Pool(mp.cpu_count()) as p:
        results = p.starmap(
            map_reads,
            [(bam, gff, rpm, extension, sense, floor, matrix, verbose)
             for bam, gff in product(bam_files, [gff_df, original_df])]
        )

    # Outputting per-bin read density
    for newGFF, (bam, gff) in zip(results, product(bam_files, [input, original])):
        bin_cols = [f"bin_{n}_{str(Path(bam).name)}" for n in range(1, int(matrix)+1)]
        out_df = pd.DataFrame(newGFF, columns=["GENE_ID", "locusLine"] + bin_cols)
        out_name = check_path(
            Path(Path(gff).parents[1], "mappedGFF", f"{Path(gff).stem}_{Path(bam).stem}_mapped.txt")
        )
        out_df.to_csv(out_name, sep="\t", index=False)

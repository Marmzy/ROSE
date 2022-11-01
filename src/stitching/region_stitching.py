#!/usr/bin/env python

import pandas as pd

from classes.locus import LocusCollection, makeTSSLocus
from typing import List, Tuple


def makeStartDict(
    annotFile: str
) -> pd.core.frame.DataFrame:
    """Create subsetted dataframe from UCSC refseq file

    Args:
        annotFile (str): Path to UCSC refseq annotation file

    Returns:
        pd.core.frame.DataFrame: Subsetted annotation dataframe
    """

    #Reading the annotation file in as a dataframe
    refseqTable = pd.read_csv(annotFile, sep="\t")

    #Remove duplicate "name" entries and extract specific data
    startDict = refseqTable[refseqTable.columns.intersection(["name", "strand", "chrom", "txStart", "txEnd", "name2"])].copy()
    startDict.drop_duplicates(subset=["name"], inplace=True)
    startDict.loc[startDict["strand"]=="-", ["txStart", "txEnd"]] = (startDict.loc[startDict["strand"] == "-", ("txEnd", "txStart")].values)
    startDict.rename({"name": "id", "strand": "sense", "chrom": "chr", "txStart": "start", "txEnd": "end", "name2": "name"}, axis=1, inplace=True)

    return startDict


def regionStitching(
    boundCollection: LocusCollection,
    stitchWindow: int,
    tssWindow: int,
    startDict: pd.core.frame.DataFrame,
    removeTSS: bool
) -> Tuple[LocusCollection, List[List[str]]]:
    """Stitch enhancer loci that are within a certain distance together

    Args:
        boundCollection (LocusCollection): Collection of enhancer loci derived from the input file
        stitchWindow (int): Distance to stitch enhancer loci together
        tssWindow (int): Max distance for stitched loci to removed from a TSS loci
        startDict (pd.core.frame.DataFrame): UCSC refseq dataframe
        removeTSS (bool): Bool to exclude stitched enhancer loci too close to TSS loci

    Returns:
        Tuple[LocusCollection, List[List[str]]]: Tuple of the stiched enhancer loci collection and list stiched loci that overlap multiple TSS
    """

    #Initialising variables
    debugOutput = []

    #Filter regions that overlap the TSS of an active gene
    if removeTSS:
        
        #Initialising variables
        removeTicker=0

        #Create loci centered around +/- tssWindow of transcribed genes
        tssLoci = [makeTSSLocus(row, tssWindow) for row in zip(*startDict.to_dict("list").values())]
        tssCollection = LocusCollection(tssLoci, 50)  ## -> isn't 500 the standard?? 

        #Get all enhancer loci
        boundLoci = boundCollection.getLoci()

        #Loop over enhancer loci and remove those that are enveloped within an active gene's
        #TSS +/- tssWindow region 
        for locus in list(boundLoci):
            if len(tssCollection.getContainers(locus, "both")) > 0:
                boundCollection.remove(locus)
                debugOutput.append([str(locus), locus._ID, "Contained"])
                removeTicker += 1

        print(f"Removed {removeTicker} loci because they contain a TSS")

    #Create a LocusCollection of stitched
    stitchedCollection = boundCollection.stitchCollection(stitchWindow, "both")

    #Remove stitched enhancer loci that overlap with 2+ TSS
    if removeTSS:

        #Initialising variables
        fixedLoci = []
        originalTicker = 0
        removeTicker = 0

        #Create loci centered around of transcribed genes
        tssLoci = [makeTSSLocus(row, 50) for row in zip(*startDict.to_dict("list").values())]
        tssCollection = LocusCollection(tssLoci, 50)  ## -> isn't 500 the standard?? 

        #Loop over stiched enhancer loci
        for stitchedLocus in stitchedCollection.getLoci():

            #Get names of transcribed genes loci that overlap with the stiched locus
            overlappingTSSLoci = tssCollection.getOverlap(stitchedLocus, "both")
            tssNames = startDict.loc[startDict["id"].isin([str(tssLocus._ID) for tssLocus in overlappingTSSLoci]), "name"].values
            tssNames = set(tssNames)
    
            #Remove stiched enhancer loci that overlap with 2+ gene loci and unstitch them
            if len(tssNames) > 2:
                originalLoci = boundCollection.getOverlap(stitchedLocus, "both")
                originalTicker += len(originalLoci)
                fixedLoci += originalLoci
                debugOutput.append([str(stitchedLocus), stitchedLocus._ID, "Multiple_TSS"])
                removeTicker+=1
            else:
                fixedLoci.append(stitchedLocus)

        print(f"Removed {removeTicker} stitched loci because they overlapped with multiple TSSs")
        print(f"Added back {originalTicker} original enhancer loci\n")
        fixedCollection = LocusCollection(fixedLoci, 50)
        return fixedCollection, debugOutput
    else:
        return stitchedCollection, debugOutput
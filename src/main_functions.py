#!/usr/bin/env python

import pandas as pd

from src.utils.locus import LocusCollection, gffToLocusCollection, makeTSSLocus
from typing import List, Tuple


# def checkRefCollection(
#     referenceCollection: LocusCollection
# ) -> None:
#     """Checking that loci in the collection have unique names

#     Args:
#         referenceCollection (LocusCollection): Collection of all .gff3 entries as loci

#     Raises:
#         ValueError: If non-unique identifiers are found among .gff3 entries
#     """

#     #Getting loci IDs
#     namesList = [locus._ID for locus in referenceCollection.getLoci()]
    
#     #Ensuring loci all have unique identifiers
#     if len(namesList) != len(set(namesList)):
#         raise ValueError("Last column (attributes column) of the .gff3 file must contain unique identifiers for each entry")
#     else:
#         print("Reference collection passes qc")


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

    #filter out all bound regions that overlap the TSS of an ACTIVE GENE
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
                removeTicker+=1

        print(f"Removed {removeTicker} loci because they contain a TSS")

    #Create a LocusCollection of stitched
    stitchedCollection = boundCollection.stitchCollection(stitchWindow, "both")

    if removeTSS:
        #now replace any stitched region that overlap 2 distinct genes
        #with the original loci that were there

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
        print(f"Added back {originalTicker} original enhancer loci")
        fixedCollection = LocusCollection(fixedLoci, 50)
        return fixedCollection, debugOutput
    else:
        return stitchedCollection, debugOutput


def mapCollection(stitchedCollection,referenceCollection,bamFileList,mappedFolder,output,refName):


    '''
    makes a table of factor density in a stitched locus and ranks table by number of loci stitched together
    '''

    
    print('FORMATTING TABLE')
    loci = stitchedCollection.getLoci()

    locusTable = [['REGION_ID','CHROM','START','STOP','NUM_LOCI','CONSTITUENT_SIZE']]

    lociLenList = []

    #strip out any that are in chrY
    for locus in list(loci):
        if locus.chr() == 'chrY':
            loci.remove(locus)
    
    for locus in loci:
        #numLociList.append(int(stitchLocus.ID().split('_')[1]))
        lociLenList.append(locus.len())
        #numOrder = order(numLociList,decreasing=True)
    lenOrder = ROSE_utils.order(lociLenList,decreasing=True)
    ticker = 0
    for i in lenOrder:
        ticker+=1
        if ticker%1000 ==0:
            print(ticker)
        locus = loci[i]

        #First get the size of the enriched regions within the stitched locus
        refEnrichSize = 0
        refOverlappingLoci = referenceCollection.getOverlap(locus,'both')
        for refLocus in refOverlappingLoci:
            refEnrichSize+=refLocus.len()

        try:
            stitchCount = int(locus.ID().split('_')[0])
        except ValueError:
            stitchCount = 1
        
        locusTable.append([locus.ID(),locus.chr(),locus.start(),locus.end(),stitchCount,refEnrichSize])
        
            

    print('GETTING MAPPED DATA')
    for bamFile in bamFileList:
        
        bamFileName = bamFile.split('/')[-1]

        print('GETTING MAPPING DATA FOR  %s' % bamFile)
        #assumes standard convention for naming enriched region gffs
        
        #opening up the mapped GFF
        print('OPENING %s%s_%s_MAPPED.gff' % (mappedFolder,refName,bamFileName))

        mappedGFF =ROSE_utils.parseTable('%s%s_%s_MAPPED.gff' % (mappedFolder,refName,bamFileName),'\t')        

        signalDict = defaultdict(float)
        print('MAKING SIGNAL DICT FOR %s' % (bamFile))
        mappedLoci = []
        for line in mappedGFF[1:]:

            chrom = line[1].split('(')[0]
            start = int(line[1].split(':')[-1].split('-')[0])
            end = int(line[1].split(':')[-1].split('-')[1])
            mappedLoci.append(ROSE_utils.Locus(chrom,start,end,'.',line[0]))
            try:
                signalDict[line[0]] = float(line[2])*(abs(end-start))
            except ValueError:
                print('WARNING NO SIGNAL FOR LINE:')
                print(line)
                continue
                
                
        
        mappedCollection = ROSE_utils.LocusCollection(mappedLoci,500)
        locusTable[0].append(bamFileName)

        for i in range(1,len(locusTable)):
            signal=0.0
            line = locusTable[i]
            lineLocus = ROSE_utils.Locus(line[1],line[2],line[3],'.')
            overlappingRegions = mappedCollection.getOverlap(lineLocus,sense='both')
            for region in overlappingRegions:
                signal+= signalDict[region.ID()]
            locusTable[i].append(signal)

    ROSE_utils.unParseTable(locusTable,output,'\t')

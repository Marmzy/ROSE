#!/usr/bin/env python

import argparse
import pandas as pd

from collections import defaultdict
from pathlib import Path
from utils.file_helper import check_file, check_path
from classes.locus import gffToLocusCollection, Locus, LocusCollection


def parseArgs() -> argparse.Namespace:
    """Parse arguments from CLI

    Returns:
        argparse.Namespace: Argparse space containing parsed arguments
    """

    parser = argparse.ArgumentParser(description="Calculate read density signals for each stitched enhancer locus")

    #Required arguments
    parser.add_argument("-s", "--stitch", type=str, help="Stitched enhancer loci .gff3 file")
    parser.add_argument("-g", "--gff", type=str, help="File (.bed, .gff or .gtf) containing binding sites to make enhancers")
    parser.add_argument("-b", "--bams", type=str, nargs="+", help="List of .bam file used in analysis")
    parser.add_argument("-d", "--dir",  type=str, help="Directory where output will be stored")

    #Printing arguments to the command line
    args = parser.parse_args()

    print("Called with args:")
    print(f"{args}\n")

    return args


def mapCollection() -> None:
    """Calculate density signal for each stitched enhancer locus
    """

    #Parse arguments from the command line
    args = parseArgs()

    #Initialising variables
    locusTable = [["REGION_ID", "CHROM", "START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE"]]
    
    #Read binding sites and stitched enhancer loci files as LocusCollection object
    referenceCollection = gffToLocusCollection(check_file(args.gff))
    stitchedCollection = gffToLocusCollection(check_file(args.stitch))
    loci = list(stitchedCollection.getLoci())

    #Remove stitched enhancer loci located on chromosome 'Y'
    for locus in list(loci):
        if locus._chr == "chrY":
            loci.remove(locus)

    #Order stitched enhancer loci by their length
    loci.sort(reverse=True, key=lambda l: (len(l) is None, len(l)))
    
    #Get the size of the enriched regions within the stitched enhancer locus
    for locus in loci:
        refEnrichSize = sum(len(refLocus) for refLocus in list(referenceCollection.getOverlap(locus, "both")))

        try:
            stitchCount = int(locus._ID.split("_")[0])
        except ValueError:
            stitchCount = 1

        locusTable.append([locus._ID, locus._chr, locus._start, locus._end, stitchCount, refEnrichSize])

    #Calculate stitched enhancer loci signal density for each .bam file
    for bam in args.bams:

        #Initialising variables
        mappedLoci = []

        #Open mapped stitched enhancer loci file as a dataframe
        mappedGFF = check_file(Path(args.dir, f"{Path(args.stitch).stem}_{Path(bam).name}_mapped.txt"))
        mappedGFF = pd.read_csv(mappedGFF, sep="\t", header=0, comment="#")

        #Calculate signal for all stitched enhancer loci
        signalDict = defaultdict(float)
        for line in zip(*mappedGFF.to_dict("list").values()):
            chrom = line[1].split("(")[0]
            start = int(line[1].split(":")[-1].split("-")[0])
            end = int(line[1].split(":")[-1].split("-")[1])
            mappedLoci.append(Locus(chrom, start, end, ".", line[0]))
            signalDict[line[0]] = float(line[2])*(abs(end-start))
        
        mappedCollection = LocusCollection(mappedLoci, 500)
        locusTable[0].append(Path(bam).name)

        #Append signal data for table in desired order
        for i in range(1, len(locusTable)):
            signal = 0.0
            line = locusTable[i]
            lineLocus = Locus(line[1], line[2], line[3], ".")
            overlappingRegions = mappedCollection.getOverlap(lineLocus, sense="both")
            for region in overlappingRegions:
                signal += signalDict[region._ID]
            locusTable[i].append(signal)

    #Outputting stitched enhancer loci signal density values per .bam file
    out_df = pd.DataFrame(locusTable[1:], columns=locusTable[0])
    out_name = check_path(Path(Path(args.dir).parents[0], f"{Path(args.stitch).stem}_enhancer_region_map.txt"))
    out_df.to_csv(out_name, sep="\t", index=False)


if __name__ == "__main__":
    mapCollection()
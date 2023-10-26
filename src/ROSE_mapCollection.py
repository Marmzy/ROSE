#!/usr/bin/env python

import pandas as pd

from collections import defaultdict
from pathlib import Path
from src.classes.locus import gffToLocusCollection, Locus, LocusCollection
from src.utils.file_helper import check_file, check_path
from typing import Union


def map_collection(
    stitch: str,
    gff: str,
    bam_files: str,
    control: Union[str, None],
    output: str,
) -> str:
    """Calculate density signal for each stitched enhancer locus

    Args:
        stitch (str): Stiched enhancers .gff3 file
        gff (str): Original .gff3 enhancers file
        bam_files (str): List of target .bam file(s)
        control (str): Control .bam file
        output (str): Output directory

    Returns:
        str: Stitched enhancers signal density values file
    """

    # Initialising variables
    locusTable = [
        ["REGION_ID", "CHROM", "START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE"]
    ]

    # Loading the .bam files from the directory
    # and making sure they are indexed
    bam_files = bam_files.copy()
    if control:
        bam_files.append(control)

    # Read binding sites and stitched enhancer loci files as
    # LocusCollection objects
    referenceCollection = gffToLocusCollection(check_file(gff))
    stitchedCollection = gffToLocusCollection(check_file(stitch))
    loci = list(stitchedCollection.getLoci())

    # Remove stitched enhancer loci located on chromosome 'Y'
    for locus in list(loci):
        if locus._chr == "chrY":
            loci.remove(locus)

    # Order stitched enhancer loci by their length
    loci.sort(reverse=True, key=lambda x: (len(x) is None, len(x)))

    # Get the size of the enriched regions within the stitched enhancer locus
    for locus in loci:
        refEnrichSize = sum(
            len(refLocus)
            for refLocus in list(referenceCollection.getOverlap(locus, "both"))
        )

        try:
            stitchCount = int(locus._ID.split("_")[0])
        except ValueError:
            stitchCount = 1

        locusTable.append(
            [locus._ID, locus._chr, locus._start,
             locus._end, stitchCount, refEnrichSize]
        )

    # Calculate stitched enhancer loci signal density for each .bam file
    for bam in bam_files:

        # Initialising variables
        mappedLoci = []

        # Open mapped stitched enhancer loci file as a dataframe
        mappedGFF = check_file(
            Path(output, f"{Path(stitch).stem}_{Path(bam).stem}_mapped.txt")
        )
        mappedGFF = pd.read_csv(mappedGFF, sep="\t", header=0, comment="#")

        # Calculate signal for all stitched enhancer loci
        signalDict = defaultdict(float)
        for line in zip(*mappedGFF.to_dict("list").values()):
            chrom = line[1].split("(")[0]
            start = int(line[1].split(":")[-1].split("-")[0])
            end = int(line[1].split(":")[-1].split("-")[1])
            mappedLoci.append(Locus(chrom, start, end, ".", line[0]))
            signalDict[line[0]] = float(line[2])*(abs(end-start))

        mappedCollection = LocusCollection(mappedLoci, 500)
        locusTable[0].append(Path(bam).name)

        # Append signal data for table in desired order
        for i in range(1, len(locusTable)):
            signal = 0.0
            line = locusTable[i]
            lineLocus = Locus(line[1], line[2], line[3], ".")
            overlappingRegions = mappedCollection.getOverlap(
                lineLocus, sense="both"
            )
            for region in overlappingRegions:
                signal += signalDict[region._ID]
            locusTable[i].append(signal)

    # Outputting stitched enhancer loci signal density values per .bam file
    out_df = pd.DataFrame(locusTable[1:], columns=locusTable[0])
    out_name = check_path(
        Path(
            Path(output).parents[0],
            f"{Path(stitch).stem}_enhancer_region_map.txt"
        )
    )
    out_df.to_csv(out_name, sep="\t", index=False)

    return out_name

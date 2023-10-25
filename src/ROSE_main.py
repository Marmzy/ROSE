#!/usr/bin/env python

import pandas as pd

from pathlib import Path
from src.classes.locus import gffToLocusCollection, locusCollectionToGFF
from src.stitching.region_stitching import makeStartDict, regionStitching
from src.utils.conversion import bed_to_gff3, check_gff, gff_to_gff3, gtf_to_gff3
from src.utils.file_helper import get_path, check_path
from typing import Tuple


def stitch_loci(
    input: str,
    output: str,
    annot: str,
    stitch: int = 12500,
    tss: int = 0,
    debug: bool = False,
    verbose: bool = True,
) -> Tuple[str, str]:
    """Stitch enhancer loci together

    Args:
        input (str): Enhancer binding sites file
        output (str): Output directory name
        annot (str): UCSC annotation file
        stitch (int, optional): Max linking distance for stitching (default=12500)
        tss (int, optional): Distance from TSS to exclude (0 = no TSS exclusion) (default=0)
        debug (bool, optional): Enhancer stitching debugging output (default=False)
        verbose (bool, optional): Verbose messages (default=True)

    Raises:
        ValueError: if input file is not a .bed, .gtf, .gff or gff3 file

    Returns:
        Tuple[str, str]: Original .gff3 enhancers file \
                         and stiched enhancers .gff3 file
    """

    #Initialising variables
    path = get_path()
    inputGFFFile = check_path(Path(path, output, "gff", Path(input).stem)) + ".gff3"
    suffix = "_TSS_distal" if bool(tss) else ""

    #Copying/creating the input .gff3 file
    if verbose:
        print(f"Converting input {Path(input).suffix} file to .gff3 format")
    match Path(input).suffix:
        case ".bed":
            bed_to_gff3(input, inputGFFFile)

        case ".gff":
            gff_to_gff3(input, inputGFFFile)

        case ".gtf":
            gtf_to_gff3(input, inputGFFFile, full=False)

        case ".gff3":
            check_gff(input, inputGFFFile)
        
        case _:
            raise ValueError(
                "Input file must be a .bed, .gtf, .gff or .gff3 file"
            )
    
    #Using the .gff3 file to define enhancers
    if verbose:
        print(f"Using {inputGFFFile} as the input .gff3 file\n")
    inputName = str(Path(inputGFFFile).stem)

    #Making the start dict
    startDict = makeStartDict(annot)

    #Loading enhancers as loci collection object
    referenceCollection = gffToLocusCollection(inputGFFFile)

    #Stitching regions together
    if verbose:
        print("Stitching regions together")
    stitchedCollection, debugOutput = regionStitching(
        referenceCollection,
        stitch,
        tss,
        startDict,
        bool(tss)
    )

    #Create a gff3 dataframe from the stitched enhancers loci collection
    if verbose:
        print("Making GFF from stitched collection")
    stitchedGFF = locusCollectionToGFF(stitchedCollection)
    
    #Defining output file names
    stitchedGFFFile = check_path(
        Path(path, output, "gff", f"{inputName}_{stitch/1000}kb_stitched{suffix}.gff3")
    )
    debugOutFile = check_path(
        Path(path, output, "gff", f"{inputName}_{stitch/1000}kb_stitched{suffix}.debug")
    )

    #Outputting the gff3 dataframe
    with open(stitchedGFFFile, "w") as f_out:
        f_out.write("##gff-version 3\n##source-version ROSE\n")
        stitchedGFF.to_csv(f_out, sep="\t", header=False, index=False, mode="a")

    #Outputting the debugging information
    if debug:
        pd.DataFrame(
            debugOutput,
            columns=["Enhancer", "Region", "Reason"]
        ).to_csv(debugOutFile, index=False, sep="\t")

    return inputGFFFile, stitchedGFFFile
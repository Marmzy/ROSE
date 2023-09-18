#!/usr/bin/env python

import argparse
import pandas as pd

from classes.locus import gffToLocusCollection, locusCollectionToGFF
from pathlib import Path
from stitching.region_stitching import makeStartDict, regionStitching
from utils.conversion import bed_to_gff3, check_gff, gff_to_gff3, gtf_to_gff3
from utils.file_helper import get_path, check_file, check_path, str2bool


def parseArgs() -> argparse.Namespace:
    """Parse arguments from CLI

    Returns:
        argparse.Namespace: Argparse space containing parsed arguments
    """

    parser = argparse.ArgumentParser(description="Stitch regions together to form enhancers, map read density to stitched regions and \
                                                  rank enhancers by read density to discover super-enhancers")

    #Required arguments
    parser.add_argument("-i", "--input", type=str, help="File (.bed, .gff or .gtf) containing binding sites to make enhancers")
    parser.add_argument("-o", "--output", type=str, help="Output directory name")
    parser.add_argument("-a", "--annot", type=str, help="UCSC (.ucsc) annotation file")

    #Optional arguments
    parser.add_argument("-s", "--stitch", type=int, nargs="?", default=12500, help="Max linking distance for stitching")
    parser.add_argument("-t", "--tss", type=int, nargs="?", default=0, help="Distance from TSS to exclude (0 = no TSS exclusion)")
    parser.add_argument("-d", "--debug", type=str2bool, nargs="?", default=False, help="Output debugging messages")
    parser.add_argument("-v", "--verbose", type=str2bool, nargs="?", const=True, default=False, help="Print verbose messages")


    #Printing arguments to the command line
    args = parser.parse_args()

    print(f"Called with args:\n{args}\n")

    #Ensuring that argument files exist
    check_file(args.input)
    check_file(args.annot)

    return args


def main() -> None:
    """Stitch enhancer loci together

    Raises:
        ValueError: if input file is not a .bed, .gtf, .gff or gff3 file
    """

    #Parse arguments from the command line
    args = parseArgs()

    #Initialising variables
    path = get_path()
    inputGFFFile = check_path(Path(path, args.output, "gff", Path(args.input).stem)) + ".gff3"
    suffix = "_TSS_distal" if bool(int(args.tss)) else ""

    #Copying/creating the input .gff3 file
    if args.verbose:
        print(f"Converting input {Path(args.input).suffix} file to .gff3 format")
    match Path(args.input).suffix:
        case ".bed":
            bed_to_gff3(args.input, inputGFFFile)

        case ".gff":
            gff_to_gff3(args.input, inputGFFFile)

        case ".gtf":
            gtf_to_gff3(args.input, inputGFFFile, full=False)

        case ".gff3":
            check_gff(args.input, inputGFFFile)
        
        case _:
            raise ValueError("Input file must be a .bed, .gtf, .gff or .gff3 file")
    
    #Using the .gff3 file to define enhancers
    if args.verbose:
        print(f"Using {inputGFFFile} as the input .gff3 file\n")
    inputName = str(Path(inputGFFFile).stem)

    #Making the start dict
    startDict = makeStartDict(args.annot)

    #Loading enhancers as loci collection object
    referenceCollection = gffToLocusCollection(inputGFFFile)

    #Stitching regions together
    if args.verbose:
        print("Stitching regions together")
    stitchedCollection, debugOutput = regionStitching(
        referenceCollection,
        int(args.stitch),
        int(args.tss),
        startDict,
        bool(int(args.tss))
    )

    #Create a gff3 dataframe from the stitched enhancers loci collection
    if args.verbose:
        print("Making GFF from stitched collection")
    stitchedGFF = locusCollectionToGFF(stitchedCollection)
    
    #Defining output file names
    stitchedGFFFile = check_path(
        Path(path, args.output, "gff", f"{inputName}_{int(args.stitch)/1000}kb_stitched{suffix}.gff3")
    )
    debugOutFile = check_path(
        Path(path, args.output, "gff", f"{inputName}_{int(args.stitch)/1000}kb_stitched{suffix}.debug")
    )

    #Outputting the gff3 dataframe
    with open(stitchedGFFFile, "w") as f_out:
        f_out.write("##gff-version 3\n##source-version ROSE\n")
        stitchedGFF.to_csv(f_out, sep="\t", header=False, index=False, mode="a")

    #Outputting the debugging information
    if args.debug:
        pd.DataFrame(debugOutput, columns=["Enhancer", "Region", "Reason"]).to_csv(debugOutFile, index=False, sep="\t")


if __name__ == "__main__":
    main()
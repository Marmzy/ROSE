#!/usr/bin/env python

'''
PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS,
AND RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS
APRIL 11, 2013
VERSION 0.1
CONTACT: youngcomputation@wi.mit.edu
'''

import argparse
import pandas as pd

from pathlib import Path
from src.main_functions import regionStitching
from src.utils.annotation import makeStartDict
from src.utils.conversion import bed_to_gff3, check_gff, gff_to_gff3, gtf_to_gff3
from src.utils.file_helper import get_path, check_file, check_path
from src.utils.locus import gffToLocusCollection, locusCollectionToGFF

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

    parser = argparse.ArgumentParser(description='Stitch regions together to form enhancers, map read density to stitched regions and \
                                                  rank enhancers by read desnity to discover super-enhancers')

    #Required arguments
    parser.add_argument('-g', '--genome', type=str, help='Genome build (MM10, MM9, MM8, HG18, HG19, HG38)')
    parser.add_argument('-i', '--input', type=str, help='File (.bed, .gff or .gtf) containing binding sites to make enhancers')
    parser.add_argument('-o', '--output', type=str, help='Output directory name')
    parser.add_argument('-r', '--rankby',  type=str, help='File (.bam) to rank enhancers by')

    #Optional arguments
    parser.add_argument('-b', '--bams', nargs='*', help="Comma separated list of additional files (.bam) to map to")
    parser.add_argument('-c', '--control',  type=str, nargs='?', help="File (.bam) to rank enhancer by")
    parser.add_argument('-s', '--stitch', type=int, nargs='?', default=12500, help="Max linking distance for stitching")
    parser.add_argument('-t', '--tss', type=int, nargs='?', default=0, help="Distance from TSS to exclude (0 = no TSS exclusion)")
    parser.add_argument('-v', '--verbose', type=str2bool, nargs='?', const=True, default=False, help='Print verbose messages')


    #Printing arguments to the command line
    args = parser.parse_args()

    print("Called with args:")
    print(f"{args}\n")

    #Ensuring that argument files exist
    check_file(args.input)
    check_file(args.rankby)
    check_file(args.control)

    return args


def main():
    '''
    main run call
    '''

    #Parse arguments from the command line
    args = parseArgs()
    path = get_path()

    #Initialising variables
    # debug = False
    genomeDict = {
        "HG18": Path(path, "data", "annotation", "hg18_refseq.ucsc"),
        "HG19": Path(path, "data", "annotation", "hg19_refseq.ucsc"),
	    "HG38": Path(path, "data", "annotation", "hg38_refseq.ucsc"),
        "MM8":  Path(path, "data", "annotation", "mm8_refseq.ucsc"),
        "MM9":  Path(path, "data", "annotation", "mm9_refseq.ucsc"),
        "MM10": Path(path, "data", "annotation", "mm10_refseq.ucsc"),
        }
    stitchWindow = int(args.stitch)

    if args.control:        
        bamFileList = [args.rankby, args.control]
    else:
        bamFileList = [args.rankby]

    if bool(int(args.tss)):
        suffix = "_TSS_distal"
    else:
        suffix = ""

    # if options.bams:
    #     bamFileList += options.bams.split(',')
    #     bamFileLIst = ROSE_utils.uniquify(bamFileList)

    #Ensuring necessary output directories exist
    output = check_path(Path(path, args.output))
    gffFolder = check_path(Path(path, args.output, "gff"))
    mappedFolder = check_path(Path(path, args.output, "mappedGFF"))

    #Copying/creating the input .gff3 file
    if Path(args.input).suffix == ".bed":
        if args.verbose:
            print("Converting input .bed file to .gff3 format")
        inputGFFFile = str(Path(path, "output", "gff", Path(args.input).stem)) + ".gff3"
        # bed_to_gff3(args.input, inputGFFFile)

    elif Path(args.input).suffix == ".gff":
        if args.verbose:
            print("Converting input .gtf file to .gff3 format")
        inputGFFFile = str(Path(path, "output", "gff", Path(args.input).stem)) + ".gff3"
        gff_to_gff3(args.input, inputGFFFile)

    elif Path(args.input).suffix == ".gtf":
        if args.verbose:
            print("Converting input .gtf file to .gff3 format")
        inputGFFFile = str(Path(path, "output", "gff", Path(args.input).stem)) + ".gff3"
        # gtf_to_gff3(args.input, inputGFFFile, full=False)

    elif Path(args.input).suffix == ".gff3":
        if args.verbose:
            print("Checking input .gff3 file")
        inputGFFFile = str(Path(path, "output", "gff", Path(args.input).stem)) + ".gff3"
        # check_gff(args.input, inputGFFFile)
        
    else:
        raise ValueError("Input file must be a .bed, .gtf, .gff or gff3 file")

    #Using the .gff3 file to define enhancers
    if args.verbose:
        print(f"Using {inputGFFFile} as the input .gff file\n")
    inputName = str(Path(inputGFFFile).stem)

    #Making the start dict
    startDict = makeStartDict(check_file(str(genomeDict[args.genome.upper()])))

    #Loading enhancers as loci collection object
    referenceCollection = gffToLocusCollection(check_file(inputGFFFile))

    #Stitching regions together
    if args.verbose:
        print("Stitching regions together")
    stitchedCollection, debugOutput = regionStitching(referenceCollection, int(args.stitch), int(args.tss), startDict, bool(int(args.tss)))

    #Create a gff3 dataframe from the stitched enhancers loci collection
    if args.verbose:
        print("Making GFF from stitched collection")
    stitchedGFF = locusCollectionToGFF(stitchedCollection)
    
    #Defining output file names
    stitchedGFFName = f"{inputName}_{stitchWindow/1000}kb_stitched{suffix}"
    stitchedGFFFile = Path(Path(gffFolder), f"{stitchedGFFName}.gff3")
    debugOutFile = Path(Path(gffFolder), f"{inputName}_{stitchWindow/1000}kb_stitched{suffix}.debug")
    outputFile1 = Path(output, f"{stitchedGFFName}_enhancer_region_map.txt")

    #Outputting the gff3 dataframe
    with open(stitchedGFFFile, "w") as f_out:
        f_out.write("##gff-version 3\n")
        f_out.write("##source-version ROSE\n")
        stitchedGFF.to_csv(f_out, sep="\t", header=False, index=False, mode="a")

    #Outputting the debug data
    # if debug:
    #     print(f"Writing debug output to {debugOutFile}")
    #     pd.DataFrame(debugOutput).to_csv(debugOutFile, sep="\t", header=False, index=False)

    # # bin for bam mapping
    # nBin = 1

    # mappedOut1 = Path(Path(mappedFolder), f"{stitchedGFFName}_{str(Path(bamFileList[0]).name)}_mapped.gff3")
    # mappedOut2 = Path(Path(mappedFolder), f"{inputGFFFile}_{str(Path(bamFileList[1]).name)}_mapped.gff3")

    # for bamFile in bamFileList:

    #     print(bamFile)
    #     bamFileName = str(Path(bamFile).name)
    #     print(bamFileName)

    #     #MAPPING TO THE STITCHED GFF
    #     mappedOut1 = Path(Path(mappedFolder), f"{stitchedGFFName}_{bamFileName}_mapped.gff3")

        #WILL TRY TO RUN AS A BACKGROUND PROCESS. BATCH SUBMIT THIS LINE TO IMPROVE SPEED
        # cmd1 = "python ROSE_bamToGFF.py -f 1 -e 200 -r -m %s -b %s -i %s -o %s &" % (nBin, bamFile, stitchedGFFFile, mappedOut1)
        # main_bam(1, )

    #     #MAPPING TO THE ORIGINAL GFF
    #     mappedOut2 ='%s%s_%s_MAPPED.gff' % (mappedFolder,inputName,bamFileName)
    #     #WILL TRY TO RUN AS A BACKGROUND PROCESS. BATCH SUBMIT THIS LINE TO IMPROVE SPEED
    #     cmd2 = "python ROSE_bamToGFF.py -f 1 -e 200 -r -m %s -b %s -i %s -o %s &" % (nBin,bamFile,inputGFFFile,mappedOut2)
    #     print(cmd2)
    #     os.system(cmd2)
        

    
    # print('PAUSING TO MAP')
    # time.sleep(10)

    # #CHECK FOR MAPPING OUTPUT
    # outputDone = False
    # ticker = 0
    # print('WAITING FOR MAPPING TO COMPLETE. ELAPSED TIME (MIN):')
    # while not outputDone:

    #     '''
    #     check every 5 minutes for completed output
    #     '''
    #     outputDone = True
    #     if ticker%6 == 0:
    #         print(ticker*5)
    #     ticker +=1
    #     #CHANGE THIS PARAMETER TO ALLOW MORE TIME TO MAP
    #     if ticker == 144:
    #         print('ERROR: OPERATION TIME OUT. MAPPING OUTPUT NOT DETECTED')
    #         exit()
    #         break
    #     for bamFile in bamFileList:
            
    #         #GET THE MAPPED OUTPUT NAMES HERE FROM MAPPING OF EACH BAMFILE
    #         bamFileName = bamFile.split('/')[-1]
    #         mappedOut1 ='%s%s_%s_MAPPED.gff' % (mappedFolder,stitchedGFFName,bamFileName)

    #         try:
    #              mapFile = open(mappedOut1,'r')
    #              mapFile.close()
    #         except IOError:
    #             outputDone = False

    #         mappedOut2 ='%s%s_%s_MAPPED.gff' % (mappedFolder,inputName,bamFileName)
            
    #         try:
    #             mapFile = open(mappedOut2,'r')
    #             mapFile.close()
    #         except IOError:
    #             outputDone = False
    #     if outputDone == True:
    #         break
    #     time.sleep(300)
    # print('MAPPING TOOK %s MINUTES' % (ticker*5))

    # print('BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS')
    # #CALCULATE DENSITY BY REGION
    # mapCollection(stitchedCollection,referenceCollection,bamFileList,mappedFolder,outputFile1,refName = stitchedGFFName)


    # time.sleep(10)

    # print('CALLING AND PLOTTING SUPER-ENHANCERS')


    # if options.control:

    #     rankbyName = options.rankby.split('/')[-1]
    #     controlName = options.control.split('/')[-1]
    #     cmd = 'R --no-save %s %s %s %s < ROSE_callSuper.R' % (outFolder,outputFile1,inputName,controlName)

    # else:
    #     rankbyName = options.rankby.split('/')[-1]
    #     controlName = 'NONE'
    #     cmd = 'R --no-save %s %s %s %s < ROSE_callSuper.R' % (outFolder,outputFile1,inputName,controlName)
    # print(cmd)
    # os.system(cmd)



if __name__ == "__main__":
    main()
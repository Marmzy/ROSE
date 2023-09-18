#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd

from pathlib import Path
from scipy.stats import rankdata
from superenhancers.output import convert_stitched_to_bed, hockey_stick_plot, write_enhancer_table
from superenhancers.super_enhancer import calculate_cutoff
from utils.file_helper import check_file, check_path


def parseArgs() -> argparse.Namespace:
    """Parse arguments from CLI

    Returns:
        argparse.Namespace: Argparse space containing parsed arguments
    """

    parser = argparse.ArgumentParser(description="")

    #Required arguments
    parser.add_argument("-o", "--output", type=str, help="Output directory name")
    parser.add_argument("-d", "--density", type=str, help="Stitched enhancer loci signal density file")
    parser.add_argument("-g", "--gff", type=str, help="File (.bed, .gff or .gtf) containing binding sites to make enhancers")
    
    #Optional arguments
    parser.add_argument("-c", "--control",  type=str, nargs="?", help="Control (.bam) file")


    #Printing arguments to the command line
    args = parser.parse_args()

    print(f"Called with args:\n{args}\n")

    return args


def main() -> None:
    """Find superenhancers among stitched enhancer loci
    """

    #Parse arguments from the command line
    args = parseArgs()

    #Read stitched enhancer loci density signal file as dataframe
    stitched_regions = pd.read_csv(check_file(args.density), sep="\t")
    rank_cols = list(stitched_regions.columns[:6])

    #Subtract control signal if control is available
    if args.control:
        rankBy = stitched_regions.iloc[:, 6:-1].sub(stitched_regions.iloc[:, -1], axis=0)
    else:
        rankBy = stitched_regions.iloc[:, 6:]
    
    #Loop over each bam file density
    for bam in rankBy.columns:

        #Dropping other bam file densities from the dataframe
        cols = rank_cols + [bam] + [stitched_regions.columns[-1]]
        stitched_regions_temp = stitched_regions.filter(cols)

        #Setting negative values to 0
        rankBy_vector = rankBy[bam]
        rankBy_vector[rankBy_vector < 0] = 0

        #Calculate the superenhancer density signal cut-off value
        y_cutoff = calculate_cutoff(np.asarray(rankBy_vector.copy()))

        #Get superenhancers based on their density signal
        superEnhancerRows = np.where(rankBy_vector > y_cutoff)[0]

        #Create output file header
        enhancerDescription = f"{Path(args.gff).name} Enhancers\nCreated from {Path(args.density).name}"
        enhancerDescription += f"\nRanked by {bam}\nUsing cutoff of {y_cutoff} for Super-Enhancers"

        #Creating hockey stick plot
        hockey_stick_plot(np.asarray(rankBy_vector).copy(), y_cutoff, superEnhancerRows, args.output, args.gff, args.control, bam)

        #Rank stitched enhancer loci by control corrected density signal and output to .bed file
        bedName = str(Path(args.output, f"{Path(args.gff).stem}_{Path(bam).stem}_enhancers_withSuper.bed"))
        convert_stitched_to_bed(stitched_regions_temp, enhancerDescription, bedName, np.asarray(rankBy_vector), superEnhancerRows, bam)

        #Calculate stitched enhancer loci rankings and super status
        enhancer_rank = len(stitched_regions)-rankdata(rankBy_vector, method="ordinal")+1
        super_status = [1 if sr in superEnhancerRows else 0 for sr in range(0, len(stitched_regions))]
        additional_data = pd.DataFrame({"enhancerRank": enhancer_rank, "isSuper": super_status})

        #Output rankings and status dataframe
        enhancer_file = check_path(Path(args.output, f"{Path(args.gff).stem}_{Path(bam).stem}_AllEnhancers.table.txt"))
        write_enhancer_table(stitched_regions_temp, enhancerDescription, enhancer_file, additional_data)

        super_file = check_path(Path(args.output, f"{Path(args.gff).stem}_{Path(bam).stem}_SuperEnhancers.table.txt"))
        write_enhancer_table(stitched_regions_temp.iloc[superEnhancerRows, :], enhancerDescription, super_file, additional_data.iloc[superEnhancerRows, :])



if __name__ == "__main__":
    main()
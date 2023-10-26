#!/usr/bin/env python

import numpy as np
import pandas as pd

from pathlib import Path
from scipy.stats import rankdata
from src.superenhancers.output import (
    convert_stitched_to_bed,
    hockey_stick_plot,
    write_enhancer_table
)
from src.superenhancers.super_enhancer import calculate_cutoff
from src.utils.file_helper import check_file, check_path
from typing import Union


def get_super(
    output: str,
    density: str,
    gff: str,
    control: Union[str, None],
) -> None:
    """Find superenhancers among stitched enhancer loci
    """

    # Read stitched enhancer loci density signal file as dataframe
    stitched_regions = pd.read_csv(check_file(density), sep="\t")
    rank_cols = list(stitched_regions.columns[:6])


    # Subtract control signal if control is available
    if control:
        rankBy = stitched_regions.iloc[:, 6:-1].sub(
            stitched_regions.iloc[:, -1],
            axis=0
        )
    else:
        rankBy = stitched_regions.iloc[:, 6:]


    # Loop over each bam file density
    for bam in rankBy.columns:

        # Dropping other bam file densities from the dataframe
        cols = rank_cols + [bam] + [stitched_regions.columns[-1]]
        stitched_regions_temp = stitched_regions.filter(cols)

        # Setting negative values to 0
        rankBy_vector = rankBy[bam].copy()
        rankBy_vector[rankBy_vector < 0] = 0

        # Calculate the superenhancer density signal cut-off value
        y_cutoff = calculate_cutoff(np.asarray(rankBy_vector.copy()))

        # Get superenhancers based on their density signal
        superEnhancerRows = np.where(rankBy_vector > y_cutoff)[0]

        # Creating hockey stick plot
        hockey_stick_plot(
            np.asarray(rankBy_vector).copy(),
            y_cutoff,
            superEnhancerRows,
            output,
            gff,
            control,
            bam
        )

        # Create output file header
        enhancerDescription = f"{Path(gff).name} Enhancers\nCreated from {Path(density).name}"
        enhancerDescription += f"\nRanked by {bam}\nUsing cutoff of {y_cutoff} for Super-Enhancers"

        # Rank stitched enhancer loci by control corrected density signal and
        # output to .bed file
        bedName = str(Path(output, f"{Path(gff).stem}_{Path(bam).stem}_enhancers_withSuper.bed"))
        convert_stitched_to_bed(
            stitched_regions_temp,
            enhancerDescription,
            bedName,
            np.asarray(rankBy_vector),
            superEnhancerRows,
            bam
        )

        # Calculate stitched enhancer loci rankings and super status
        enhancer_rank = len(stitched_regions)-rankdata(rankBy_vector, method="ordinal")+1
        super_status = [
            1 if sr in superEnhancerRows else 0
            for sr in range(0, len(stitched_regions))
        ]
        additional_data = pd.DataFrame(
            {"enhancerRank": enhancer_rank, "isSuper": super_status}
        )

        # Output rankings and status dataframe
        enhancer_file = check_path(
            Path(output, f"{Path(gff).stem}_{Path(bam).stem}_AllEnhancers.table.txt")
        )
        write_enhancer_table(
            stitched_regions_temp,
            enhancerDescription,
            enhancer_file,
            additional_data
        )

        super_file = check_path(
            Path(output, f"{Path(gff).stem}_{Path(bam).stem}_SuperEnhancers.table.txt")
        )
        write_enhancer_table(
            stitched_regions_temp.iloc[superEnhancerRows, :],
            enhancerDescription,
            super_file,
            additional_data.iloc[superEnhancerRows, :]
        )

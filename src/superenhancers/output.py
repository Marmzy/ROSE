#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from datetime import datetime
from pathlib import Path
from scipy.stats import rankdata
from typing import Union
from utils.file_helper import check_path


def convert_stitched_to_bed(
    stitchedRegions: pd.core.frame.DataFrame,
    trackDescription: str,
    output: str,
    densitySignal: np.ndarray,
    superRows: np.ndarray,
    bam: str
) -> None:
    """Create a .bed file ranking stitched enhancer loci by their (corrected) density signal

    Args:
        stitchedRegions (pd.core.frame.DataFrame): Stitched enhancer loci signal density dataframe
        trackDescription (str): Analysis description
        output (str): Output file name
        densitySignal (np.ndarray): (Control corrected) density signal values
        superRows (np.ndarray): Vector of superenhancer indices
        bam (str): .bam file name
    """

    #Create output dataframe
    df = stitchedRegions.loc[:, ["CHROM", "START", "STOP", "REGION_ID"]]

    #Rank density signal-corrected stitched enhancer loci from smallest to largest
    if len(densitySignal) == len(stitchedRegions):
        densitySignal = rankdata(densitySignal, method="ordinal")
        densitySignal = len(densitySignal) - densitySignal +1
        df["RANK"] = densitySignal

    #Create a description track and output the dataframe to a .bed file
    trackDescription = f"{trackDescription}\nCreated on {datetime.now().strftime('%Y-%m-%d')}".replace("\n", "\t")
    with open(output, "w") as f_out:
        f_out.write(f"track name='{bam}_Enhancers' description='{trackDescription}' itemRGB=On color=0,0,0\n")
        df.iloc[superRows, :].to_csv(f_out, sep="\t", header=False, index=False, mode="a")


def hockey_stick_plot(
    vector: np.ndarray,
    y_cutoff: float,
    super_enhancer_rows: np.ndarray,
    out: str,
    gff: str,
    control: Union[str, None],
    bam: str
) -> None:
    """Visualise stitched enhancer categorisation with a hockey-stick plot

    Args:
        vector (np.ndarray): Array of density signal values
        y_cutoff (np.float): Density signal cut-off value
        super_enhancer_rows (np.ndarray): Array of super-enhancer row indices
        out (str): Output directory name
        gff (str): .gff file name
        control (Union[str, None]): Control file name
        bam (str): .bam file name
    """

    #Ordering the stitched enhancer loci signal density vector from low to high
    vector.sort()

    #Creating the hockey stick plot
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(list(range(len(vector))), vector, "o-", color="red")
    ax.axhline(y_cutoff, color="grey", linestyle="dashed")
    ax.axvline(len(vector)-len(super_enhancer_rows), color="grey", linestyle="dashed")
    ax.set_xlabel(f"{bam} Stitched Enhancers")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(1,4))
    if control:
        ax.set_ylabel(f"{bam} - {Path(control).name} Signal")
    else:
        ax.set_ylabel(f"{bam} Signal")
    ax.set_title(f"Cut-off used: {y_cutoff}\nSuper-Enhancers identified: {len(super_enhancer_rows)}")
    fig.tight_layout()
    fig.savefig(check_path(Path(out, f"{Path(gff).stem}_{Path(bam).stem}_plot_points.png")))
    

def write_enhancer_table(
    enhancer_df: pd.core.frame.DataFrame,
    description: str,
    output: str,
    additional_data: pd.core.frame.DataFrame
) -> None:
    """Output sorted stitched enhancer loci regions data

    Args:
        enhancer_df (pd.core.frame.DataFrame): Stitched enhancer loci dataframe
        description (str): Data description
        output (str): Output file name
        additional_data (pd.core.frame.DataFrame): Stitched enhancer loci rankings and super status dataframe

    Raises:
        ValueError: If the lengths of the two dataframes don't match
    """

    #Creating a descriptive header
    description = f"#{description}\nCreated on {datetime.now().strftime('%Y-%m-%d')}".replace("\n", "\n#")

    #Merge stitched enhancer loci regions data with superenhancer ranking data
    if len(additional_data) == len(enhancer_df):
        out_df = pd.concat([enhancer_df, additional_data], axis=1).sort_values("enhancerRank")
        with open(output, "w") as f_out:
            f_out.write(f"{description}\n")
            out_df.to_csv(f_out, sep="\t", index=False, mode="a")
    else:
        raise ValueError("Stitched enhancer loci dataframe and additional data have differing lengths")
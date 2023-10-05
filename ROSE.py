#!/usr/bin/env python3

import argparse

from src.ROSE_bamToGFF import calc_read_density
from src.ROSE_main import stitch_loci
from src.utils.file_helper import check_file, read_yaml


def parseArgs() -> argparse.Namespace:
    """Parse arguments from CLI

    Returns:
        argparse.Namespace: Argparse space containing parsed arguments
    """

    parser = argparse.ArgumentParser(description="Run ROSE from start to finish")

    #Required arguments
    parser.add_argument("-c", "--config", type=str, help="Configuration .yaml file")

    #Printing arguments to the command line
    args = parser.parse_args()
    check_file(args.config)

    return args


def main():

    #Parse arguments from the command line
    args = parseArgs()

    #Read the input configuration file
    conf = read_yaml(args.config)

    #Stitch enhancer loci
    original, stitched = stitch_loci(
        input = check_file(conf["data"]["input"]),
        output = conf["data"]["output"],
        annot = check_file(conf["data"]["annotation"]),
        **conf["stitching"],
        verbose = conf["verbose"],
    )

    #Map reads and calculate read density per stitched
    #enhancer locus bin
    calc_read_density(
        conf["data"]["rankby"],
        stitched,
        original,
        check_file(conf["data"]["control"]),
        **conf["mapping"],
        verbose = conf["verbose"],
    )


if __name__ == "__main__":
    main()
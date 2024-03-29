#!/usr/bin/env python3

import argparse

from pathlib import Path
from src.ROSE_bamToGFF import calc_read_density
from src.ROSE_callSuper import get_super
from src.ROSE_main import stitch_loci
from src.ROSE_mapCollection import map_collection
from src.utils.file_helper import check_file, get_path, read_yaml


def parseArgs() -> argparse.Namespace:
    """Parse arguments from CLI

    Returns:
        argparse.Namespace: Argparse space containing parsed arguments
    """

    parser = argparse.ArgumentParser(description="Run ROSE from start to finish")

    # Required arguments
    parser.add_argument("-c", "--config", type=str, help="Configuration .yaml file")

    # Printing arguments to the command line
    args = parser.parse_args()
    check_file(args.config)

    return args


def main():

    # Parse arguments from the command line
    args = parseArgs()

    # Read the input configuration file
    path = get_path()
    conf = read_yaml(args.config)

    # Check existence of input files
    check_file(conf["data"]["input"])
    check_file(conf["data"]["annotation"])
    if conf["data"]["control"]:
        check_file(conf["data"]["control"])

    # Stitch enhancer loci
    original, stitched = stitch_loci(
        input=conf["data"]["input"],
        output=conf["data"]["output"],
        annot=conf["data"]["annotation"],
        **conf["stitching"],
        verbose=conf["verbose"],
    )

    # Map reads to each stitched enhancer locus bin
    calc_read_density(
        conf["data"]["rankby"],
        stitched,
        original,
        conf["data"]["control"],
        **conf["mapping"],
        verbose=conf["verbose"],
    )

    # Calculate read density signal for each stitched enhancer locus
    density = map_collection(
        stitched,
        original,
        conf["data"]["rankby"],
        conf["data"]["control"],
        Path(path, conf["data"]["output"], "mappedGFF")
    )

    # Identifying and visualising superenhancers
    get_super(
        Path(path, conf["data"]["output"]),
        density,
        conf["data"]["input"],
        conf["data"]["control"]
    )


if __name__ == "__main__":
    main()

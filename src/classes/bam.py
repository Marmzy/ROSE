#!/usr/bin/env python

import re
import subprocess

from src.classes.locus import Locus
from typing import List


class Bam:

    def __init__(
        self,
        file: str
    ) -> None:
        """Create Bam class object

        Args:
            file (str): .bam file
        """

        self._bam = file

    def checkChrStatus(
        self
    ) -> None:
        """Check the chromosome naming in the .bam file
        """

        # Run samtools view and capture output
        cmd = f"samtools view {self._bam} | head -n 1"
        stdout = subprocess.run(cmd, capture_output=True, shell=True, check=True, text=True).stdout

        # Capture chromosome naming from output
        self._chr = 1 if "chr" in stdout.split()[2] else 0

    def convertBitwiseFlag(
        self,
        flag: str
    ) -> str:
        """Get read strand information from FLAG field

        Args:
            flag (str): .bam file FLAG field

        Returns:
            str: Read strand
        """

        return "-" if flag == "16" else "+"

    def getTotalReads(
        self
    ) -> int:
        """Get total mapped alignments from the .bam file

        Returns:
            int: Total mapped alignments
        """

        # Run samtools flagstat and capture output
        cmd = f"samtools flagstat {self._bam}"
        stdout = subprocess.run(
            cmd, capture_output=True, shell=True, check=True, text=True
        ).stdout.splitlines()

        # Return total mapped alignments from output
        for line in stdout:
            if "mapped (" in line:
                return int(line.split(" ")[0])

    def getRawReads(
        self,
        locus: Locus
    ) -> List[List[str]]:
        """Get all reads in the extended sitched enhancer locus region

        Args:
            locus (Locus): Extended stiched enhancer locus

        Returns:
            List[List[str]]: List of all reads within the region
        """

        # Run samtools view and capture output
        locusLine = f"{locus._chr}:{str(locus._start)}-{str(locus._end)}"
        cmd = f"samtools view {self._bam} {locusLine}"
        stdout = subprocess.run(
            cmd, capture_output=True, shell=True, check=True, text=True
        ).stdout.splitlines()

        # Get all reads in the stitched enhancer locus region
        return [line.split("\t") for line in stdout]

    def readsToLoci(
        self,
        reads: List[List[str]]
    ) -> List[Locus]:
        """Convert reads to Locus objects

        Args:
            reads (List[List[str]]): List of reads within the extended
                                     stitched enhancer region

        Returns:
            List[Locus]: List of Locus objects
        """

        # Initialising variables
        loci = []

        # Looping over the reads and converting them to Locus objects
        for read in reads:
            strand = self.convertBitwiseFlag(read[1])

            # Remove skipped nucleotides in read
            if "N" in read[5]:
                length, total = 0, 0
                for c in re.findall(r"\d+\w", read[5]):
                    if c[-1] != "N":
                        length += int(c[:-1])
                    else:
                        loci.append(
                            Locus(read[2], int(read[3])+total,
                                  int(read[3])+total+length, strand, "")
                        )
                        total = total + length + int(c[:-1])
                        length = 0
                if c[-1] != "N":
                    loci.append(
                        Locus(read[2], int(read[3])+total, int(read[3])+total+length, strand, "")
                    )
            else:
                loci.append(Locus(read[2], int(read[3]), int(read[3])+len(read[9]), strand, ""))

        return loci

    def getReadsLocus(
        self,
        locus: Locus
    ) -> List[Locus]:
        """Get read loci within the extended stitched enhancer locus region

        Args:
            locus (Locus): Extended stitched enhancer locus

        Returns:
            List[Locus]: List of read loci
        """

        # Get reads within the extended stitched enhancer locus region
        reads = self.getRawReads(locus)

        # Convert reads to loci
        loci = self.readsToLoci(reads)

        return loci

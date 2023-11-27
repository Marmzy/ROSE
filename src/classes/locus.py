#!/usr/bin/env python
from __future__ import annotations

import pandas as pd

from typing import List, Tuple, TypedDict, Union


class Locus:
    def __init__(
        self,
        chr: str,
        start: int,
        end: int,
        sense: str,
        ID: str = ""
    ) -> None:
        """Create locus region object

        Args:
            chr (str): Chromosome name
            start (int): Locus start coordinate
            end (int): Locus end coordinate
            sense (str): Locus' strand
            ID (str, optional): Unique locus ID. Defaults to "".
        """

        self._chr = chr
        self._sense = sense
        self._start = start
        self._end = end
        self._ID = ID

    def contains(
        self,
        otherLocus: Locus
    ) -> bool:
        """Check if a locus lies completely within an other locus

        Args:
            otherLocus (Locus): Locus to compare with

        Returns:
            bool: Bool delineating if one locus lies entirely within another
        """

        if self._chr != otherLocus._chr:
            return False
        elif self._sense != "." and otherLocus._sense != "." and self._sense != otherLocus._sense:
            return False
        elif self._start > otherLocus._start or otherLocus._end > self._end:
            return False
        else:
            return True

    def getAntisenseLocus(
        self
    ) -> Locus:
        """Create an antisense Locus for the given Locus

        Returns:
            Locus: Locus object with an antisense strand
        """

        if self._sense == ".":
            return self
        switch = {"+": "-", "-": "+"}
        return Locus(self._chr, self._start, self._end, switch[self._sense], self._ID)

    def overlaps(
        self,
        otherLocus: Locus
    ) -> bool:
        """Check if two loci overlapping coordinates

        Args:
            otherLocus (Locus): Locus to compare coordinates with

        Returns:
            bool: Bool delineating if two loci overlap
        """

        if self._chr != otherLocus._chr:
            return False
        elif not self._sense == "." or otherLocus._sense == "." or self._sense == otherLocus._sense:
            return False
        elif self._start > otherLocus._end or otherLocus._start > self._end:
            return False
        else:
            return True

    def __eq__(
        self,
        other: Locus
    ) -> bool:
        """Check if two loci have the same chromosome, strand, start
           and end positions

        Args:
            other (Locus): Locus to compare with

        Returns:
            bool: Bool whether two loci are "equal"
        """

        if self.__class__ != other.__class__:
            return False
        if self._chr != other._chr:
            return False
        if self._start != other._start:
            return False
        if self._end != other._end:
            return False
        if self._sense != other._sense:
            return False
        return True

    def __hash__(
        self
    ) -> int:
        """Create hash for tracking entries

        Returns:
            int: Hash ID
        """

        return self._start + self._end

    def __len__(
        self
    ) -> int:
        """Calculate locus length

        Returns:
            int: Length of the locus
        """

        return self._end - self._start + 1

    def __ne__(
        self,
        other: Locus
    ) -> bool:
        """Check if two loci do not have the same chromosome, strand, start
           and end positions

        Args:
            other (Locus): Locus to compare with

        Returns:
            bool: Bool whether two loci are not "equal"
        """

        return not self.__eq__(other)

    def __str__(
        self
    ) -> str:
        """Return locus name

        Returns:
            str: Locus name
        """

        return f"{self._chr}({self._sense}):{self._start}-{self._end}"


class LocusCollection:
    def __init__(
        self,
        loci: List[Locus],
        windowSize: int
    ) -> None:
        """init

        Args:
            loci (List[Locus]): List of locus regions
            windowSize (int): Size that regions will be broken down into
        """

        self._chrToCoordToLoci = {}
        self._loci = {}
        self._winSize = windowSize
        for locus in loci:
            self.__addLocus(locus)

    def append(
        self,
        new: Locus
    ) -> None:
        """Add new locus to collection

        Args:
            new (Locus): Locus to add
        """

        self.__addLocus(new)

    def remove(
        self,
        old: Locus
    ) -> None:
        """Remove locus from loci the LocusCollection

        Args:
            old (Locus): Locus to remove

        Raises:
            ValueError: If locus is not in the loci collection
        """

        # Check that locus is actually present
        if old not in self._loci.keys():
            raise ValueError("Locus to remove is not in the collection")

        # Removing the locus from the loci dictionary
        del self._loci[old]

        # Removing the locus from the KeyRange dictionary
        senseList = ["+", "-"] if old._sense == "." else [old._sense]

        for k in self.__getKeyRange(old):
            for sense in senseList:
                self._chrToCoordToLoci[old._chr+sense][k].remove(old)

    def getContainers(
        self,
        locus: Locus,
        sense: str = "sense"
    ) -> TypedDict:
        """Get TSS loci that envelop the enhancer locus

        Args:
            locus (Locus): Enhancer locus
            sense (str, optional): TSS loci strands. Defaults to "sense".

        Returns:
            TypedDict: Dictionary keys, with TSS loci that envelop the enhancer locus as keys
        """

        # Get all loci from the TSS collection that (partly) overlap
        # with the enhancer locus
        matches = self.__subsetHelper(locus, sense)

        # Get all TSS loci that envelop the enhancer locus
        if sense in ["sense", "both"]:
            realMatches = {match: None for match in matches if match.contains(locus)}
        if sense in ["antisense", "both"]:
            realMatches = {
                match: None for match in matches if match.getAntisenseLocus().contains(locus)
            }

        return realMatches.keys()

    def getLoci(
        self
    ) -> TypedDict:
        """Get all loci in the collection

        Returns:
            TypedDict: Dict keys, with loci as keys
        """

        return self._loci.keys()

    def getOverlap(
        self,
        locus: Locus,
        sense: str = "sense"
    ) -> set:
        """Get enhancer loci that overlap with a given enhancer locus

        Args:
            locus (Locus): Locus against which overlapping loci are to be searched for
            sense (str, optional): Locus strand. Defaults to "sense".

        Returns:
            set: Set of overlapping loci
        """

        # Get all loci from the enhancer collection that overlap with
        # the enhancer locus
        matches = self.__subsetHelper(locus, sense)

        # Remove loci that don't really overlap
        if sense in ["sense", "both"]:
            realMatches = {lcs for lcs in matches if lcs.overlaps(locus)}
        if sense in ["antisense", "both"]:
            realMatches = {lcs for lcs in matches if lcs.getAntisenseLocus().overlaps(locus)}

        return realMatches

    def stitchCollection(
        self,
        stitchWindow: int = 1,
        sense: str = "both"
    ) -> LocusCollection:
        """Stitch enhancer loci together into a LocusCollection

        Args:
            stitchWindow (int, optional): Window size to stitch enhancers
                                          together. Defaults to 1.
            sense (str, optional): LocusCollection strand. Defaults to "both".

        Returns:
            LocusCollection: Collection of stichted enhancer loci
        """

        # Initialising variables
        locusList = self.getLoci()
        oldCollection = LocusCollection(locusList, 500)
        stitchedCollection = LocusCollection([], 500)

        # Loop over all enhancers
        for locus in locusList:
            if locus in oldCollection._loci:
                oldCollection.remove(locus)
                overlappingLoci = oldCollection.getOverlap(
                    Locus(locus._chr, locus._start-stitchWindow, locus._end+stitchWindow,
                          locus._sense, locus._ID),
                    sense
                )

                # Loop over enhancers that overlap with the target enhancer and
                # create stitched enhancer locus
                stitchTicker = 1
                while len(overlappingLoci) > 0:
                    stitchTicker += len(overlappingLoci)
                    overlapCoords = [locus._start, locus._end]

                    for overlappingLocus in overlappingLoci:
                        overlapCoords += [overlappingLocus._start, overlappingLocus._end]
                        oldCollection.remove(overlappingLocus)
                    if sense == "both":
                        locus = Locus(locus._chr, min(overlapCoords), max(overlapCoords),
                                      ".", locus._ID)
                    else:
                        locus = Locus(locus._chr, min(overlapCoords), max(overlapCoords),
                                      locus._sense, locus._ID)
                    overlappingLoci = oldCollection.getOverlap(
                        Locus(locus._chr, locus._start-stitchWindow, locus._end+stitchWindow,
                              locus._sense),
                        sense
                    )

                # Add the stiched enhancer locus to a new LocusCollection
                locus._ID = f"{stitchTicker}_{locus._ID}_lociStitched"
                stitchedCollection.append(locus)

        return stitchedCollection

    def __addLocus(
        self,
        locus: Locus
    ) -> None:
        """Add locus data to dictionary in the following order:
           chromosome+strand, region, Locus object

        Args:
            locus (Locus): Enhancer locus region
        """

        # Add unique enhancer locus regions
        if locus not in self._loci:
            self._loci[locus] = None
            if locus._sense == ".":
                chrKeyList = [f"{locus._chr}+", f"{locus._chr}-"]
            else:
                chrKeyList = [f"{locus._chr}{locus._sense}"]
            for chrKey in chrKeyList:
                if chrKey not in self._chrToCoordToLoci:
                    self._chrToCoordToLoci[chrKey] = dict()
                for n in self.__getKeyRange(locus):
                    if n not in self._chrToCoordToLoci[chrKey]:
                        self._chrToCoordToLoci[chrKey][n] = []
                    self._chrToCoordToLoci[chrKey][n].append(locus)

    def __getKeyRange(
        self,
        locus: Locus
    ) -> range:
        """Break locus region down into a key range

        Args:
            locus (Locus): Enhancer locus region

        Returns:
            range: Locus region key range
        """

        # Create range based on the amount of window sizes that can
        # fit in the locus region
        start = locus._start // self._winSize
        end = locus._end // self._winSize + 1
        return range(start, end)

    def __len__(
        self
    ) -> int:
        """Calculate number of loci in collection

        Returns:
            int: Number of loci
        """

        return len(self._loci)

    def __subsetHelper(
        self,
        locus: Locus,
        sense: str
    ) -> TypedDict:
        """Fetch list of loci, given input

        Args:
            locus (Locus): Enhancer locus object
            sense (str): Strand on which the locus is located

        Raises:
            ValueError: If sense value is impossible

        Returns:
            TypedDict[Locus]: Subsetted list of loci
        """

        # Initialising variables
        matches = dict()
        senses = ["+", "-"]

        # Create filter to select enhancer loci
        if locus._sense == "." or sense == "both":
            lamb = lambda s: True
        elif sense == "sense":
            lamb = lambda s: s == locus._sense
        elif sense == "antisense":
            lamb = lambda s: s != locus._sense
        else:
            raise ValueError(f"Inappropriate sense value: '{sense}'")

        # Select enhancer loci based on filter
        for s in filter(lamb, senses):
            chrKey = locus._chr + s
            if chrKey in self._chrToCoordToLoci:
                for n in self.__getKeyRange(locus):
                    if n in self._chrToCoordToLoci[chrKey]:
                        for lcs in self._chrToCoordToLoci[chrKey][n]:
                            matches[lcs] = None

        return matches.keys()


def gffToLocusCollection(
    infile: str
) -> LocusCollection:
    """Create a collection of Locus objects for each entry in the .gff3 file

    Args:
        infile (str): Input .ggf3 file path

    Raises:
        ValueError: If non-unique identifiers are found among .gff3 entries

    Returns:
        LocusCollection: Collection of all .gff3 entries as loci
    """

    # Reading the input file as a dataframe
    gff_df = pd.read_csv(infile, sep="\t", header=None, comment='#')
    gff = gff_df.iloc[:, [0, 3, 4, 6, 8]].copy()

    # Converting each entry to a Locus object
    lociList = [Locus(*row) for row in zip(*gff.to_dict("list").values())]

    # Checking that entries have unique names
    if len(gff.iloc[:, 4]) != len(set(gff.iloc[:, 4])):
        raise ValueError(
            "Last column (attributes column) of the .gff3 file must contain \
                unique identifiers for each entry"
        )

    return LocusCollection(lociList, 500)


def locusCollectionToGFF(
    locusCollection: LocusCollection
) -> pd.DataFrame:
    """Create .gff3 file from stitched enhancer loci collection

    Args:
        locusCollection (LocusCollection): Stitched enhancer loci collection

    Returns:
        pd.DataFrame: .gff3 formatted dataframe
    """

    lociList = locusCollection.getLoci()
    gff_df = pd.DataFrame(
        [(locus._chr, "ROSE", "stitched_enhancer_locus", locus._start,
          locus._end, ".", locus._sense, ".", locus._ID)
         for locus in lociList
         ]
    )

    return gff_df


def makeTSSLocus(
    entry: Tuple[Union[int, str]],
    window: int
) -> Locus:
    """Create locus object for genes

    Args:
        entry (Tuple[Union[int, str]]): Gene entry
        window (int): Window to expand gene site with TSS

    Returns:
        Locus: Gene locus
    """

    # Creating the locus object
    return Locus(entry[1], entry[3]-window, entry[3]+window, entry[2], entry[0])

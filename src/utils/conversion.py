#!/usr/bin/env python

import pandas as pd
import re
import shutil


def bed_to_gff3(
    input: str,
    output: str
) -> None:
    """Convert .bed file to .gff3 file (modeled after R's rtracklayer)

    Args:
        input (str): Input .bed file path
        output (str): Output .gff3 file path
    """

    #Reading the .bed file as a dataframe
    bed_df = pd.read_csv(input, sep="\t", header=None, comment="#")

    #Converting the .bed dataframe to a .gff3 dataframe
    gff_df = pd.DataFrame({
        "seqid": bed_df.iloc[:, 0],
        "source": ["ROSE"]*len(bed_df),
        "type": ["sequence_feature"]*len(bed_df),
        "start": bed_df.iloc[:, 1]+1,
        "end": bed_df.iloc[:, 2],
        "score": bed_df.iloc[:, 4],
        "strand": bed_df.iloc[:, 5],
        "phase": ["."]*len(bed_df),
        "attributes": "name="+bed_df.iloc[:, 3]
    })

    #Outputting the gff3 dataframe
    with open(output, "w") as f_out:
        f_out.write("##gff-version 3\n##source-version ROSE\n")
        gff_df.to_csv(f_out, sep="\t", header=False, index=False, mode="a")


def check_gff(
    input: str,
    output: str
) -> None:
    """Check if .gff3 file's seqid column values start with 'chr'

    Args:
        input (str): Input .ggf3 file path
        output (str): Output .gff3 file path

    Raises:
        ValueError: ValueError if seqid column's values don't start with 'chr', nor are Human NCBI refseq accesion numbers
    """

    #Initialising variables
    RefseqAccHuman = {
        "NC_000001": "chr1",
        "NC_000002": "chr2",
        "NC_000003": "chr3",
        "NC_000004": "chr4",
        "NC_000005": "chr5",
        "NC_000006": "chr6",
        "NC_000007": "chr7",
        "NC_000008": "chr8",
        "NC_000009": "chr9",
        "NC_000010": "chr10",
        "NC_000011": "chr11",
        "NC_000012": "chr12",
        "NC_000013": "chr13",
        "NC_000014": "chr14",
        "NC_000015": "chr15",
        "NC_000016": "chr16",
        "NC_000017": "chr17",
        "NC_000018": "chr18",
        "NC_000019": "chr19",
        "NC_000020": "chr20",
        "NC_000021": "chr21",
        "NC_000022": "chr22",
        "NC_000023": "chrX",
        "NC_000024": "chrY",
    }

    #Reading the .gff3 file as a dataframe
    df = pd.read_csv(input, sep="\t", header=None, comment="#")

    #If .gff3 seqid column values do not start with "chr", see if they are Human NCBI refseq accession numbers and replace them
    if "chr" not in df.iloc[:, 0].values:
        df.iloc[:, 0] = df.iloc[:, 0].replace(RefseqAccHuman, regex=True)
        chr_df = df.loc[df.iloc[:, 0].str.startswith('chr', na=False)]

        #Outputting the gff3 dataframe
        if len(df) == len(chr_df):
            with open(output, "w") as f_out:
                f_out.write("##gff-version 3\n##source-version ROSE\n")
                df.to_csv(f_out, sep="\t", header=False, index=False, mode="a")
        else:
            raise ValueError(f"Input file {input}'s seqid column values must start with 'chr'")
    else:
        #Copy the .gff3 file if chromosome names are in correct format
        shutil.copyfile(input, output)


def gff_to_gff3(
    input: str,
    output: str
) -> None:

    #Reading the .gff file as a dataframe
    df = pd.read_csv(input, sep="\t", header=None, comment="#")

    #Filling missing values 
    df.iloc[:, 1] = ["ROSE"] * len(df)
    df.iloc[:, 2].fillna("sequence_feature", inplace=True)
    df.iloc[:, 5].fillna(".", inplace=True)
    df.iloc[:, 7].fillna(".", inplace=True)
    df.iloc[:, 8].fillna(f"{df.iloc[:, 0]}:{df.iloc[:, 6]}:{df.iloc[:, 3]}-{df.iloc[:, 4]}", inplace=True)

    #Outputting the gff3 dataframe
    with open(output, "w") as f_out:
        f_out.write("##gff-version 3\n##source-version ROSE\n")
        df.to_csv(f_out, sep="\t", header=False, index=False, mode="a")

    
def gtf_to_gff3(
    input: str,
    output: str,
    full: bool = False
) -> None:
    """Convert .gtf file to .gff3 file (modeled after R's rtracklayer)

    Args:
        input (str): Input .bed file path
        output (str): Output .gff3 file path
        full (bool, optional): Boolean to fully convert the .gtf atttributes to .gff3 format. Defaults to False.
    """

    #Initialising variables
    gtf_attributes = ["gene_id", "db_xref", "gbkey", "gene", "gene_biotype", "transcript_id",
                      "model_evidence", "product", "exon_number", "protein_id", "anticodon",
                      "inference", "note", "exception", "transl_except", "pseudo", "partial"]

    #Reading the .gtf file as a dataframe
    df = pd.read_csv(input, sep="\t", header=None, comment="#")
    df.iloc[:, 1] = ["ROSE"] * len(df)
    if full:
        df.iloc[:, 8] = [";".join([a for a in [search(attr, annot) for attr in gtf_attributes] if a is not None]) for annot in df.iloc[:, 8].values]
    else:
        df.iloc[:, 8] = [";".join([a for a in [search(attr, annot) for attr in ["gene_id", "transcript_id"]] if a is not None]) for annot in df.iloc[:, 8].values]

    #Outputting the gff3 dataframe
    with open(output, "w") as f_out:
        f_out.write("##gff-version 3\n##source-version ROSE\n")
        df.to_csv(f_out, sep="\t", header=False, index=False, mode="a")


def search(
    attr: str,
    entry: str
) -> str:
    """Checks if an attribute is present for a given .gtf entry

    Args:
        attr (str): .gtf attribute
        entry (str): .gtf attributes entry

    Returns:
        str: Attribute=value concatenation
    """

    #Check if attribute is present for given entry, and if it is, return it's value
    if re.search(f'{attr} "', entry):
        return "{}={}".format(attr, re.findall(rf"(?<={attr} \")[^\"]+", entry)[0])
## ROSE

### Description

Rank Ordering of Super-Enhancers aka ROSE is a tool for identifying super-enhancers. It does this by separating super-enhancers from typical enhancers using sequencing data (.bam) given a file of previsously identified constituent enhancers (.gff). The original ROSE tool was developed by Charles Y. Lin, David A. Orlando and Brian J. Abraham at Young Lab Whitehead Institute/MIT. This new ROSE version is an attempt to update the code from Python 2 to 3, use newer versions of tools, make the code more readable to allow for better in-depth understanding of the algorithm and to increase the computational speed.

---

### Citation

*Master Transcription Factors and Mediator Establish Super-Enhancers at Key Cell Identity Genes*
Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young. [Cell](https://www.sciencedirect.com/science/article/pii/S0092867413003929) 153, 307-319, April 11, 2013

and

*Selective Inhibition of Tumor Oncogenes by Disruption of Super-enhancers* 
Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando, Christopher R. Vakoc, James E. Bradner, Tong Ihn Lee, and Richard A. Young. [Cell](https://www.sciencedirect.com/science/article/pii/S0092867413003930) 153, 320-334, April 11, 2013

---

### Requirements

- .bam files of sequencing reads for factor of interest (reads for control is recommended, but optional).
	- .bam files must have chromosome IDs starting with "chr"
	- .bam files must be sorted and indexed using [SAMtools](http://www.htslib.org/doc/samtools.html)

- .gff file of constituent enhancers previously identified.
	- Input .gff file can be [.gff3 format](https://asia.ensembl.org/info/website/upload/gff3.html) (recommended), [.gtf format](https://asia.ensembl.org/info/website/upload/gff.html), .gff format or [.bed format](https://asia.ensembl.org/info/website/upload/bed.html)
	- If .gff format, the file must have the following columns:
		* chromosome (chr#)
		* unique ID for each constituent enhancer region
		* start of constituent
		* end of constituent
		* strand (+,-,.)
		* unique ID for each constituent enhancer region
		
- .ucsc annotation file in UCSC table track format (https://genome.ucsc.edu/cgi-bin/hgTables)

---

### Usage

```
Usage: ROSE.sh [-h help] [-g genome] [-i input] [-o output] [-r rankby] [optional flags]
 #Required arguments
 -g, --genome     Genome build (MM8, MM9, MM10, HG18, HG19, HG38)
 -i, --input      File (.bed, .gff or .gtf) containing enhancer binding sites
 -o, --output     Name of output directory where data will be stored
 -r, --rankby     .bam file to rank enhancers by

 #Additional arguments
 -v, --verbose    Print verbose messages (default=true)

 #Additional arguments for ROSE_main.py
 -c, --control    .bam file to rank enhancers by
 -s, --stitch     Max linking distance for stitching (default=12500)
 -t, --tss        Distance from TSS to exclude (0 = no TSS exclusion) (default=0)

 #Additional arguments for ROSE_bamToGFF.py
 -n, --sense      Strand to map to (default='both')
 -f, --floor      Read floor threshold necessary to count towards density (default=1)
 -x, --extension  Extends reads by n bp (default=200)
 -p, --rpm        Normalizes density to reads per million (rpm) (default=true)
 -m, --matrix     Variable bin sized matrix (default=1)

Example: ROSE.sh -g hg18 -i ./data/HG18_MM1S_MED1.gff -o output -r ./data/MM1S_MED1.hg18.bwt.sorted.bam -c ./data/MM1S_WCE.hg18.bwt.sorted.bam -s 12500 -t 2500 -n both -x 200 -p true -m 1 -v true
```

`ROSE.sh` will run ROSE from start to end. It will (amongst others) call the following scripts:
	- `ROSE_main.py`: Stitches regions together to form enhancers
	- `ROSE_bamToGFF.py`: Map .bam reads to stitched enhancers and calculate read density
	- `ROSE_mapCollection.py`: Calculate stitched enhancers' read density signal
	- `ROSE_callSuper.R`: Rank regions by their density signal and create cutoff to separate super-enhancers from typical enhancers

Example ROSE data is proved by Young lab and can be downloaded from: <https://shorturl.at/nouzR>

---

### 3. USAGE

Program is run by calling `ROSE_main.py`

From within root directory: 
`python ROSE_main.py -g GENOME_BUILD -i INPUT_CONSTITUENT_GFF -r RANKING_BAM -o OUTPUT_DIRECTORY`
`[optional: -s STITCHING_DISTANCE -t TSS_EXCLUSION_ZONE_SIZE -c CONTROL_BAM]`

Required parameters:

`GENOME_BUILD`: one of hg18, hg19, hg38, mm8, mm9, or mm10 referring to the UCSC genome build used for read mapping

`INPUT_CONSTITUENT_GFF`: .gff file (described above) of regions that were previously calculated to be enhancers. i.e. Med1-enriched regions identified using MACS.

`RANKING_BAM`: .bam file to be used for ranking enhancers by density of this factor. i.e. Med1 ChIP-Seq reads.
OUTPUT_DIRECTORY: directory to be used for storing output.

Optional parameters:

`STITCHING_DISTANCE`: maximum distance between two regions that will be stitched together (Default: 12.5kb)

`TSS_EXCLUSION_ZONE_SIZE`: exclude regions contained within +/- this distance from TSS in order to account for promoter biases (Default: 0; recommended if used: 2500). If this value is 0, will not look for a gene file.

`CONTROL_BAM`: .bam file to be used as a control. Subtracted from the density of the `RANKING_BAM`. i.e. Whole cell extract reads.

### 4. CODE PROCEDURE

`ROSE_main.py` will:

- format output directory hierarchy
Root name of input .gff (`[input_enhancer_list].gff`) used as naming root for output files.
- stitch enhancer constituents in `INPUT_CONSTITUENT_GFF` based on `STITCHING_DISTANCE` and make .gff and .bed of stitched collection. TSS exclusion, if not zero, is attempted before stitching. Names of stitched regions start with number of regions stitched followed by leftmost constituent ID
- call `bamToGFF.py` to get density of `RANKING_BAM` and `CONTROL_BAM` in stitched regions and constituents.
Maximum time to wait for `bamToGFF.py` is 12h but can be changed -- quits if running too long.
- call `callSuper.R` to sort stitched enhancers by their background-subtracted density of `RANKING_BAM` and separate into two groups

### 5. OUTPUT:

All file names begin with the root of `INPUT_CONSTITUENT_GFF`

`**OUTPUT_DIRECTORY/gff/`

.gff: copied .gff file of `INPUT_CONSTITUENT_GFF` 
(chrom, name, [blank], start, end, [blank], [blank], strand, [blank], [blank], name)
`STITCHED.gff`: regions created by stitching together `INPUT_CONSTITUENT_GFF` at `STITCHING_DISTANCE`
(chrom, name, [blank], start, end, [blank], [blank], strand, [blank], [blank], name) Name is number of constituents stitched together followed by ID of leftmost constituent.

`**OUTPUT_DIRECTORY/mappedGFF/`

`*_MAPPED.gff`: output of bamToGFF using each bam file containing densities of factor in each constituent
(constituent ID, region tested, average read density in units of reads-per-million-mapped per bp of constituent)

`*_STITCHED*_MAPPED.gff`: output of bamToGFF using each bam file containing densities of factor in each stitched enhancer (stitched enhancer ID, region tested, average read density in units of reads-per-million-mapped per bp of stitched enhancer)

`**OUTPUT_DIRECTORY/`

`STITCHED_ENHANCER_REGION_MAP.txt`: all densities from bamToGFF calculated in stitched enhancers 
(stitched enhancer ID, chromosome, stitched enhancer start, stitched enhancer end, number of constituents stitched, rank of `RANKING_BAM` signal, signal of `RANKING_BAM`). Signal of `RANKING_BAM` is density times length

`*_AllEnhancers.table.txt`: Rankings and super status for each stitched enhancer
(stitched enhancer ID, chromosome, stitched enhancer start, stitched enhancer end, number of constituents stitched, size of constituents that were stitched together, signal of `RANKING_BAM`, rank of `RANKING_BAM`, binary of super-enhancer (1) vs. typical (0)). Signal of `RANKING_BAM` is density times length.

`*_SuperEnhancers.table.txt`: Rankings and super status for super-enhancers 
(stitched enhancer ID, chromosome, stitched enhancer start, stitched enhancer end, number of constituents stitched, size of constituents that were stitched together, signal of `RANKING_BAM`, rank of `RANKING_BAM`, binary of super-enhancer (1) vs. typical (0)). Signal of `RANKING_BAM` is density times length.

`*_Enhancers_withSuper.bed`: .bed file to be loaded into the UCSC browser to visualize super-enhancers and typical enhancers.
(chromosome, stitched enhancer start, stitched enhancer end, stitched enhancer ID, rank by `RANKING_BAM` signal)

`*_Plot_points.png`: visualization of the ranks of super-enhancers and the two groups. Stitched enhancers are ranked by their `RANKING_BAM` signal and their ranks are along the X axis. Corresponding `RANKING_BAM` signal on the Y axis.

---

Developed using `Python 3.8.10`, `R 4.2.1`, and `SAMtools 1.10`

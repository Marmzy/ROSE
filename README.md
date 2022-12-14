## ROSE

### Description

Rank Ordering of Super-Enhancers aka ROSE is a tool for identifying super-enhancers. It does this by separating super-enhancers from typical enhancers using sequencing data (.bam) given a file of previsously identified constituent enhancers (.gff). The original ROSE tool was developed by Charles Y. Lin, David A. Orlando and Brian J. Abraham at Young Lab Whitehead Institute/MIT. This new ROSE version is an attempt to update the code from Python 2 to 3, convert the few R code to Python, use newer versions of tools, make the code more readable to allow for better in-depth understanding of the algorithm and to increase the computational speed.

This version of ROSE was developed using `Python 3.8.10`, and `SAMtools 1.10`

---

### Citation

*Master Transcription Factors and Mediator Establish Super-Enhancers at Key Cell Identity Genes*
Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young. [Cell](https://www.sciencedirect.com/science/article/pii/S0092867413003929) 153, 307-319, April 11, 2013

and

*Selective Inhibition of Tumor Oncogenes by Disruption of Super-enhancers* 
Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando, Christopher R. Vakoc, James E. Bradner, Tong Ihn Lee, and Richard A. Young. [Cell](https://www.sciencedirect.com/science/article/pii/S0092867413003930) 153, 320-334, April 11, 2013

---

### Requirements

- .bam file of sequencing reads for factor of interest (reads for control is recommended, but optional).
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
 -a, --annot      UCSC table track annotation file

 #Additional arguments
 -v, --verbose    Print verbose messages (default=true)

 #Additional arguments for ROSE_main.py
 -c, --control    .bam file to rank enhancers by
 -s, --stitch     Max linking distance for stitching (default=12500)
 -t, --tss        Distance from TSS to exclude (0 = no TSS exclusion) (default=0)
 -d, --debug      Enhancer stitching debugging output (default=False)

 #Additional arguments for ROSE_bamToGFF.py
 -n, --sense      Strand to map to (default='both')
 -f, --floor      Read floor threshold necessary to count towards density (default=1)
 -x, --extension  Extends reads by n bp (default=200)
 -p, --rpm        Normalizes density to reads per million (rpm) (default=true)
 -m, --matrix     Variable bin sized matrix (default=1)

Example: ROSE.sh -g hg18 -i ./data/HG18_MM1S_MED1.gff -o output -r ./data/MM1S_MED1.hg18.bwt.sorted.bam -a ./data/annotation/hg18_refseq.ucsc -c ./data/MM1S_WCE.hg18.bwt.sorted.bam -s 12500 -t 2500 -n both -x 200 -p true -m 1 -v true
```

`ROSE.sh` will run ROSE from start to end. It will (amongst others) call the following scripts:

- `ROSE_main.py`: Stitches regions together to form stitched enhancers
- `ROSE_bamToGFF.py`: Map .bam reads to stitched enhancers and calculate read density
- `ROSE_mapCollection.py`: Calculate stitched enhancers' read density signal
- `ROSE_callSuper.R`: Rank regions by their density signal and create cutoff to separate super-enhancers from typical enhancers

Example ROSE data is provided by Young lab and can be downloaded from: <https://shorturl.at/nouzR>

---

### Output

Explanation of the ROSE output files

`~/output_dir/gff/`

- `*.gff3`: .gff3 formatted copy of --input file
- `*_stitched.gff3`: .gff3 file of stitched enhancer loci. Names reflect the number of enhancers stitched together and contain leftmost enhancer ID.
- `*_stitched.debug`: List of enhancer loci that have not been stitched together and their reasons why

`~/output_dir/mappedGFF/`

- `*_mapped.txt`: Density of reads for each enhancer
- `*_stitched_*_mapped.txt`: Density of reads for each stitched enhancer locus

`~/output_dir/`

- `*_stitched_*_enhancer_region_map.txt`: Density signal for all stitched enhancer loci
- `*_AllEnhancers.table.txt`: Ranking by (control corrected) read density signal and classification of all stitched enhancer loci
- `*_SuperEnhancers.table.txt`: Ranking by (control corrected) read density signal of super-enhancer stitched enhancer loci only
- `*_Enhancers_withSuper.bed`: .bed file to be loaded into the UCSC browser to visualize super-enhancers and typical enhancers
- `*_Plot_points.png`: Visualisation of the stitched enhancer loci read density signals

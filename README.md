## ROSE

### Description

Rank Ordering of Super-Enhancers aka ROSE is a tool for identifying super-enhancers. It does this by separating super-enhancers from typical enhancers using sequencing data (.bam) given a file of previsously identified constituent enhancers (.gff). The original ROSE tool was developed by Charles Y. Lin, David A. Orlando and Brian J. Abraham at Young Lab Whitehead Institute/MIT. This new ROSE version is an attempt to update the code from Python 2 to 3, convert the few R code to Python, use newer versions of tools, make the code more readable to allow for better in-depth understanding of the algorithm and to increase the computational speed.

This version of ROSE was developed using `Python 3.10.13`, and `SAMtools 1.16`

---

### Citation

*Master Transcription Factors and Mediator Establish Super-Enhancers at Key Cell Identity Genes*
Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young. [Cell](https://www.sciencedirect.com/science/article/pii/S0092867413003929) 153, 307-319, April 11, 2013

and

*Selective Inhibition of Tumor Oncogenes by Disruption of Super-enhancers* 
Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando, Christopher R. Vakoc, James E. Bradner, Tong Ihn Lee, and Richard A. Young. [Cell](https://www.sciencedirect.com/science/article/pii/S0092867413003930) 153, 320-334, April 11, 2013

---

### Requirements

- .bam file(s) of sequencing reads for factor(s) of interest (reads for control are recommended, but optional).
	- .bam file(s) must have chromosome IDs starting with "chr"
	- .bam file(s) must be sorted and indexed using [SAMtools](http://www.htslib.org/doc/samtools.html)
	- If a control .bam file is used, the target .bam file(s) must be located in a different directory

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

- .yaml format configuration file 

---

### Usage

As ROSE contains a lot of optional parameters, rather than passing them as additional options on the command line, the newest version opts to use .yaml configuration files. The .yaml file allows for easier tracking of parameters used during an analysis. Below is an example of a .yaml file to be used on the example data provided by Young lab. The exmaple .yaml file uses default setting recommended by Young lab.

```
data:
  annotation: "./data/annotation/hg18_refseq.ucsc"	# UCSC table track annotation file
  control: "./data/MM1S_WCE.hg18.bwt.sorted.bam"	  # Control .bam file to rank enhancers by
  input: "./data/HG18_MM1S_MED1.gff"				        # File (.bed, .gff or .gtf) containing enhancer binding sites
  output: "output"									                # Output directory name
  rankby: "./data/bams"								              # List of .bam files to rank enhancers by

mapping:
  extension: 200									                  # Extends reads by n bp
  floor: 1											                    # Read floor threshold necessary to count towards density
  matrix: 1											                    # Variable bin sized matrix
  rpm: True											                    # Normalizes density to reads per million (rpm)
  sense: "both"										                  # Strand to map to

stitching:
  debug: True										                    # Enhancer stitching debugging output
  stitch: 12500										                  # Max linking distance for stitching
  tss: 2500											                    # Distance from TSS to exclude

verbose: True										                    # Print verbose messages
```

Once the .yaml file has been set, ROSE can easily be run using the following command: `python3 ROSE.py -c ./config/example.yaml`

`ROSE.py` will run ROSE from start to end. It will (amongst others) call the following scripts:

- `ROSE_main.py`: Stitches regions together to form stitched enhancers
- `ROSE_bamToGFF.py`: Map .bam reads to stitched enhancers and calculate read density
- `ROSE_mapCollection.py`: Calculate stitched enhancers' read density signal
- `ROSE_callSuper.py`: Rank regions by their density signal and create cutoff to separate super-enhancers from typical enhancers

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
- `*_plot_points.png`: Visualisation of the stitched enhancer loci read density signals

---

### Docker

[nottuh/rose](https://hub.docker.com/r/nottuh/rose) is the official Docker image for this version ROSE. It contains all necessary tools and packages to run ROSE smoothly. After downloading the latest image from Docker Hub, ROSE can easily be run by mounting the directory containing the input data to the image. Below is a snippet of how to run ROSE with the Docker image, using the example data provided by Young lab:

`docker run --volume $PWD:$PWD --workdir $PWD nottuh/rose:1.3.0 python3 ROSE.py -c example.yaml`

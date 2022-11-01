#!/usr/bin/env Rscript

suppressMessages(library(data.table))
suppressMessages(library(docstring))
suppressMessages(library(optparse))


get_path <- function() {
	#' Get path to ROSE tool directory
	#' 
	#' @return Path to ROSE tool directory

	path <- strsplit(getwd(), "/")
	return(paste(unlist(path)[1:which(unlist(path)=="ROSE")], collapse="/"))
}


option_list <- list(
	make_option(c("-o", "--output"), action="store", type="character", default=NULL, help="Output directory name"),
	make_option(c("-d", "--density"), action="store", type="character", default=NULL, help="Stitched enhancer loci signal density file"),
	make_option(c("-g", "--gff"), action="store", type="character", default=NULL, help="Constituent enhancers peak file"),
	make_option(c("-c", "--control"), action="store", type="character", default=NULL, help="Control .bam file")
)

#Parse arguments from command line
opt <- parse_args(OptionParser(option_list=option_list))

#Initialising variables
path <- get_path()
densityName <- tail(unlist(strsplit(opt$density, "/")))

#Load necessary functions
source(paste(c(path, "src/superenhancers/output.R"), collapse="/"))
source(paste(c(path, "src/superenhancers/super_enhancer.R"), collapse="/"))

#Read stitched enhancer loci density signal file as dataframe
stitched_regions <- fread(opt$density, sep="\t")

#Subtract background signal if control is available
if (is.null(opt$control)) {
	rankBy_vector = stitched_regions[,7]
} else {
	rankBy_vector = stitched_regions[,7] - stitched_regions[,8]
}

#Setting negative values to 0
rankBy_vector[rankBy_vector < 0] <- 0

#Calculating the superenhancer signal cut-off value
cutoff_options <- calculate_cutoff(as.matrix(rankBy_vector))

#Get superenhancers based on cut-off
superEnhancerRows <- which(rankBy_vector > cutoff_options$absolute)
typicalEnhancers <- setdiff(1:nrow(stitched_regions), superEnhancerRows)
enhancerDescription <- paste(basename(opt$gff)," Enhancers\nCreated from ", tail(unlist(strsplit(opt$density, "/")), n=1), "\nRanked by ", colnames(stitched_regions)[7], "\nUsing cutoff of ", cutoff_options$absolute, " for Super-Enhancers", sep="", collapse="")

#Get base .gff file name
gff <- tools::file_path_sans_ext(basename(opt$gff))

#Creating hockey stick plot
hockey_stick_plot(rankBy_vector, cutoff_options, superEnhancerRows, opt$output, gff, opt$control, colnames(stitched_regions)[7])

#Rank stitched enhancer loci by background signal corrected density signal and output to .bed file
bedFileName = paste(opt$output, "/", gff, "_Enhancers_withSuper.bed", sep="")
convert_stitched_to_bed(stitched_regions, paste(colnames(stitched_regions)[7], "Enhancers"), enhancerDescription, bedFileName, densitySignal=as.matrix(rankBy_vector), splitSuper=TRUE, superRows=superEnhancerRows)

#Output stitched enhancer loci ranking
convert_stitched_to_gateway_bed(stitched_regions, paste(opt$output, "/", gff, sep=""), splitSuper=TRUE, densitySignal=as.matrix(rankBy_vector), superRows=superEnhancerRows)

#Select super enhancers
true_super_enhancers <- stitched_regions[superEnhancerRows,]

#Create matrix of stitched enhancer loci rankings and super status
data = c(nrow(stitched_regions)-rank(rankBy_vector, ties.method="first")+1, ifelse(1:nrow(stitched_regions) %in% superEnhancerRows, 1, 0))
additionalTableData <- matrix(data=data, ncol=2, nrow=nrow(stitched_regions))
colnames(additionalTableData) <- c("enhancerRank", "isSuper")

#Output rankings and status matrix
enhancerTableFile = paste(opt$output, "/", gff, "_AllEnhancers.table.txt", sep="")
superTableFile = paste(opt$output, "/", gff, "_SuperEnhancers.table.txt", sep="")

writeSuperEnhancer_table(stitched_regions, enhancerDescription, enhancerTableFile, additionalData=additionalTableData)
writeSuperEnhancer_table(true_super_enhancers, enhancerDescription, superTableFile, additionalData=additionalTableData[superEnhancerRows,])








	

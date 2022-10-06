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


convert_stitched_to_bed <- function(inputStitched,trackName,trackDescription,outputFile,splitSuper=TRUE,score=c(),superRows=c(),baseColor="0,0,0",superColor="255,0,0"){
	outMatrix <- matrix(data="",ncol=4+ifelse(length(score)==nrow(inputStitched),1,0),nrow=nrow(inputStitched))
	
	outMatrix[,1] <- as.character(inputStitched$CHROM)
	outMatrix[,2] <- as.character(inputStitched$START)
	outMatrix[,3] <- as.character(inputStitched$STOP)
	outMatrix[,4] <- as.character(inputStitched$REGION_ID)
	if(length(score)==nrow(inputStitched)){
		score <- rank(score,ties.method="first")
		score <- length(score)-score+1  #Stupid rank only does smallest to largest. 
		outMatrix[,5] <- as.character(score)
	}
	trackDescription <- paste(trackDescription,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
	trackDescription <- gsub("\n","\t", trackDescription)
	tName <- gsub(" ","_",trackName)
	cat('track name="', tName,'" description="', trackDescription,'" itemRGB=On color=',baseColor,"\n",sep="",file=outputFile)
	write.table(file= outputFile,outMatrix,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	if(splitSuper==TRUE){
		cat("\ntrack name=\"Super_", tName,'" description="Super ', trackDescription,'" itemRGB=On color=', superColor,"\n",sep="",file=outputFile,append=TRUE)
		write.table(file= outputFile,outMatrix[superRows,],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	}
}



convert_stitched_to_gateway_bed <- function(inputStitched,outputFileRoot,splitSuper=TRUE,score=c(),superRows=c()){
	outMatrix <- matrix(data="",ncol=6,nrow=nrow(inputStitched))
	
	outMatrix[,1] <- as.character(inputStitched$CHROM)
	outMatrix[,2] <- as.character(inputStitched$START)
	outMatrix[,3] <- as.character(inputStitched$STOP)
	outMatrix[,4] <- as.character(inputStitched$REGION_ID)
	
	if(length(score)==nrow(inputStitched)){
		score <- rank(score,ties.method="first")
		score <- length(score)-score+1  #Stupid rank only does smallest to largest. 
		outMatrix[,5] <- as.character(score)
	}
	
	outMatrix[,6] <- as.character(rep('.',nrow(outMatrix)))
	
	
	outputFile1 = paste(outputFileRoot,'_Gateway_Enhancers.bed',sep='')
	write.table(file= outputFile1,outMatrix,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	if(splitSuper==TRUE){
		outputFile2 = paste(outputFileRoot,'_Gateway_SuperEnhancers.bed',sep='')

		write.table(file= outputFile2,outMatrix[superRows,],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	}
}




writeSuperEnhancer_table <- function(superEnhancer,description,outputFile,additionalData=NA){
	description <- paste("#",description,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
	description <- gsub("\n","\n#",description)
	cat(description,"\n",file=outputFile)
	if(is.matrix(additionalData)){
		if(nrow(additionalData)!=nrow(superEnhancer)){
			warning("Additional data does not have the same number of rows as the number of super enhancers.\n--->>> ADDITIONAL DATA NOT INCLUDED <<<---\n")
		}else{
			superEnhancer <- cbind(superEnhancer,additionalData)
			superEnhancer = superEnhancer[order(superEnhancer$enhancerRank),]
			
		}
	}
	write.table(file=outputFile,superEnhancer,sep="\t",quote=FALSE,row.names=FALSE,append=TRUE)
}



#============================================================================
#===================SUPER-ENHANCER CALLING AND PLOTTING======================
#============================================================================

#Parse arguments from command line
option_list <- list(
	make_option(c("-o", "--output"), action="store", type="character", default=NULL, help="Output directory name"),
	make_option(c("-d", "--density"), action="store", type="character", default=NULL, help="Stitched enhancer loci signal density file"),
	make_option(c("-g", "--gff"), action="store", type="character", default=NULL, help=""),
	make_option(c("-c", "--control"), action="store", type="character", default=NULL, help="Control .bam file")
)

opt <- parse_args(OptionParser(option_list=option_list))

#Initialising variables
path <- get_path()
densityName <- tail(unlist(strsplit(opt$density, "/")))

# outFolder <- opt$output
# enhancerFile <- opt$density
# enhancerName <- opt$gff
# wceName <- opt$control

#Load necessary functions
source(paste(c(path, "src/utils/output.R"), collapse="/"))
source(paste(c(path, "src/utils/super_enhancer.R"), collapse="/"))

#Read stitched enhancer loci density signal file as dataframe
stitched_regions <- fread(opt$density, sep="\t")

# rankBy_factor = colnames(stitched_regions)[7]
# prefix = unlist(strsplit(rankBy_factor, "_"))[1]

# print(rankBy_factor)
# print(prefix)

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
enhancerDescription <- paste(opt$gff," Enhancers\nCreated from ", tail(unlist(strsplit(opt$density, "/")), n=1), "\nRanked by ", colnames(stitched_regions)[7], "\nUsing cutoff of ", cutoff_options$absolute, " for Super-Enhancers", sep="", collapse="")

#Creating hockey stick plot
# hockey_stick_plot(opt$output, opt$gff)
quit()




#Writing a bed file
bedFileName = paste(outFolder,enhancerName,'_Enhancers_withSuper.bed',sep='')
convert_stitched_to_bed(stitched_regions,paste(rankBy_factor,"Enhancers"), enhancerDescription,bedFileName,score=rankBy_vector,splitSuper=TRUE,superRows= superEnhancerRows,baseColor="0,0,0",superColor="255,0,0")


bedFileRoot = paste(outFolder,enhancerName,sep='')

convert_stitched_to_gateway_bed(stitched_regions,bedFileRoot,splitSuper=TRUE,score=rankBy_vector,superRows = superEnhancerRows)


#This matrix is just the super_enhancers
true_super_enhancers <- stitched_regions[superEnhancerRows,]

additionalTableData <- matrix(data=NA,ncol=2,nrow=nrow(stitched_regions))
colnames(additionalTableData) <- c("enhancerRank","isSuper")
additionalTableData[,1] <- nrow(stitched_regions)-rank(rankBy_vector,ties.method="first")+1
additionalTableData[,2] <- 0
additionalTableData[superEnhancerRows,2] <- 1


#Writing enhancer and super-enhancer tables with enhancers ranked and super status annotated
enhancerTableFile = paste(outFolder,enhancerName,'_AllEnhancers.table.txt',sep='')
writeSuperEnhancer_table(stitched_regions, enhancerDescription,enhancerTableFile, additionalData= additionalTableData)

superTableFile = paste(outFolder,enhancerName,'_SuperEnhancers.table.txt',sep='')
writeSuperEnhancer_table(true_super_enhancers, enhancerDescription,superTableFile, additionalData= additionalTableData[superEnhancerRows,])








	

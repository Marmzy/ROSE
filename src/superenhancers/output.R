#!/usr/bin/env Rscript

suppressMessages(library(docstring))


convert_stitched_to_bed <- function(stitchedRegions, trackName, trackDescription, output, splitSuper=TRUE, densitySignal=c(), superRows=c()) {
	#' Create a .bed file ranking stitched enhancer loci by corrected density signal
    #' 
    #' @param stitchedRegions list Stitched enhancer loci signal density dataframe
    #' @param trackName character Target .bam file
    #' @param trackDescription character Analysis description
    #' @param output character Output file name
    #' @param splitSuper logical Boolean whether to add a superenhancers only section to the .bed file 
    #' @param densitySignal double Background signal corrected density signal values
    #' @param superRows integer Vector of superenancer indices
    #' 
    #' @return .bed file containing stitched enhancer loci information and their ranking

    #Create output matrix
    outMatrix <- matrix(data="", ncol=4+ifelse(nrow(densitySignal)==nrow(stitchedRegions),1,0), nrow=nrow(stitchedRegions))
	outMatrix[,1] <- as.character(stitchedRegions$CHROM)
	outMatrix[,2] <- as.character(stitchedRegions$START)
	outMatrix[,3] <- as.character(stitchedRegions$STOP)
	outMatrix[,4] <- as.character(stitchedRegions$REGION_ID)

    #Rank density signal-corrected stitched enhancer loci from smallest to largest
	if (nrow(densitySignal) == nrow(stitchedRegions)) {
		densitySignal <- rank(densitySignal, ties.method="first")
		densitySignal <- length(densitySignal)-densitySignal+1
		outMatrix[,5] <- as.character(densitySignal)
	}

    #Create description track for output .bed file
	trackDescription <- paste(trackDescription, "\nCreated on ", format(Sys.time(), "%b %d %Y"), collapse="", sep="")
	trackDescription <- gsub("\n", "\t", trackDescription)
	tName <- gsub(" ", "_", trackName)
	cat('track name="', tName, '" description="', trackDescription, '" itemRGB=On color=0,0,0\n', sep="", file=output)
	fwrite(file=output, outMatrix, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

	if (splitSuper==TRUE) {
		cat("\ntrack name=\"Super_", tName, '" description="Super ', trackDescription, '" itemRGB=On color=255,0,0\n', sep="", file=output, append=TRUE)
		fwrite(file=output, outMatrix[superRows,], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
	}
}


convert_stitched_to_gateway_bed <- function(stitchedRegions, output, splitSuper=TRUE, densitySignal=c(), superRows=c()) {
    #' Create a .bed file ranking stitched enhancer loci by corrected density signal
    #' 
    #' @param stitchedRegions list Stitched enhancer loci signal density dataframe
    #' @param output character Output file name
    #' @param splitSuper logical Boolean whether to add a superenhancers only section to the .bed file 
    #' @param densitySignal double Background signal corrected density signal values
    #' @param superRows integer Vector of superenancer indices
	#' 
    #' @return .bed file containing stitched enhancer loci information and their ranking

    #Create output matrix
    outMatrix <- matrix(data="", ncol=6, nrow=nrow(stitchedRegions))
	outMatrix[,1] <- as.character(stitchedRegions$CHROM)
	outMatrix[,2] <- as.character(stitchedRegions$START)
	outMatrix[,3] <- as.character(stitchedRegions$STOP)
	outMatrix[,4] <- as.character(stitchedRegions$REGION_ID)
	
    #Rank stitched enhancer loci corrected density signal from smallest to largest
	if (nrow(densitySignal) == nrow(stitchedRegions)) {
		densitySignal <- rank(densitySignal, ties.method="first")
		densitySignal <- length(densitySignal)-densitySignal+1
		outMatrix[,5] <- as.character(densitySignal)
	}
	outMatrix[,6] <- as.character(rep(".", nrow(outMatrix)))
	
    #Output matrix to .bed file
	outputFile1 = paste(output, "_Gateway_Enhancers.bed", sep="")
	fwrite(file=outputFile1, outMatrix, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)

	if (splitSuper==TRUE) {
		outputFile2 = paste(output, "_Gateway_SuperEnhancers.bed", sep="")
		fwrite(file=outputFile2, outMatrix[superRows,], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
	}
}


hockey_stick_plot <- function(vector, cutoff_options, super_enhancer_rows, out_folder, density, control, bam) {
    #' Visualise superenhancer and regular enhancer signals
    #' 
    #' @param vector list Background subtracted stitched enhancer loci density signal list
    #' @param cutoff_options list Superenhancer signal cutoff values list
    #' @param super_enhancer_rows integer Superenhancer row integers
    #' @param out_folder character Output folder name
    #' @param density character Stitched enahncer loci signal density file name
    #' @param control character Control .bam file name
    #' @param bam character Target .bam file name
    #' 
    #' @return Hockey stick plot

    #Ordering stitched enhancer loci signal density vector from low to high
    vector <- as.matrix(vector[order(vector)])

    #Creating the hockey stick plot
    png(filename=paste(out_folder, "/", density, "_Plot_points.png", sep=""), height=600, width=600)
    if (is.null(control)) {
        plot(1:length(vector), vector, col="red", xlab=paste(bam, "_enhancers"), ylab=paste(bam, " Signal"), pch=19, cex=2)	
    } else {
        plot(1:length(vector), vector, col="red", xlab=paste(bam, "_enhancers"), ylab=paste(bam, " Signal - ", control), pch=19, cex=2)
    }
    abline(h=cutoff_options$absolute, col="grey", lty=2)
    abline(v=length(vector)-length(super_enhancer_rows), col="grey", lty=2)
    lines(1:length(vector), vector, lwd=4, col="red")
    text(0, 0.8*max(vector), paste("Cutoff used: ", cutoff_options$absolute, "\n", "Super-Enhancers identified: ", length(super_enhancer_rows)), pos=4)
    dev.off()
}


writeSuperEnhancer_table <- function(superEnhancer, description, outputFile, additionalData=NA) {
    #' Output stitched enhancer loci regions data
    #' 
    #' @param superEnhancer list Stitched enhancer loci data table
    #' @param description character Data description
    #' @param outputFile character Output file name
    #' @param additionalData double Matrix of stitched enhancer loci rankings and super status
    #' 
    #' @return Output regions data to tab-delimited file

    #Create descriptive header
	description <- paste("#", description, "\nCreated on ", format(Sys.time(), "%b %d %Y"), collapse="", sep="")
	description <- gsub("\n", "\n#", description)
	cat(description, "\n", file=outputFile)

    #Merge stitched enhancer loci regions data with superenhancer rankings data
	if (is.matrix(additionalData)) {
		if (nrow(additionalData)!=nrow(superEnhancer)) {
			warning("Additional data does not have the same number of rows as the number of super enhancers.\n--->>> ADDITIONAL DATA NOT INCLUDED <<<---\n")
		} else {
			superEnhancer <- cbind(superEnhancer, additionalData)
			superEnhancer = superEnhancer[order(superEnhancer$enhancerRank),]
		}
	}
	fwrite(file=outputFile, superEnhancer, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
}
#!/usr/bin/env Rscript

suppressMessages(library(data.table))
suppressMessages(library(docstring))


 calculate_cutoff <- function(vector) {
    #' Calculate superenhancer signal cut-off by sliding a diagonal line until the tangent is found
    #' 
    #' @param vector matrix Background corrected stitched enhancer signal values matrix 
    #' 
    #' @return List of superenhancer signal cut-off values


	#Get the slope of the line to slide
 	vector <- vector[order(vector)]
	slope <- (max(vector) - min(vector)) / length(vector)

	#Find the x-axis point where a line passing through that point has the minimum number of points below it (ie. tangent)
	xPt <- floor(optimize(numPts_below_line, lower=1, upper=length(vector), vector=vector, slope=slope)$minimum)
	y_cutoff <- vector[xPt]

	return(list(absolute=y_cutoff, overMedian=y_cutoff/median(vector), overMean=y_cutoff/mean(vector)))
}


numPts_below_line <- function(vector, slope, x) {
    #' Calculate linear equation of diagonal at point [x,y] and the number of signal values below the diagonal
    #' 
    #' @param vector matrix Background corrected stitched enhancer signal values matrix 
    #' @param slope double Diagonal slope
    #' @param x double X point at which to calculate the linear equation
    #' 
    #' @return Number of signal values below the diagonal at point [x,y]
    
    #Calculate lienar equation (y = ax+b)
	y <- vector[x]
	b <- y-(slope*x)

    #Calculate number of signal values below diagonal
	xPts <- 1:length(vector)
	return(sum(vector <= (xPts*slope+b)))
}

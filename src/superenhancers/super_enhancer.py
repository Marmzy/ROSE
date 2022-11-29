#!/usr/bin/env python

import math
import numpy as np

from scipy.optimize import minimize_scalar


def calculate_cutoff(
    vector: np.ndarray
) -> np.float:
    """Calculate superenhancer signal density cutoff value

    Args:
        vector (np.ndarray): Array of (control corrected) stitched enhancer loci density signals

    Returns:
        np.float: Density signal cutoff value to delineate superenhancers from normal enhaners 
    """

    #Get the slope of the line to slide
    vector.sort()
    slope = (max(vector) - min(vector)) / len(vector)

    #Minimising the (control corrected) stitched enhancer loci density signal function (aka finding the tangent of the function)
    x_min = math.floor(minimize_scalar(numPts_below_line, bounds=(1, len(vector)), args=(vector, slope), method="bounded")["x"])
    y_cutoff = vector[x_min]

    return y_cutoff


def numPts_below_line(
    x: float,
    vector: np.ndarray,
    slope: np.float64
) -> int:
    """Stitched enhancer loci density signal function to minimise

    Args:
        x (float): Signal density value to perform the calculation at
        vector (np.ndarray): Signal density values vector
        slope (np.float64): Slope coefficient

    Returns:
        int: Number of signal density values equal to or smaller than a line going through point x 
    """

    #Calculate line equation at point x   (y = ax+b)
    y = vector[int(x)-1]
    b = y-(slope*x)

    #Calculate y values of points on the line and the number of points that are larger than those of the vector
    xPts = np.array(range(1, len(vector)+1))
    return sum(vector <= (xPts*slope+b)) 
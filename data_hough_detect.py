
from visualize import visualize
from sim import wave2d
from skimage.transform import hough
import scipy
import numpy as np
from sunpy.time import TimeRange
from sunpy.time import parse_time
import os
import eitwaveutils
import util
import copy
from datetime import timedelta

def main():

    # type of data we want to use
    data_type = '.jp2'

    # where to store those data
    data_storage = "~/Data/eitwave/jp2/AGU/"

    # time range we want
    time_range = TimeRange('2011/06/01','2011/06/02')


    m2deg = 360./(2*3.1415926*6.96e8)
    
    # acquire the data
    hek_result, filelist = eitwaveutils.acquire_data(data_storage, data_type, time_range)

    for flare in hek_result:
        
        # Define the transform parameters
        params = eitwaveutils.params(flare)
        
    
    # load in the data with a single EIT wave
    #filelist = eitwaveutils.loaddata(data_storage, data_type)

    # read in files and accumulate them
    maps = eitwaveutils.accumulate(filelist[0:20], accum=1, super=4, verbose=True)

    # Unravel the maps
    new_maps = eitwaveutils.map_unravel(maps, params, verbose=True)
    
    # calculate the differences
    diffs = eitwaveutils.map_diff(new_maps)

    #determine the threshold to apply to the difference maps.
    #diffs > diff_thresh will be True, otherwise False.
    threshold_maps = eitwaveutils.map_threshold(new_maps,factor=0.7) 

    # transform difference maps into binary maps
    binary_maps = eitwaveutils.map_binary(diffs, threshold_maps)
    
    # detection based on the hough transform
    detection = eitwaveutils.hough_detect(binary_maps, vote_thresh=12)
    
    # detection based on the probabilistic hough transform.  Takes the
    # keywords of the probabilistic hough transform - see the documentation
    # of skimage.transform.probabilistic_hough (scikit-image.org) 
    #detection = eitwaveutils.prob_hough_detect(binary_maps,threshold=10)
    
    
    detection = eitwaveutils.cleanup(detection,
                                     size_thresh=50,
                                     inv_thresh=8)

    #If there is anything left in 'detection', fit a function to the original
    #diffmaps in the region defined by 'detection'. Simplest case: fit a Gaussian
    #in the y-direction for some x or range of x.
    #eitwaveutils.fit_wavefront should probably take the arguments of fitfunc.
    #use 'detection' to guess starting fit parameters?

    #get just the positive elements of the difference map. Perform fitting on these positive diffmaps.
    posdiffs=copy.deepcopy(diffs)
    for i in range(0,len(diffs)):
        temp= diffs[i] < 0
        posdiffs[i][temp] = 0

    #fit a function to the difference maps in the cases where there has been a detection
    wavefront = eitwaveutils.fit_wavefront(posdiffs, detection)
    
    #strip out the velocity information from the wavefront fitting
    velocity = eitwaveutils.wavefront_velocity(wavefront[0])

    #strip out the position and width information from the wavefront fitting
    pos_width = eitwaveutils.wavefront_position_and_width(wavefront[0])
    
    visualize(detection)
    return maps, new_maps, diffs, threshold_maps, binary_maps, detection, wavefront, velocity, pos_width

if __name__ == '__main__':
    main()


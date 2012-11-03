
from visualize import visualize
from sim import wave2d
from skimage.transform import hough
import scipy
import numpy as np
import sunpy
import os
import eitwaveutils



def main():

    m2deg = 360./(2*3.1415926*6.96e8)
    params = {
              "epi_lat": 40., #degrees, HG latitude of wave epicenter
              "epi_lon": -20., #degrees, HG longitude of wave epicenter
              #HG grid, probably would only want to change the bin sizes
              "lat_min": -90.,
              "lat_max": 90.,
              "lat_bin": 0.2,
              "lon_min": -180.,
              "lon_max": 180.,
              "lon_bin": 5.,
              #    #HPC grid, probably would only want to change the bin sizes
              "hpcx_min": -1025.,
              "hpcx_max": 1023.,
              "hpcx_bin": 2.,
              "hpcy_min": -1025.,
              "hpcy_max": 1023.,
              "hpcy_bin": 2.,
              "hglt_obs": 0,
              "rotation": 360./(27.*86400.), #degrees/s, rigid solar rotation
              }
    
    #params = {
    #    "cadence": 12., #seconds
    #    
    #    "hglt_obs": 0., #degrees
    #    "rotation": 360./(27.*86400.), #degrees/s, rigid solar rotation
    #   
    #    #Wave parameters that are initial conditions
    #    "direction": 25., #degrees, measured CCW from HG +latitude
    #    "epi_lat": 30., #degrees, HG latitude of wave epicenter
    #    "epi_lon": 45., #degrees, HG longitude of wave epicenter
    #    
    #    #Wave parameters that can evolve over time
    #    #The first element is constant in time
    #    #The second element (if present) is linear in time
    #    #The third element (if present) is quadratic in time
    #    #Be very careful of non-physical behavior
    #    "width": [90., 1.5], #degrees, full angle in azimuth, centered at 'direction'
    #    "wave_thickness": [6.0e6*m2deg,6.0e4*m2deg], #degrees, sigma of Gaussian profile in longitudinal direction
    #    "wave_normalization": [1.], #integrated value of the 1D Gaussian profile
    #    "speed": [9.33e5*m2deg, -1.495e3*m2deg], #degrees/s, make sure that wave propagates all the way to lat_min for polynomial speed
    #    
    #    #Noise parameters
    #    "noise_type": "Poisson", #can be None, "Normal", or "Poisson"
    #    "noise_scale": 0.3,
    #    "noise_mean": 1.,
    #    "noise_sdev": 1.,
    #    
    #    "max_steps": 20,
    #    
    #    #HG grid, probably would only want to change the bin sizes
    #    "lat_min": -90.,
    #    "lat_max": 90.,
    #    "lat_bin": 0.2,
    #    "lon_min": -180.,
    #    "lon_max": 180.,
    #    "lon_bin": 5.,
    #    
    #    #HPC grid, probably would only want to change the bin sizes
    #    "hpcx_min": -1025.,
    #    "hpcx_max": 1023.,
    #    "hpcx_bin": 2.,
    #    "hpcy_min": -1025.,
    #    "hpcy_max": 1023.,
    #    "hpcy_bin": 2.
    #}
    
    # Lots of big images.  Need to be smart about how to handle the data
    
    # load in the data with a single EIT wave
    filelist = eitwaveutils.loaddata("~/Data/eitwave_data/jp2/20110601_02_04/",
                                     '.jp2')

    # read in files and accumulate them
    maps = eitwaveutils.accumulate(filelist, accum=2, super=4)

    # Unravel the maps
    new_maps = eitwaveutils.map_unravel(maps, params)
    
    # calculate the differences
    diffs = eitwaveutils.map_diff(new_maps, diff_thresh=100)
    
    detection = eitwaveutils.hough_detect(diffs,
                                          vote_thresh=12,
                                          inv_thresh=8)
    
    detection = eitwaveutils.cleanup(detection,
                                     size_thresh=50)
    

    # shape of the data
    imgShape = new_maps[0].shape
    
    # storage for the detection
    detection = []
    
    invThresh = 8
    
    # areas that have a size lower than threshold are exclused
    sizeThresh = 50
    
    
    visualize(detection)
    return maps, diffs, detection

if __name__ == '__main__':
    main()


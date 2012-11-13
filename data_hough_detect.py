
from visualize import visualize
from sim import wave2d
from skimage.transform import hough
import scipy
import numpy as np
from sunpy.net import hek
from sunpy.net import helioviewer
from sunpy.time import TimeRange
from sunpy.time import parse_time
import os
import eitwaveutils
from datetime import timedelta

def main():

    # Time range we are interested in
    time_range = TimeRange('2011/06/01','2011/06/02')
    # Query the HEK for flare information we need
    client = hek.HEKClient()
    hek_result = client.query(hek.attrs.Time('2011/06/01','2011/06/02'),
                              hek.attrs.EventType('FL'),
                              hek.attrs.FRM.Name=='SEC standard')
    
    #vals = eitwaveutils.goescls2number( [hek['fl_goescls'] for hek in hek_result] )
    #flare_strength_index = sorted(range(len(vals)), key=vals.__getitem__)

    # Download all the JP2 files for the duration of the event
    hv = helioviewer.HelioviewerClient()
    for flare in hek_result:
        start_time = parse_time(flare['event_starttime'])
        end_time = start_time + timedelta(minutes=60)
        jp2_list = []
        this_time = start_time
        while this_time <= end_time:
            jp2 = hv.download_jp2(this_time, observatory='SDO', 
                                  instrument='AIA', detector='AIA',
                                  measurement='193',
                                  directory = '~/Data/eitwave/jp2/AGU/')
            if not(jp2 in jp2_list):
                jp2_list.append(jp2)
                
            this_time = this_time + timedelta(seconds = 6)


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
    maps = eitwaveutils.accumulate(filelist[0:10], accum=2, super=4, verbose=True)

    # Unravel the maps
    new_maps = eitwaveutils.map_unravel(maps, params, verbose=True)
    
    # calculate the differences
    diffs = eitwaveutils.map_diff(new_maps, diff_thresh=100)
    
    # detection based on the hough transform
    #detection = eitwaveutils.hough_detect(diffs, vote_thresh=12)
    
    # detection based on the probabilistic hough transform.  Takes the
    # keywords of the probabilistic hough transform - see the documentation
    # of skimage.transform.probabilistic_hough (scikit-image.org) 
    detection = eitwaveutils.prob_hough_detect(diffs)
    
    
    detection = eitwaveutils.cleanup(detection,
                                     size_thresh=50,
                                     inv_thresh=8)
    

    
    
    visualize(detection)
    return maps, diffs, detection

if __name__ == '__main__':
    main()


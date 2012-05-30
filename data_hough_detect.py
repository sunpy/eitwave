
from visualize import visualize
from sim import wave2d
from skimage.transform import hough
import numpy as np
import sunpy
import os

def htLine(distance,angle,img):
    shape = img.shape
    ny = shape[0]
    nx = shape[1]
    eps = 1.0/float(ny)

    if abs(np.sin(angle)) > eps:
        gradient = - np.cos(angle) / np.sin(angle)
        constant = distance / np.sin(angle)
        for x in range(0,nx):
            y = gradient*x + constant
            if y <= ny-1 and y >= 0:
                img[y,x] = 255
            else:
                img[:,distance] = 255

    return img

def main():

    m2deg = 360./(2*3.1415926*6.96e8)
    params = {
              "epi_lat": 30., #degrees, HG latitude of wave epicenter
              "epi_lon": 45., #degrees, HG longitude of wave epicenter
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
    directory = "/home/ireland/eitwave_data/jp2/20110601_02_04/"
    extension = ".jp2"
    # get the file list
    lll = []
    for f in os.listdir(directory):
        if f.endswith(extension):
            lll.append(f)
    
    # accumulate every "naccum" files
    naccum = 2
    
    # accumulate using superpixels
    nsuper = 4
    
    # number of files to consider
    #nfiles = len(lll)
    
    nfiles = 10
    
    # storage for all the maps
    maps = []
    
    # accumulate naccum files, super-pixel them, then add them together.
    j = 0
    while j+naccum <= nfiles:
        i = 0
        print('\n Starting new accumulation:')
        while i < naccum:
            filename = os.path.join(directory,lll[j+i])
            print('  File %(#)i out of %(nfiles)i' % {'#':i+j, 'nfiles':nfiles})
            print('  Reading in file '+filename)
            map1 = (sunpy.make_map(filename)).superpixel((nsuper,nsuper))
            if i == 0:
                m = map1
            else:
                m = m + map1
            i = i + 1
        j = j + naccum
        maps.append(m)
    print 1
    z = wave2d.transform(params,maps[0])
    
    # number of running differences
    ndiff = len(maps)-1
    
    # Each JP2 file has a maximum of 255
    maxval = 255 * nsuper * nsuper * naccum
    
    # difference threshold as a function of the maximum value
    diffthresh = 0.01 * maxval #300
    
    # Hough transform voting threshold
    votethresh = 5000
    
    # shape of the data
    imgShape = maps[0].shape
    
    # storage for the detection
    detection = []
    diffs = []
     
    # calculate running difference images
    for i in range(0,ndiff):
        
        # take the difference
        diffmap = ( maps[i+1] -maps[i] )  > diffthresh
    
        # keep
        diffs.append(diffmap)
    
        # extract the image from the storage array
        img = diffmap
    
        # Perform the hough transform on each of the difference maps
        transform,theta,d = hough(img)
    
        # Filter the hough transform results and find the best lines
        # in the data
        indices =  (transform >votethresh).nonzero()
        distances = d[indices[0]]
        theta = theta[indices[1]]
        n =len(indices[1])
        print n
    
        # Perform the inverse transform to get a series of rectangular
        # images that show where the wavefront is.
        invTransform = sunpy.make_map(maps[i+1])
        invTransform.data = np.zeros(imgShape)
        for i in range(0,n):
            nextLine = htLine( distances[i],theta[i], np.zeros(shape=imgShape) )
            invTransform = invTransform + nextLine
    
        # Dump the inverse transform back into a series of maps
        detection.append(invTransform)
    
    
    visualize(detection)

if __name__ == '__main__':
    main()


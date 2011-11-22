from sim import wave2d
from visualize import visualize
from scikits.image.transform import hough
from scikits.image.morphology import greyscale_dilate
import numpy as np
import pylab as plt

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

m2deg = 360./(2*3.1415926*6.96e8)

params = {
    "cadence": 12., #seconds
    
    "hglt_obs": 0., #degrees
    "rotation": 360./(27.*86400.), #degrees/s, rigid solar rotation
    
    #Wave parameters that are initial conditions
    "direction": 25., #degrees, measured CCW from HG +latitude
    "epi_lat": 30., #degrees, HG latitude of wave epicenter
    "epi_lon": 45., #degrees, HG longitude of wave epicenter
    
    #Wave parameters that can evolve over time
    #The first element is constant in time
    #The second element (if present) is linear in time
    #The third element (if present) is quadratic in time
    #Be very careful of non-physical behavior
    "width": [90., 1.5], #degrees, full angle in azimuth, centered at 'direction'
    "wave_thickness": [6.0e6*m2deg,6.0e4*m2deg], #degrees, sigma of Gaussian profile in longitudinal direction
    "wave_normalization": [1.], #integrated value of the 1D Gaussian profile
    "speed": [9.33e5*m2deg, -1.495e3*m2deg], #degrees/s, make sure that wave propagates all the way to lat_min for polynomial speed
    
    #Noise parameters
    "noise_type": "Poisson", #can be None, "Normal", or "Poisson"
    "noise_scale": 0.3,
    "noise_mean": 1.,
    "noise_sdev": 1.,
    
    "max_steps": 20,
    
    #HG grid, probably would only want to change the bin sizes
    "lat_min": -90.,
    "lat_max": 90.,
    "lat_bin": 0.2,
    "lon_min": -180.,
    "lon_max": 180.,
    "lon_bin": 5.,
    
    #HPC grid, probably would only want to change the bin sizes
    "hpcx_min": -1025.,
    "hpcx_max": 1023.,
    "hpcx_bin": 2.,
    "hpcy_min": -1025.,
    "hpcy_max": 1023.,
    "hpcy_bin": 2.
}

wave_maps = wave2d.simulate_raw(params)
#wave_maps = wave2d.simulate(params, verbose = True)
#visualize(wave_maps)
#
# Use Albert's wavemaps to test the hough transform as a means
# of detecting EIT waves

#
# Initial detection is based on Hough transform of absolute
# value of the running difference.
#
# Possible algorithm outline
#
# (1) Hough transform (HT) of absolute value of the running difference.
# (2) Threshold HT and transform back
# (3) Remove areas which are 'small', keep areas which are large
# (4) Use the remaining detected areas as a mask in the original data
# (5) Apply HT to masked original data
# (6) Threshold and transform back
# (7) Remove areas which are 'small', keep areas which are large
# (8) This localises the wavefront
#

ndiff = len(wave_maps)-1

# difference threshold
diffthresh = 0.2

# Hough transform voting threshold
votethresh = 20

# shape of the data
imgShape = wave_maps[0].shape

# storage for the detection
detection = []

for i in range(0,ndiff):
    diffmap = 255*(abs(wave_maps[i+1] - wave_maps[i]) > diffthresh)
    # extract the image
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
    invTransform = wave_maps[i+1]
    invTransform.data = np.zeros(imgShape)
    for i in range(0,n):
        nextLine = htLine( distances[i],theta[i], np.zeros(shape=imgShape) )
        invTransform = invTransform + nextLine

    # Dump the inverse transform back into a series of maps
    detection.append(invTransform)

visualize(detection)


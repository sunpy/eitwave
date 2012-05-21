from sim import wave2d
from visualize import visualize
from scikits.image.transform import hough
from scikits.image.morphology import greyscale_dilate
import numpy as np
import pylab as plt
import sunpy

def reshaper(img,sx=1,sy=1):
    """ Sum the input image into sx by sy superpixels.  Taken from
    http://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html
    """
    if img.ndim != 2:
        print('Input must have exactly two dimensions; input passed to function has %i dimension(s)', img.ndim)
        return None
    
    sz = np.shape(img)
    if sz[0] % sx != 0:
        print('Sum value "sx" must divide exactly into the image x-dimension size')
        return None
    if sz[1] % sy != 0:
        print('Sum value "sy" must divide exactly into the image y-dimension size')
        return None
   
    # Get the size of the new super-pixel array
    nx = sz[0]/sx
    ny = sz[1]/sy
    # Perform the summing by reshaping up to a higher dimensional array and
    # summing along the higher dimensions
    reshaped = img.reshape(nx,sx,ny,sy)
    return reshaped


def sum(img,sx=1,sy=1):
    """ Sum the input image into sx by sy superpixels.  Taken from
    http://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html
    """
    y = reshaper(img,sx=sx,sy=sy)
    summed = y.sum(axis=3).sum(axis=1)
    return summed

def avg(img,sx=1,sy=1):
    """ Average the input image into sx by sy superpixels.  Taken from
    http://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html
    """
    y = reshaper(img,sx=sx,sy=sy)
    average = (y.sum(axis=3).sum(axis=1))/(np.float32(sx*sy))
    return average


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

# load in all the data

# normalize all the images

# number of scales
nscale = 6

# number of running differences
ndiff = len(maps)-1

# difference threshold
diffthresh = 0.01

# Hough transform voting threshold
votethresh = 10

# shape of the data
imgShape = wave_maps[0].shape

# storage for the detection
detection = []
diffs = []

# Sum in space over multiple lengthscales
for j in range(0,nscale):
    lengthscale = 2^j
    
    # calculate running difference images
    for i in range(0,ndiff):
        
        # get the summed images
        map_after = maps[i+1].sum(lengthscale,lengthscale)
        map_now   = maps[i].sum(lengthscale,lengthscale)
    
        # take the difference
        diffmap = abs(map_after - map_now) > diffthresh

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
        invTransform = sunpy.map.BaseMap(wave_maps[i+1])
        invTransform.data = np.zeros(imgShape)
        for i in range(0,n):
            nextLine = htLine( distances[i],theta[i], np.zeros(shape=imgShape) )
            invTransform = invTransform + nextLine

        
        # Dump the inverse transform back into a series of maps
        detection.append(invTransform.resample([4096,4096]))


visualize(detection)


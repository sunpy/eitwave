#
# Utilities that implement functions to enable EIT wave detection
#
from visualize import visualize
from sim import wave2d
from skimage.transform import hough
import scipy
import numpy as np
import sunpy
import os
import util2


def loaddata(directory, extension):
    """ get the file list and sort it.  For well behaved file names the file
    name list is returned ordered by time"""
    lst = []
    dir = os.path.expanduser(directory)
    for f in os.listdir(dir):
        if f.endswith(extension):
            lst.append(os.path.join(dir,f))
    return sorted(lst)

def accumulate(filelist, accum=2, super=4, verbose=False):
    """Add up data in time and space. Accumulate 'accum' files in time, and 
    then form the images into super by super superpixels."""
    # counter for number of files.
    j = 0
    # storage for the returned maps
    maps = []
    nfiles = len(filelist)
    while j+accum <= nfiles:
        i = 0
        while i < accum:
            filename = filelist[i]
            if verbose:
                print('  File %(#)i out of %(nfiles)i' % {'#':i+j, 'nfiles':nfiles})
                print('  Reading in file '+filename)
            map1 = (sunpy.make_map(filename)).superpixel((super,super))
            if i == 0:
                m = map1
            else:
                m = m + map1
            i = i + 1
        j = j + accum
        maps.append(m)
    return maps

def map_unravel(maps, params, verbose=False):
    """ Unravel the maps into a rectangular image. """
    new_maps =[]
    for m in maps:
        if verbose:
            print("Unraveling map at "+str(m.date))
        unraveled = util2.map_hpc_to_hg_rotate(m,
                                               epi_lon=params.get('epi_lon'),
                                               epi_lat=params.get('epi_lat'),
                                               xbin=5,
                                               ybin=0.2)
        unraveled[np.isnan(unraveled)]=0.0
        new_maps += [unraveled]
    return new_maps

def linesampleindex(a, b, np=1000):
    """ Get the indices in an array along a line"""
    x, y = np.linspace(a[0],b[0],np), np.linspace(a[1],b[1],np)
    xi = x.astype(np.int)
    yi = y.astype(np.int)
    return xi, yi

def map_diff(maps, diff_thresh=100):
    """ calculate running difference images """
    diffs = []
    for i in range(0,len(maps)-1):
        # take the difference
        diffmap = abs(maps[i+1] - maps[i])>diff_thresh
        # keep
        diffs.append(diffmap)
    return diffs

def hough_detect(diffs, vote_thresh=12):
    """ Use the Hough detection method to detect lines in the data.
    With enough lines, you can fill in the wave front."""
    detection=[]
    
    for i in range(0,len(diffs)):

        # extract the image from the storage array
        img = diffs[i]

        # Perform the hough transform on each of the difference maps
        transform, theta, d = hough(img)
    
        # Filter the hough transform results and find the best lines
        # in the data
        indices =  (transform>vote_thresh).nonzero()
        distances = d[indices[0]]
        theta = theta[indices[1]]
        n = len(indices[1])
    
        # Perform the inverse transform to get a series of rectangular
        # images that show where the wavefront is.
        invTransform = sunpy.make_map(np.zeros(img.shape),
                                      diffs[i]._original_header)
        invTransform.data = np.zeros(img.shape)
        
        # Add in all the detected lines
        for i in range(0,n):
            nextLine = htLine(distances[i], theta[i], np.zeros(shape=img.shape))
            invTransform = invTransform + nextLine
            
        detection.append(invTransform)

    return detection

def cleanup(detection, size_thresh=50, inv_thresh=8):
    """Clean up the detection"""
    cleaned=[]
    
    for d in detection:
        d[(d<inv_thresh).nonzero()] = 0.0
        #invTransform[(invTransform>=invThresh).nonzero()] = 1.0
        labeled_array, num_features = scipy.ndimage.measurements.label(d)
        for j in range(1,num_features):
            region = (labeled_array == j).nonzero()
            if np.size( region ) <= size_thresh:
                d[region] = 0
                
        # Dump the inverse transform back into a series of maps
        cleaned.append(d)
 
    return cleaned


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
                img[y,x] = 1
    else:
        img[:,distance] = 1

    return img



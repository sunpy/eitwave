#
# Utilities that implement functions to enable EIT wave detection
#
from visualize import visualize
from sim import wave2d
from skimage.transform import hough
from skimage.transform import probabilistic_hough
import scipy
import numpy as np
import sunpy
import os
import util2
import types

from sunpy.net import hek
from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange
from sunpy.lightcurve import LogicalLightCurve
from sunpy.lightcurve import LightCurve
from matplotlib import pyplot as plt
import datetime
import pickle
import pandas


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
            filename = filelist[i+j]
            if verbose:
                print('File %(#)i out of %(nfiles)i' % {'#':i+j, 'nfiles':nfiles})
                print('Reading in file '+filename)
            map1 = (sunpy.make_map(filename)).superpixel((super,super))
            if i == 0:
                m = map1
            else:
                m = m + map1
            i = i + 1
        j = j + accum
        maps.append(m)
        if verbose:
            print('Accumulated map List has length %(#)i' % {'#':len(maps)} )
    return maps

def map_unravel(maps, params, verbose=False):
    """ Unravel the maps into a rectangular image. """
    new_maps =[]
    for index, m in enumerate(maps):
        if verbose:
            print("Unraveling map %(#)i of %(n)i " % {'#':index+1, 'n':len(maps)})
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
    
    for img in diffs:

        # Perform the hough transform on each of the difference maps
        transform, theta, d = hough(img)
    
        # Filter the hough transform results and find the best lines in the
        # data.  Keep detections that exceed the Hough vote threshold.
        indices =  (transform>vote_thresh).nonzero()
        distances = d[indices[0]]
        theta = theta[indices[1]]
    
        # Perform the inverse transform to get a series of rectangular
        # images that show where the wavefront is.
        # Create a map which is the same as the 
        invTransform = sunpy.make_map(np.zeros(img.shape), img._original_header)
        invTransform.data = np.zeros(img.shape)
        
        # Add up all the detected lines over each other.  The idea behind
        # adding up all the lines on top of each other is that pixels that
        # have larger number of detections are more likely to be in the
        # wavefront.  Note that we are using th Hough transform - which is used
        # to detect lines - to detect and fill in a region.  You might see this
        # as an abuse of the Hough transform!
        for i in range(0,len(indices[1])):
            nextLine = htLine(distances[i], theta[i], np.zeros(shape=img.shape))
            invTransform = invTransform + nextLine

        detection.append(invTransform)

    return detection

def prob_hough_detect(diffs, **ph_kwargs):
    """Use the probabilistic hough transform to detect regions in the data
    that we will flag as being part of the EIT wave front."""
    detection=[]
    for img in diffs:
        lines = probabilistic_hough(img, ph_kwargs)
        if lines is not None:
            for line in lines:
                pass
    return detection


def cleanup(detection, size_thresh=50, inv_thresh=8):
    """Clean up the detection.  The original detection is liable to be quite
    noisy.  There are many different ways of cleaning it up."""
    cleaned=[]
    
    for d in detection:
        # Remove points from the detections that have less than 'inv_thresh'
        # detections
        d[(d<inv_thresh).nonzero()] = 0.0
        
        #
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

def goescls2number(gcls):
    """Convert GOES classes into number to aid size comparison"""
    def calc(gcls):
        powers_of_ten = {'A':1, 'B':10, 'C':100, 'M':1000, 'X':10000}
        power = gcls[0].upper()
        if power in powers_of_ten:
            return powers_of_ten[power] * float(gcls[1:])
        else:
            return None

    if isinstance(gcls, types.StringType):
        return calc(gcls)
    if isinstance(gcls, types.ListType):
        return [calc(x) for x in gcls]


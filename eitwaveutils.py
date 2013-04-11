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
import util
from sunpy.net import helioviewer
from sunpy.time import TimeRange
from sunpy.time import parse_time
from datetime import timedelta, datetime

def params(flare,**kwargs):
    
    m2deg = 360./(2*3.1415926*6.96e8)

    """ Define the parameters we will use for the unraveling of the maps"""
    params = {"epi_lat": flare['event_coord1'], #30., #degrees, HG latitude of wave epicenter
              "epi_lon": flare['event_coord2'], #45., #degrees, HG longitude of wave epicenter
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

    return params


def acquire_data(directory, extension, flare, duration=60, verbose=True):

    # vals = eitwaveutils.goescls2number( [hek['fl_goescls'] for hek in
    # hek_result] )
    # flare_strength_index = sorted(range(len(vals)), key=vals.__getitem__)
    # Get the data for each flare.
    if verbose:
        print('Event start time: ' + flare['event_starttime'])
        print('GOES Class: ' + flare['fl_goescls'])
    data_range = TimeRange(parse_time(flare['event_starttime']),
                           parse_time(flare['event_starttime']) +
                           timedelta(minutes=duration))
    if extension.lower() == '.jp2':
        data = acquire_jp2(directory, data_range)
    if extension.lower() in ('.fits', '.fts'):
        data = []
    # Return the flare list from the HEK and a list of files for each flare in
    # the HEK flare list
    return data


def listdir_fullpath(d):
    dd = os.path.expanduser(d)
    return [os.path.join(dd, f) for f in os.listdir(dd)]


def get_jp2_dict(directory):
    directory_listing = {}
    l = os.listdir(os.path.expanduser(directory))
    for f in l:
        try:
            ymd = f.split('__')[0]
            hmsbit = f.split('__')[1]
            hms = hmsbit.split('_')[0] + '_' + hmsbit.split('_')[1] + '_' + hmsbit.split('_')[2]
            dt = datetime.strptime(ymd + '__' + hms, '%Y_%m_%d__%H_%M_%S')
            directory_listing[dt] = os.path.join(os.path.expanduser(directory), f)
        except:
            print("Non Helioviewer-timestamp formatted file present.")
    return directory_listing


def acquire_jp2(directory, time_range, observatory='SDO', instrument='AIA',
                detector='AIA', measurement='193', verbose=True):
    """Acquire Helioviewer JPEG2000 files between the two specified times"""

    # Create a Helioviewer Client
    hv = helioviewer.HelioviewerClient()

    # Start the search
    jp2_list = []
    this_time = time_range.t1
    while this_time <= time_range.t2:
        # update the directory dictionary with the latest contents
        directory_dict = get_jp2_dict(directory)

        # find what the closest image to the requested time is
        response = hv.get_closest_image(this_time, observatory=observatory,
                              instrument=instrument, detector=detector,
                              measurement=measurement)

        # if this date is not already present, download it
        if not(response["date"] in directory_dict):
            if verbose:
                print('Downloading new file:')
            jp2 = hv.download_jp2(this_time, observatory=observatory,
                              instrument=instrument, detector=detector,
                              measurement=measurement, directory=directory,
                              overwrite=True)
        else:
            # Otherwise, get its location
            jp2 = directory_dict[response["date"]]
        # Only one instance of this file should exist
        if not(jp2 in jp2_list):
            jp2_list.append(jp2)
            if verbose:
                print('Found file ' + jp2 + '. Total found: ' + str(len(jp2_list)))

        # advance the time
        this_time = this_time + timedelta(seconds=6)
    return jp2_list


def loaddata(directory, extension):
    """ get the file list and sort it.  For well behaved file names the file
    name list is returned ordered by time"""
    lst = []
    loc = os.path.expanduser(directory)
    for f in os.listdir(loc):
        if f.endswith(extension):
            lst.append(os.path.join(loc, f))
    return sorted(lst)


def accumulate(filelist, accum=2, nsuper=4, verbose=False):
    """Add up data in time and space. Accumulate 'accum' files in time, and
    then form the images into super by super superpixels."""
    # counter for number of files.
    j = 0
    # storage for the returned maps
    maps = []
    nfiles = len(filelist)
    while j + accum <= nfiles:
        i = 0
        while i < accum:
            filename = filelist[i + j]
            if verbose:
                print('File %(#)i out of %(nfiles)i' % {'#': i+j, 'nfiles':nfiles})
                print('Reading in file '+filename)
            map1 = (sunpy.make_map(filename)).superpixel((nsuper, nsuper))
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
        unraveled = util.map_hpc_to_hg_rotate(m,
                                               epi_lon=params.get('epi_lon'),
                                               epi_lat=params.get('epi_lat'),
                                               xbin=5,
                                               ybin=0.2)
        unraveled[np.isnan(unraveled)]=0.0
        new_maps += [unraveled]
    return new_maps


def check_dims(new_maps):
    """ Check the dimensions of unravelled maps for any inconsistencies. Perform a resampling
    if necessary to maintain consistent dimensions."""
    #sometimes unravelling maps leads to slight variations in the unraveeled image dimensions.
    #check dimensions of maps and resample to dimensions of first image in sequence if need be.
    #note that maps.shape lists the dimensions as (y,x) but maps.resample takes the arguments
    #as (x,y).
    ref_dim = new_maps[0].shape[::-1]
    resampled_maps = []
    for i, Map in enumerate(new_maps):
        if Map.shape[::-1] != ref_dim:
            tmp = Map.resample(ref_dim, method='linear', minusone=True)
            print('Notice: resampling performed on frame ' + str(i) +
                  ' to maintain consistent dimensions.')
            resampled_maps.append(tmp)
        else:
            resampled_maps.append(Map)
    return resampled_maps


def linesampleindex(a, b, np=1000):
    """ Get the indices in an array along a line"""
    x, y = np.linspace(a[0], b[0], np), np.linspace(a[1], b[1], np)
    xi = x.astype(np.int)
    yi = y.astype(np.int)
    return xi, yi


def map_diff(maps):
    """ calculate running difference images """
    diffs = []
    for i in range(0,len(maps)-1):
        # take the difference
        diffmap = (maps[i+1] - maps[i])
        # keep
        diffs.append(diffmap)
    return diffs

def map_threshold(maps, factor):
    threshold_maps = []
    for i in range(1,len(maps)):
        sqrt_map = np.sqrt(maps[i]) * factor
        threshold_maps.append(sqrt_map)
    return threshold_maps

def map_binary(diffs, threshold_maps):
    """turn difference maps into binary images"""
    binary_maps = []
    for i in range(0,len(diffs)):
        #for values > threshold_map in the diffmap, return True, otherwise False
        filtered_map = diffs[i] > threshold_maps[i]
        binary_maps.append(filtered_map)
    return binary_maps


'''Ideas

extract a submap that is where we expect the wave to be and just concentrate on
that region
 - will speed up processing

 adaptive thresholding; as time increases the threshold decreases, anticipating
 the decreasing amplitude of the wave.

'''

def hough_detect(diffs, vote_thresh=12):
    """ Use the Hough detection method to detect lines in the data.
    With enough lines, you can fill in the wave front."""
    detection=[]
    print("Performing hough transform on binary maps...")
    for img in diffs:
        # Perform the hough transform on each of the difference maps
        transform, theta, d = hough(img)

        # Filter the hough transform results and find the best lines in the
        # data.  Keep detections that exceed the Hough vote threshold.
        indices = (transform>vote_thresh).nonzero()
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
        invTransform = sunpy.make_map(np.zeros(img.shape), img._original_header)
        lines = probabilistic_hough(img, ph_kwargs)
        if lines is not None:
            for line in lines:
                pos1=line[0]
                pos2=line[1]
                fillLine(pos1,pos2,invTransform)
        detection.append(invTransform)
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

def check_fit(result):
    """Remove bad fit results that are not caught by the fitting flag. Returns
    a blank list if the fit is deemed to be bad, otherwise returns the input unchanged."""
    #check that the location of the wave lies within +90 and -90 degrees
    if result[0][1] > 90.0 or result[0][1] < -90.0:
        result=[]
        return result
    #check that the width of the wave is not too large (> 15)
    if result[0][2] > 15.0:
        result=[]
        return result
    return result

def fit_wavefront(diffs, detection):
    """Fit the wavefront that has been detected by the hough transform.
    Simplest case is to fit along the y-direction for some x or range of x."""
    dims=diffs[0].shape
    answers=[]
    wavefront_maps=[]
    for i in range (0, len(diffs)):
        if (detection[i].max() == 0.0):
            #if the 'detection' array is empty then skip this image
            fit_map=sunpy.make_map(np.zeros(dims),diffs[0]._original_header)
            print("Nothing detected in image " + str(i) + ". Skipping.")
            answers.append([])
            wavefront_maps.append(fit_map)
        else:
            #if the 'detection' array is not empty, then fit the wavefront in the image
            img = diffs[i]
            fit_map=np.zeros(dims)

            #get the independent variable for the columns in the image
            x=(np.linspace(0,dims[0],num=dims[0])*img.scale['y']) + img.yrange[0]
            
            #use 'detection' to guess the centroid of the Gaussian fit function
            guess_index=detection[i].argmax()
            guess_index=np.unravel_index(guess_index,detection[i].shape)
            guess_position=x[guess_index[0]]
            
            print("Analysing wavefront in image " + str(i))
            column_fits=[]
            #for each column in image, fit along the y-direction a function to find wave parameters
            for n in range (0,dims[1]):
                #guess the amplitude of the Gaussian fit from the difference image
                guess_amp=np.float(img[guess_index[0],n])
                
                #put the guess input parameters into a vector
                guess_params=[guess_amp,guess_position,5]

                #get the current image column
                y=img[:,n]
                y=y.flatten()                
                #call Albert's fitting function
                result = util.fitfunc(x,y,'Gaussian',guess_params)

                #define a Gaussian function. Messy - clean this up later
                gaussian = lambda p,x: p[0]/np.sqrt(2.*np.pi)/p[2]*np.exp(-((x-p[1])/p[2])**2/2.)
                
                #Draw the Gaussian fit for the current column and save it in fit_map
                #save the best-fit parameters in column_fits
                #only want to store the successful fits, discard the others.
                #result contains a pass/fail integer. Keep successes ( ==1).
                if result[1] == 1:
                    #if we got a pass integer, perform some other checks to eliminate unphysical values
                    result=check_fit(result)    
                    column_fits.append(result)
                    if result != []:
                        fit_column = gaussian(result[0],x)
                    else:
                        fit_column = np.zeros(len(x))
                else:
                    #if the fit failed then save as zeros/null values
                    result=[]
                    column_fits.append(result)
                    fit_column = np.zeros(len(x))
                
                #draw the Gaussian fit for the current column and save it in fit_map
                #gaussian = lambda p,x: p[0]/np.sqrt(2.*np.pi)/p[2]*np.exp(-((x-p[1])/p[2])**2/2.)
                    
                #save the drawn column in fit_map
                fit_map[:,n] = fit_column
            #save the fit parameters for the image in 'answers' and the drawn map in 'wavefront_maps'
            fit_map=sunpy.make_map(fit_map,diffs[0]._original_header)
            answers.append(column_fits)
            wavefront_maps.append(fit_map)

    #now get the mean values of the fitted wavefront, averaged over all x
    #average_fits=[]
    #for ans in answers:
    #   cleaned_answers=[]
    #  for k in range(0,len(ans)):
    #      #ans[:,1] contains a pass/fail integer. Keep successes (==1), discard the rest
    #      if ans[k][1] == 1:
    #          tmp=ans[k][0]
    #          cleaned_answers.append(tmp)
    #      else:
    #          cleaned_answers.append([])
    #  #get the mean of each fit parameter for this image and store it
    #  #average_fits.append(np.mean(g,axis=0))
        
    return answers, wavefront_maps

def wavefront_velocity(answers):
    """calculate wavefront velocity based on fit parameters for each column of an image or set of images"""
    velocity=[]
    for i in range(0,len(answers)):
        v=[]
        if i==0:
            velocity.append([])
        else:
            #skip blank entries of answers
            if answers[i] == [] or answers[i-1] == []:
                velocity.append([])
            else:
                for j in range(0,len(answers[i])):
                    #want to ignore null values for wave position
                    if answers[i][j] == [] or answers[i-1][j] == []:
                        vel=[]
                    else:             
                        vel=answers[i][j][0][1] - answers[i-1][j][0][1]
                    v.append(vel)
                velocity.append(v)
    return velocity

def wavefront_position_and_width(answers):
    """get wavefront position and width based on fit parameters for each column of an image or set of images"""
    position=[]
    width=[]
    for i in range(0,len(answers)):
        p=[]
        w=[]
        if answers[i] == []:
            position.append([])
            width.append([])
        else:
            for j in range(0,len(answers[i])):
                #want to ignore null values for wave position
                if answers[i][j] == []:
                    pos=[]
                    wid=[]
                else:
                    pos=answers[i][j][0][1]
                    wid=answers[i][j][0][2]
                p.append(pos)
                w.append(wid)
            position.append(p)
            width.append(w)
    return position,width

def fillLine(pos1,pos2,img):
    shape=img.shape
    ny = shape[0]
    nx = shape[1]
    if pos2[0] == pos1[0]:
        m = 9999
    else:
        m = (pos2[1] - pos1[1]) / (pos2[0] - pos1[0])
        
    constant = (pos2[1] - m*pos2[0])
    
    for x in range(pos1[0],pos2[0]):
        y = m*x + constant
        if y <= ny-1 and y>= 0:
            img[y,x] = 255

    return img

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


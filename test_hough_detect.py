from sim import wave2d
from visualize import visualize
from skimage.transform import hough
from skimage.transform import probabilistic_hough
from skimage.morphology import greyscale_dilate
import numpy as np
import pylab as plt
import sunpy

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

    return img,y
        

m2deg = 360./(2*3.1415926*6.96e8)

max_steps = 20

waveparams = {
    # "cadence": 12., #seconds
    "cadence": 30., #seconds
    
    "hglt_obs": 0., #degrees
    "rotation": 360./(27.*86400.), #degrees/s, rigid solar rotation
    
    #Wave parameters that are initial conditions
    "direction": 65., #degrees, measured CCW from HG +latitude
    "epi_lat": 30., #degrees, HG latitude of wave epicenter
    "epi_lon": 45., #degrees, HG longitude of wave epicenter
    
    #Wave parameters that can evolve over time
    #The first element is constant in time
    #The second element (if present) is linear in time
    #The third element (if present) is quadratic in time
    #Be very careful of non-physical behavior
    "width": [90., 0.0], #degrees, full angle in azimuth, centered at 'direction'
    "wave_thickness": [6.0e6*m2deg,2.0e4*m2deg], #degrees, sigma of Gaussian profile in longitudinal direction
    "wave_normalization": [1.], #integrated value of the 1D Gaussian profile
    #"speed": [9.33e5*m2deg, -1.495e3*m2deg], #degrees/s, make sure that wave propagates all the way to lat_min for polynomial speed
    "speed": [9.33e5*m2deg, 0],
    
    #Random noise parameters
    "noise_type": "Poisson", #can be None, "Normal", or "Poisson"
    "noise_scale": 0.05,
    "noise_mean": 1.,
    "noise_sdev": 1.,
    
    #Structured noise parameters
    "struct_type": None, #can be None, "Arcs", or "Random"
    "struct_scale": 5.,
    "struct_num": 10,
    "struct_seed": 13092,
    
    "max_steps": max_steps,
    
    "clean_nans": True,
    
    #HG grid, probably would only want to change the bin sizes
    "lat_min": -90.,
    "lat_max": 90.,
    "lat_bin": 0.2,
    "lon_min": -180.,
    "lon_max": 180.,
    "lon_bin": 5.,
    
    #HPC grid, probably would only want to change the bin sizes
    "hpcx_min": -1228.8,
    "hpcx_max": 1228.8,
    "hpcx_bin": 2.4,
    "hpcy_min": -1228.8,
    "hpcy_max": 1228.8,
    "hpcy_bin": 2.4
}

wave_maps = wave2d.simulate(params, verbose = True)
    
    #To get simulated HG' maps (centered at wave epicenter):
wave_maps_raw = wave2d.simulate_raw(params)
wave_maps_raw_noise = wave2d.add_noise(params, wave_maps_raw)
    
visualize(wave_maps)
    
import util
    
new_wave_maps = []
    
for wave in wave_maps:
    print("Unraveling map at "+str(wave.date))
    new_wave_maps += [util.map_hpc_to_hg_rotate(wave, epi_lon = params.get('epi_lon'), epi_lat = params.get('epi_lat'), xbin = 5, ybin = 0.2)]

input_maps = new_wave_maps

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

ndiff = len(input_maps)-1

# difference threshold
diffthresh = 0.2

# Hough transform voting threshold
votethresh = 10

# shape of the data
imgShape = input_maps[0].shape

# storage for the detection
detection = []
diffs = []

#temp = 255*(abs(input_maps[14] - input_maps[13]) > diffthresh)

for i in range(0,ndiff):
    # difference map - create separate maps for +ve and -ve differences
    diffmap_plus = 255*((input_maps[i+1] - input_maps[i]) > diffthresh)
    diffmap_minus = 255*((input_maps[i] - input_maps[i+1]) > diffthresh)
    # keep
    diffs.append(diffmap_plus)

    # extract the image
    img = diffmap_plus
    img2 = diffmap_minus

    # Perform the hough transform on the positive and negative difference maps separately
    transform,theta,d = hough(img)
    transform2,theta2,d2 = hough(img2)
    
    # Filter the hough transform results and find the best lines
    # in the data
    #indices =  (transform >votethresh).nonzero()
    #indices = (transform == transform.max()).nonzero()
    #instead of getting all lines above some threshold, just get the *strongest* line only
    #from the positive diffmap and the negative diffmap. May get more than 2 lines due to ties
    #in the accumulator
    indices=((transform == transform.max()) + (transform2 == transform2.max())).nonzero()
    distances = d[indices[0]]
    theta = theta[indices[1]]
    n =len(indices[1])
    print("Found " + str(n) + " lines.")

    # Perform the inverse transform to get a series of rectangular
    # images that show where the wavefront is.
    invTransform = sunpy.make_map(np.zeros(imgShape),input_maps[i+1]._original_header)
    # invTransform.data = np.zeros(imgShape)
    
    for i in range(0,n):
        nextLine = htLine( distances[i],theta[i], np.zeros(shape=imgShape) )
        invTransform = invTransform + nextLine

    # Dump the inverse transform back into a series of maps
    detection.append(invTransform)


visualize(diffs)
visualize(detection)

from matplotlib import cm
from matplotlib import colors

wmap = sunpy.make_map([wave_maps[max_steps/2], wave_maps[0]], type = "composite")
wmap.set_colors(1, cm.Reds)
wmap.set_alpha(1,0.1)
#wmap.set_norm(1, colors.Normalize(0.1,1))
wmap.show()

pmap = sunpy.make_map([detection[max_steps/2],input_maps[max_steps/2]], type ="composite")
pmap.set_alpha(1,0.6)
pmap.set_colors(0, cm.Blues)
pmap.set_colors(1, cm.Reds)
pmap.show()


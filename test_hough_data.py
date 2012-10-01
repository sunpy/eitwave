from sim import wave2d
from visualize import visualize
from scikits.image.transform import hough
from scikits.image.morphology import greyscale_dilate
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

    return img

m2deg = 360./(2*3.1415926*6.96e8)

cube = sunpy.make_map("/Users/schriste/Downloads/eitdata_19970512/*.fits", type = "cube")
dmap = cube[2] - cube[1]
dmap.show()

# need an even number of maps so get rid of one
cube = cube[0:4]

import util

tmap = util.map_hpc_to_hg(dmap)

ttmap = util.map_hpc_to_hg_rotate(dmap, epi_lon = 9.5, epi_lat = 20.44)
input_maps = []

for map in cube:
    print("Unraveling map at "+str(map.date))
    input_maps += [util.map_hpc_to_hg_rotate(map, epi_lon = 9.5, epi_lat = 20.44)]

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
diffthresh = 130

# Hough transform voting threshold
votethresh = 30

# shape of the data
imgShape = input_maps[0].shape

# storage for the detection
detection = []
diffs = []

diffmap = 255*(abs(input_maps[1] - input_maps[0]) > diffthresh)

for i in range(0,ndiff):
    # difference map
    diffmap = 255*(abs(input_maps[i+1] - input_maps[i]) > diffthresh)

    # keep
    diffs.append(diffmap)

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
    print("Found " + str(n) + " lines.")

    # Perform the inverse transform to get a series of rectangular
    # images that show where the wavefront is.
    invTransform = sunpy.map.BaseMap(input_maps[i+1])
    invTransform.data = np.zeros(imgShape)
    for i in range(0,n):
        nextLine = htLine( distances[i],theta[i], np.zeros(shape=imgShape) )
        invTransform = invTransform + nextLine

    # Dump the inverse transform back into a series of maps
    detection.append(invTransform)


visualize(diffs)
visualize(detection)

from matplotlib import cm
from matplotlib import colors

#wmap = sunpy.make_map(input_maps[max_steps/2], wave_maps[0], type = "composite")
#wmap.set_colors(1, cm.Reds)
#wmap.set_alpha(1,0.1)
#wmap.set_norm(1, colors.Normalize(0.1,1))
#wmap.show()

#pmap = sunpy.make_map(detection[max_steps/2],input_maps[max_steps/2], type ="composite")
#pmap.set_alpha(1,0.6)
#pmap.set_colors(0, cm.Blues)
#pmap.set_colors(1, cm.Reds)
#pmap.show()
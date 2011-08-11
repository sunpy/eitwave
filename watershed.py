#!/usr/bin/env python
#-*- coding:utf-8 -*-
import os
import sunpy
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from sunpy.map.sources.sdo import AIAMap
"""
Prototype code
"""

#
# Test data:
#    http://www.sunpy.org/research/agu-2011/testdata1.zip
#
basedir = "" # Location of extracted test data

files = ['AIA20101016_191255_0193.fits',
         'AIA20101016_191246_0193.fits',
         'AIA20101016_191231_0193.fits',
         'AIA20101016_191222_0193.fits',
         'AIA20101016_191207_0193.fits']

filepaths = [os.path.join(basedir, i) for i in files]

# Creates maps
maps = []

for path in filepaths:
    maps.append(sunpy.Map(path)[512:2048, 2048:3584])

# Blur maps
blurred = []
for map_ in maps:
    blurred.append(AIAMap(ndimage.gaussian_filter(map_, 10), map_.header))

# Plot map 3 - map 1
diff = blurred[2] - blurred[0]

# Distance transform
#dist = ndimage.distance_transform_edt(diff > 20000)
#dist = dist.max() - dist
#dist -= dist.min()
#dist = dist / float(dist.ptp()) * 65535
#dist = dist.astype(np.uint16).clip(0, 40000)

# Labeling
labels, nr_nuclei = ndimage.label(diff.clip(diff.min(), 0))

print("Number of nuclei found: %d" % nr_nuclei)

# Watershed
areas = np.array([(labels == s).sum() for s in np.arange(nr_nuclei) + 1])
max_area_label = areas.argmax() + 1
#ax = pylab.imshow(labels == max_area_label)

# Plot
fig = plt.figure()
axes = plt.imshow(labels == max_area_label, origin='lower')
plt.colorbar(axes)
plt.show()

########################################

#ax = pylab.imshow(labels)
#pylab.colorbar(ax)
#pylab.show()

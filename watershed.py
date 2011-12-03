#!/usr/bin/env python
#-*- coding:utf-8 -*-
import sunpy
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from sunpy.map.sources.sdo import AIAMap
import util
"""
Prototype code
"""

cube = sunpy.make_map("/Users/schriste/Downloads/data2/AIA2010*", type = "cube")[:, 512:2048, 2048:3584]

# Blur maps
blurred = []
for map_ in cube:
    blurred.append(AIAMap(ndimage.gaussian_filter(map_, 10), map_.header))

# Plot map 3 - map 1
diff = blurred[2] - blurred[4]

tmap = util.map_hpc_to_hg(diff)

# Labeling
labels, nr_nuclei = ndimage.label(diff.clip(diff.min(), 0))

print("Number of nuclei found: %d" % nr_nuclei)

# Watershed
areas = np.array([(labels == s).sum() for s in np.arange(nr_nuclei) + 1])
max_area_label = areas.argmax() + 1

# Plot
fig = plt.figure()
axes = plt.imshow(labels == max_area_label, origin='lower')
plt.colorbar(axes)
plt.show()


# cube = sunpy.MapCube("~/Dropbox/python/eitwave/data/data/ssw_*")
cube = sunpy.make_map("/Users/schriste/Downloads/data/ssw_cutout_20101016_*_.fts", type = "cube")

scube = [map.submap([300,800], [-550,-200]) for map in cube]
scube = [map.submap([1948,2781], [531,1114], units = 'pixels') for map in cube]

over_expmap = scube[::2]
under_expmap = scube[1::2]

i = 0
(scube[i+2] - scube[i]).show()

# Blur maps
blurred = []
for map_ in over_expmap:
    blurred.append(AIAMap(ndimage.gaussian_filter(map_, 10), map_.header))

# Plot map 3 - map 1
diff = blurred[2] - blurred[4]

# Labeling
labels, nr_nuclei = ndimage.label(diff.clip(diff.min(), 0))

print("Number of nuclei found: %d" % nr_nuclei)

# Watershed
areas = np.array([(labels == s).sum() for s in np.arange(nr_nuclei) + 1])
max_area_label = areas.argmax() + 1

# Plot
fig = plt.figure()
axes = plt.imshow(labels == max_area_label, origin='lower')
plt.colorbar(axes)
plt.show()

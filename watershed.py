#!/usr/bin/env python
#-*- coding:utf-8 -*-
import sunpy
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from sunpy.map.sources.sdo import AIAMap
"""
Prototype code
"""
cube = sunpy.MapCube("data/AIA2010*")[:, 512:2048, 2048:3584]

# Blur maps
blurred = []
for map_ in cube:
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

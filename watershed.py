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

# the following code will convert an AIA image to HG coordinates.
import numpy as np
import sunpy
aia = sunpy.Map("data/AIA20101016_191207_0193.fits")
saia = aia.resample([500,500])

lon_bin = 1
lat_bin = 1
lon = np.arange(-90,90, lon_bin)
lat = np.arange(-90, 90, lat_bin)
lon_grid, lat_grid = np.meshgrid(lon, lat)

from sunpy.wcs  import wcs as wcs
from scipy.interpolate import griddata

newgrid = wcs.convert_hg_hpc(saia.header, lon_grid, lat_grid, units = 'arcsec')

x,y = wcs.convert_pixel_to_data(saia.header)

points = np.vstack((x.ravel(), y.ravel())).T
values = np.array(saia).ravel()

grid = griddata(points, values,newgrid, method="linear")

import matplotlib.pyplot as plt
plt.subplot(111)
plt.imshow(grid, extent = [-90,90,-90,90])
plt.show()

dict_header = {
        "cdelt1": lon_bin,
        "naxis1": len(lon),
        "crval1": lon.min(),
        "crpix1": 0,
        "cunit1": "deg",
        "ctype1": "HG",
        "cdelt2": lat_bin,
        "naxis2": len(lat),
        "crval2": lat.min(),
        "crpix2": 0,
        "cunit2": "deg",
        "ctype2": "HG",
        "hglt_obs": 0,
        "hgln_obs": 0,
        "rsun_obs": aia.header.get('rsun_obs'),
        "rsun_ref": aia.header.get('rsun_ref'),
        "dsun_obs": aia.header.get('dsun_obs'),
        'wavelnth': aia.header.get('wavelnth')
}
    
header = sunpy.map.MapHeader(dict_header)
transformed_map = sunpy.map.BaseMap(grid, header)


import sunpy
from sunpy.wcs import wcs
import util

aia = sunpy.Map(sunpy.AIA_171_IMAGE).resample([500,500])

tmap = util.map_hpc_to_hg(aia)
tmap.show()

tmap_rot = util.map_hpc_to_hg_rotate(aia, epi_lon = 0, epi_lat = 0)
tmap_rot.show()
# BUG at the edges where the Sun wraps around, HPC to HG gives the same longitude
# for points behind and in front of the limb.

ttmap = util.map_hg_to_hpc(tmap)
ttmap.show()

# WARNING!
# if the following line is run after the previousl line than it gives an error
# but not otherwise...not sure why...
ttmap = util.map_hg_to_hpc(tmap, xbin = 5, ybin = 5)
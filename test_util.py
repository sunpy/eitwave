import sunpy
import util
from sunpy.wcs import wcs

aia = sunpy.Map(sunpy.AIA_171_IMAGE).resample([500,500])

tmap = util.map_hpc_to_hg(aia)

ttmap = util.map_hg_to_hpc(tmap)

# WARNING!
# if the following line is run after the previously line than it gives an error
# but not otherwise...not sure why...
ttmap = util.map_hg_to_hpc(tmap, xbin = 5, ybin = 5)
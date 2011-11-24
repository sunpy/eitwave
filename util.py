import sunpy
from sunpy.wcs  import wcs as wcs
from scipy.interpolate import griddata
import numpy as np
    
def map_to_hg(map, xbin = 1, ybin = 1):
    """Take a map (like an AIA map) and convert it from HPC to HG."""

    x,y = wcs.convert_pixel_to_data(map.header)
    lon_map, lat_map = wcs.convert_hpc_hg(map.header, x, y)

    lon_bin = xbin
    lat_bin = ybin 
    lon_range = (np.nanmin(lon_map), np.nanmax(lon_map))
    lat_range = (np.nanmin(lat_map), np.nanmax(lat_map))

    lon = np.arange(lon_range[0], lon_range[1], lon_bin)
 
    lat = np.arange(lat_range[0], lat_range[1], lat_bin)
    lon_grid, lat_grid = np.meshgrid(lon, lat)

    newgrid = wcs.convert_hg_hpc(map.header, lon_grid, lat_grid, units = 'arcsec')

    points = np.vstack((x.ravel(), y.ravel())).T
    values = np.array(map).ravel()
    newdata = griddata(points, values,newgrid, method="linear")

    header = map.header
    header['cdelt1'] = lon_bin
    header['naxis1'] = len(lon)
    header['crval1'] = lon.min()
    header['crpix1'] = 0
    header['crpix2'] = 0
    header['cunit1'] = "deg"
    header['ctype1'] = "HG"
    header['cdelt2'] = lat_bin
    header['naxis2'] = len(lat)
    header['crval2'] = lat.min()
    header['cunit2'] = "deg"
    header['ctype2'] = "HG"

    transformed_map = sunpy.map.BaseMap(newdata, header)

    transformed_map.cmap = map.cmap
    transformed_map.name = map.name
    transformed_map.date = map.date
    transformed_map.center = {
        "x": wcs.get_center(header, axis='x'),
        "y": wcs.get_center(header, axis='y')}

    return transformed_map

def map_to_hpc(map, xbin = 10, ybin = 10):
    """Take a map and convert it from HG to HPC."""
    
    lon,lat = wcs.convert_pixel_to_data(map.header)
    x_map, y_map = wcs.convert_hg_hpc(map.header, lon, lat)

    x_range = (np.nanmin(x_map)*60*60, np.nanmax(x_map)*60*60)
    y_range = (np.nanmin(y_map)*60*60, np.nanmax(y_map)*60*60)

    x = np.arange(x_range[0], x_range[1], xbin)
    y = np.arange(y_range[0], y_range[1], ybin)
    xgrid, ygrid = np.meshgrid(x, x)

    newgrid = wcs.convert_hpc_hg(map.header, xgrid, ygrid)

    points = np.vstack((lon.ravel(), lat.ravel())).T
    values = np.array(map).ravel()
    newdata = griddata(points, values, newgrid, method="linear")

    # now grab the original map header and update it
    header = map.header
    header["cdelt1"] = xbin
    header["naxis1"] = len(x)
    header["crval1"] = x.min()
    header["crpix1"] = 0
    header["cunit1"] = "arcsec"
    header["ctype1"] = "HPLT-TAN"
    header["cdelt2"] = ybin
    header["naxis2"] = len(y)
    header["crval2"] = lat.min()
    header["crpix2"] = 0
    header["cunit2"] = "arcsec"
    header["ctype2"] = "HPLT-TAN"
    
    transformed_map = sunpy.map.BaseMap(newdata, header)

    transformed_map.cmap = map.cmap
    transformed_map.name = map.name
    transformed_map.date = map.date
    
    transformed_map.center = {
        "x": wcs.get_center(header, axis='x'),
        "y": wcs.get_center(header, axis='y')}

    return transformed_map

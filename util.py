import sunpy
from sunpy.wcs  import wcs as wcs
from scipy.interpolate import griddata
import numpy as np
from sim.wave2d.wave2d import euler_zyz

__authors__ = ["Steven Christe"]
__email__ = "steven.d.christe@nasa.gov"

def map_hpc_to_hg(map, xbin = 1, ybin = 1):
    """Take a map (like an AIA map) and convert it from HPC to HG."""

    x,y = wcs.convert_pixel_to_data(map.header)
    lon_map, lat_map = wcs.convert_hpc_hg(map.header, x, y)

    lon_bin = xbin
    lat_bin = ybin 
    lon_range = (np.nanmin(lon_map), np.nanmax(lon_map))
    lat_range = (np.nanmin(lat_map), np.nanmax(lat_map))

    lon = np.arange(lon_range[0], lon_range[1], lon_bin)
    lat = np.arange(lat_range[0], lat_range[1], lat_bin)
    newgrid = np.meshgrid(lon, lat)

    # newgrid = wcs.convert_hg_hpc(map.header, lon_grid, lat_grid, units = 'arcsec')
    points = np.vstack((lon_map.ravel(), lat_map.ravel())).T
    values = np.array(map).ravel()

    # get rid of all of the bad (nan) indices (i.e. those off of the sun)
    index = np.isfinite(points[:,0]) * np.isfinite(points[:,1])
    points = np.vstack((points[index,0], points[index,1])).T
    values = values[index]

    newdata = griddata(points, values, newgrid, method="cubic")

    header = map.header
    header['CDELT1'] = lon_bin
    header['NAXIS1'] = len(lon)
    header['CRVAL1'] = lon.min()
    header['CRPIX1'] = 1
    header['CRPIX2'] = 1
    header['CUNIT1'] = "deg"
    header['CTYPE1'] = "HG"
    header['CDELT2'] = lat_bin
    header['NAXIS2'] = len(lat)
    header['CRVAL2'] = lat.min()
    header['CUNIT2'] = "deg"
    header['CTYPE2'] = "HG"

    transformed_map = sunpy.map.BaseMap(newdata, header)

    transformed_map.cmap = map.cmap
    transformed_map.name = map.name
    transformed_map.date = map.date
    transformed_map.center = {
        "x": wcs.get_center(header, axis='x'),
        "y": wcs.get_center(header, axis='y')}

    return transformed_map

def map_hg_to_hpc(map, xbin = 10, ybin = 10):
    """Take a map in heliographic coordinates (HG) and convert it to 
    helioprojective cartesian coordinates (HPC)."""

    lon,lat = wcs.convert_pixel_to_data(map.header)
    x_map, y_map = wcs.convert_hg_hpc(map.header, lon, lat, units ='arcsec')
    
    x_range = (np.nanmin(x_map), np.nanmax(x_map))
    y_range = (np.nanmin(y_map), np.nanmax(y_map))
    
    x = np.arange(x_range[0], x_range[1], xbin)
    y = np.arange(y_range[0], y_range[1], ybin)
    newgrid = np.meshgrid(x, y)
    
    # newgrid = wcs.convert_hpc_hg(map.header, xgrid/(3600), ygrid/(3600))
    
    points = np.vstack((x_map.ravel(), y_map.ravel())).T
    values = np.array(map).ravel()
    newdata = griddata(points, values, newgrid, method="linear")
    
    # now grab the original map header and update it
    header = map.header
    header["CDELT1"] = xbin
    header["NAXIS1"] = len(x)
    header["CRVAL1"] = x.min()
    header["CRPIX1"] = 1
    header["CUNIT1"] = "arcsec"
    header["CTYPE1"] = "HPLT-TAN"
    header["CDELT2"] = ybin
    header["NAXIS2"] = len(y)
    header["CRVAL2"] = y.min()
    header["CRPIX2"] = 1
    header["CUNIT2"] = "arcsec"
    header["CTYPE2"] = "HPLT-TAN"
    
    transformed_map = sunpy.map.BaseMap(newdata, header)
    
    transformed_map.cmap = map.cmap
    transformed_map.name = map.name
    transformed_map.date = map.date

    transformed_map.center = {
        "x": wcs.get_center(header, axis='x'),
        "y": wcs.get_center(header, axis='y')}


    return transformed_map

def map_hpc_to_hg_rotate(map, epi_lon = 0, epi_lat = 0, xbin = 1, ybin = 1):
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

    x, y, z = sunpy.wcs.convert_hg_hcc_xyz(map.header,
                                           lon_grid, lat_grid)

    zpp, xpp, ypp = euler_zyz((z, x, y), (epi_lon, 90.-epi_lat, 0.))
    
    newgrid = sunpy.wcs.convert_hcc_hpc(current_wave_map.header, xpp, ypp)

    points = np.vstack((x.ravel(), y.ravel())).T
    values = np.array(map).ravel()
    newdata = griddata(points, values, newgrid, method="linear")

    header = map.header
    header['CDELT1'] = lon_bin
    header['NAXIS1'] = len(lon)
    header['CRVAL1'] = lon.min()
    header['CRPIX1'] = 0
    header['CRPIX2'] = 0
    header['CUNIT1'] = "deg"
    header['CTYPE1'] = "HG"
    header['CDELT2'] = lat_bin
    header['NAXIS2'] = len(lat)
    header['CRVAL2'] = lat.min()
    header['CUNIT2'] = "deg"
    header['CTYPE2'] = "HG"

    transformed_map = sunpy.map.BaseMap(newdata, header)

    transformed_map.cmap = map.cmap
    transformed_map.name = map.name
    transformed_map.date = map.date
    transformed_map.center = {
        "x": wcs.get_center(header, axis='x'),
        "y": wcs.get_center(header, axis='y')}

    return transformed_map

import sunpy
from sunpy import wcs
from scipy.interpolate import griddata
import numpy as np
#from sim.wave2d.wave2d import euler_zyz
#from matplotlib import colors

__authors__ = ["Steven Christe"]
__email__ = "steven.d.christe@nasa.gov"

def map_hpc_to_hg(smap, xbin = 1, ybin = 1):
    """Take a map (like an AIA map) and convert it from HPC to HG."""

    #x,y = wcs.convert_pixel_to_data(map.header)
    x,y = wcs.convert_pixel_to_data(smap.shape[1],
                                    smap.shape[0],
                                    smap.scale['x'], 
                                    smap.scale['y'],
                                    smap.center['x'],
                                    smap.center['y'],   
                                    smap.reference_coordinate['x'],
                                    smap.reference_coordinate['y'],
                                    smap.coordinate_system['x'])
    
    #lon_map, lat_map = wcs.convert_hpc_hg(map.header, x, y)
    lon_map, lat_map = wcs.convert_hpc_hg(smap.rsun_meters,
                                          smap.dsun,
                                          smap.scale['x'],
                                          smap.scale['y'],
                                          smap.heliographic_latitude,
                                          smap.carrington_longitude,
                                          x, y)
    
    lon_bin = xbin
    lat_bin = ybin 
    lon_range = (np.nanmin(lon_map), np.nanmax(lon_map))
    lat_range = (np.nanmin(lat_map), np.nanmax(lat_map))

    lon = np.arange(lon_range[0], lon_range[1], lon_bin)
    lat = np.arange(lat_range[0], lat_range[1], lat_bin)
    newgrid = np.meshgrid(lon, lat)

    # newgrid = wcs.convert_hg_hpc(map.header, lon_grid, lat_grid, units = 'arcsec')
    points = np.vstack((lon_map.ravel(), lat_map.ravel())).T
    values = np.array(smap).ravel()

    # get rid of all of the bad (nan) indices (i.e. those off of the sun)
    index = np.isfinite(points[:,0]) * np.isfinite(points[:,1])
    points = np.vstack((points[index,0], points[index,1])).T
  
    values = values[index]
    
    newdata = griddata(points, values, newgrid, method="linear")

    header = smap.header.copy()
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

def map_hg_to_hpc(smap, xbin = 10, ybin = 10):
    """Take a map in heliographic coordinates (HG) and convert it to 
    helioprojective cartesian coordinates (HPC)."""

    lon,lat = wcs.convert_pixel_to_data(smap.header)
    x_map, y_map = wcs.convert_hg_hpc(smap.header, lon, lat, units ='arcsec')
    
    x_range = (np.nanmin(x_map), np.nanmax(x_map))
    y_range = (np.nanmin(y_map), np.nanmax(y_map))
    
    x = np.arange(x_range[0], x_range[1], xbin)
    y = np.arange(y_range[0], y_range[1], ybin)
    newgrid = np.meshgrid(x, y)
    
    # newgrid = wcs.convert_hpc_hg(map.header, xgrid/(3600), ygrid/(3600))
    
    points = np.vstack((x_map.ravel(), y_map.ravel())).T
    values = np.array(smap).ravel()
    newdata = griddata(points, values, newgrid, method="linear")
    
    # now grab the original map header and update it
    header = smap.header.copy()
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
    
    transformed_map.cmap = smap.cmap
    transformed_map.name = smap.name
    transformed_map.date = smap.date

    transformed_map.center = {
        "x": wcs.get_center(header, axis='x'),
        "y": wcs.get_center(header, axis='y')}

    return transformed_map

def map_hpc_to_hg_rotate(smap, epi_lon = 0, epi_lat = 0, xbin = 1, ybin = 1):
    """Take a map (like an AIA map) and convert it from HPC to HG."""

    #import sunpy
    #import util
    #from sunpy import wcs
    #import numpy as np
    #from scipy.interpolate import griddata
    from sim.wave2d.wave2d import euler_zyz
    #from matplotlib import colors
    
    # epi_lon = -10
    # epi_lat = 0
    
    #aia = sunpy.Map(sunpy.AIA_171_IMAGE).resample([500,500])
    # tmap = util.map_hpc_to_hg(aia)
    # tmap.show()
    
    #map = aia
    
    #x, y = wcs.convert_pixel_to_data(map.header)
    x, y = wcs.convert_pixel_to_data(smap.shape[1],
                                     smap.shape[0],
                                     smap.scale['x'], 
                                     smap.scale['y'],
                                     smap.reference_pixel['x'],
                                     smap.reference_pixel['y'],   
                                     smap.reference_coordinate['x'],
                                     smap.reference_coordinate['y'],
                                     smap.coordinate_system['x'])
    
    #hccx, hccy, hccz = wcs.convert_hpc_hcc_xyz(map.header, x, y)
    hccx, hccy, hccz = wcs.convert_hpc_hcc_xyz(smap.rsun_meters,
                                               smap.dsun,
                                               smap.units['x'],
                                               smap.units['y'],
                                               x,
                                               y)
    
    # rot_hccz, rot_hccy, rot_hccx = euler_zyz((hccz, hccx, hccy), (epi_lon, 90.-epi_lat, 0.))
    rot_hccz, rot_hccx, rot_hccy = euler_zyz((hccz, hccx, hccy), (0., epi_lat-90., -epi_lon))
    # zpp, xpp, ypp = euler_zyz(zxy_p, (0., hglt_obs, total_seconds*rotation))

    #lon_map, lat_map = wcs.convert_hcc_hg(map.header, rot_hccx, rot_hccy, z = rot_hccz)
    lon_map, lat_map = wcs.convert_hcc_hg(smap.rsun_meters,
                                          smap.heliographic_latitude,
                                          smap.carrington_longitude,
                                          rot_hccx, rot_hccy, z = rot_hccz)
    
    lon_bin = xbin
    lat_bin = ybin 
    lon_range = (np.nanmin(lon_map), np.nanmax(lon_map))
    lat_range = (np.nanmin(lat_map), np.nanmax(lat_map))
    
    lon = np.arange(lon_range[0], lon_range[1], lon_bin)
    lat = np.arange(lat_range[0], lat_range[1], lat_bin)
    newgrid = np.meshgrid(lon, lat)
    
    #This extra conversion and rotation back are needed to determine where to
    #mask out points that can't have corresponding data
    #ng_xyz = wcs.convert_hg_hcc_xyz(map.header, newgrid[0], newgrid[1])
    ng_xyz = wcs.convert_hg_hcc_xyz(smap.rsun_meters,
                                    smap.heliographic_latitude,
                                    smap.carrington_longitude,
                                    newgrid[0], newgrid[1])
    
    ng_zp, ng_xp, ng_yp = euler_zyz((ng_xyz[2], ng_xyz[0], ng_xyz[1]),
                                    (epi_lon, 90.-epi_lat, 0.))
    
    
    points = np.vstack((lon_map.ravel(), lat_map.ravel())).T
    values = np.array(smap).ravel()
        
    # get rid of all of the bad (nan) indices (i.e. those off of the sun)
    index = np.isfinite(points[:,0]) * np.isfinite(points[:,1])
    #points = np.vstack((points[index,0], points[index,1])).T
    points = points[index]
    values = values[index]
    
    newdata = griddata(points, values, newgrid, method="cubic")
    newdata[ng_zp < 0] = np.nan

    header = smap._original_header.copy()
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
    
    transformed_map.cmap = smap.cmap
    transformed_map.name = smap.name
    transformed_map.date = smap.date
    transformed_map.center['x'] = wcs.get_center(smap.shape[1], smap.scale['x'], smap.reference_coordinate['x'],smap.reference_pixel['x'])
    transformed_map.center['y'] = wcs.get_center(smap.shape[0], smap.scale['y'], smap.reference_coordinate['y'],smap.reference_pixel['y'])
    
    #transformed_map.show()
    
    return transformed_map
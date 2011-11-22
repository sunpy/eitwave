
def map_to_hg(map):
    """Take a map (like an AIA map) and convert it from HPC to HG."""
    import sunpy
    from sunpy.wcs  import wcs as wcs
    from scipy.interpolate import griddata
    import numpy as np

    # should determine the following quantities automatically
    lon_bin = 1
    lat_bin = 1
    lon_range = (-90.0,90.0)
    lat_range = (-90.0,90.0)
    lon = np.arange(lon_range[0], lon_range[1], lon_bin)
    lat = np.arange(lat_range[0], lat_range[1], lat_bin)
    lon_grid, lat_grid = np.meshgrid(lon, lat)

    newgrid = wcs.convert_hg_hpc(map.header, lon_grid, lat_grid, units = 'arcsec')

    x,y = wcs.convert_pixel_to_data(map.header)
    points = np.vstack((x.ravel(), y.ravel())).T
    values = np.array(map).ravel()
    newdata = griddata(points, values,newgrid, method="linear")

    # should grab the original maps header and replace only the following
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
        "rsun_obs": map.header.get('rsun_obs'),
        "rsun_ref": map.header.get('rsun_ref'),
        "dsun_obs": map.header.get('dsun_obs'),
        'wavelnth': map.header.get('wavelnth')
    }
    
    header = sunpy.map.MapHeader(dict_header)
    transformed_map = sunpy.map.BaseMap(newdata, header)

    transformed_map.cmap = map.cmap
    transformed_map.name = map.name
    transformed_map.date = map.date

    return transformed_map

def map_to_hpc(map):
    """Take a map and convert it from HG to HPC."""
    import sunpy
    from sunpy.wcs  import wcs as wcs
    from scipy.interpolate import griddata
    import numpy as np

    # should determine the following quantities automatically
    x_bin = 10
    y_bin = 10
    x_range = (-90.0,90.0)
    y_range = (-90.0,90.0)
    x = np.arange(lon_range[0], lon_range[1], lon_bin)
    y = np.arange(lat_range[0], lat_range[1], lat_bin)
    x_grid, x_grid = np.meshgrid(x, x)

    newgrid = wcs.convert_hpc_hg(map.header, x_grid, y_grid)

    lon,lat = wcs.convert_pixel_to_data(map.header)
    points = np.vstack((x.ravel(), y.ravel())).T
    values = np.array(map).ravel()
    newdata = griddata(points, values,newgrid, method="linear")

    # should grab the original maps header and replace only the following
    dict_header = {
        "cdelt1": x_bin,
        "naxis1": len(x),
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
        "rsun_obs": map.header.get('rsun_obs'),
        "rsun_ref": map.header.get('rsun_ref'),
        "dsun_obs": map.header.get('dsun_obs'),
        'wavelnth': map.header.get('wavelnth')
    }
    
    header = sunpy.map.MapHeader(dict_header)
    transformed_map = sunpy.map.BaseMap(newdata, header)

    transformed_map.cmap = map.cmap
    transformed_map.name = map.name
    transformed_map.date = map.date

    return transformed_map

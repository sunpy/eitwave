import sunpy
from sunpy import wcs
from scipy.interpolate import griddata
from scipy import optimize
import numpy as np
#from sim.wave2d.wave2d import euler_zyz
#from matplotlib import colors

__authors__ = ["Steven Christe"]
__email__ = "steven.d.christe@nasa.gov"

def map_hpc_to_hg(map, lon_bin = 1, lat_bin = 1):
    """
    Transform raw data in HPC coordinates to HG coordinates
    """
    return map_hpc_to_hg_rotate(map, lon_bin = lon_bin, lat_bin = lat_bin)

def map_hg_to_hpc(map, xbin = 1, ybin = 1):
    """
    Transform raw data in HG coordinates to HPC coordinates
    """
    return map_hg_to_hpc_rotate(map, xbin = xbin, ybin = ybin)

def map_hpc_to_hg_rotate(map, epi_lon = 0, epi_lat = 90, lon_bin = 1, lat_bin = 1):
    """
    Transform raw data in HPC coordinates to HG' coordinates

    HG' = HG, except center at wave epicenter
    """
    x, y = sunpy.wcs.convert_pixel_to_data([map.shape[1], map.shape[0]],
                                           [map.scale['x'], map.scale['y']],
                                           [map.reference_pixel['x'], map.reference_pixel['y']],
                                           [map.reference_coordinate['x'],map.reference_coordinate['y']])

    hccx, hccy, hccz = wcs.convert_hpc_hcc(x, y, angle_units=map.units['x'], z=True)

    rot_hccz, rot_hccx, rot_hccy = euler_zyz((hccz, hccx, hccy), (0., epi_lat-90., -epi_lon))

    lon_map, lat_map = wcs.convert_hcc_hg(rot_hccx, rot_hccy, b0_deg=map.heliographic_latitude,
                                          l0_deg=map.heliographic_longitude,z = rot_hccz)

    lon_range = (np.nanmin(lon_map), np.nanmax(lon_map))
    lat_range = (np.nanmin(lat_map), np.nanmax(lat_map))

    lon = np.arange(lon_range[0], lon_range[1], lon_bin)
    lat = np.arange(lat_range[0], lat_range[1], lat_bin)
    x_grid, y_grid = np.meshgrid(lon, lat)

    ng_xyz = wcs.convert_hg_hcc(x_grid, y_grid,b0_deg=map.heliographic_latitude,
                                l0_deg=map.heliographic_longitude,z=True)

    ng_zp, ng_xp, ng_yp = euler_zyz((ng_xyz[2], ng_xyz[0], ng_xyz[1]),
                                        (epi_lon, 90.-epi_lat, 0.))

        #ravel flattens the data into a 1D array
    points = np.vstack((lon_map.ravel(), lat_map.ravel())).T
    values = np.array(map).ravel()

    # get rid of all of the bad (nan) indices (i.e. those off of the sun)
    index = np.isfinite(points[:,0]) * np.isfinite(points[:,1])
    #points = np.vstack((points[index,0], points[index,1])).T
    points = points[index]
    values = values[index]

    newdata = griddata(points, values, (x_grid,y_grid), method="linear")
    newdata[ng_zp < 0] = np.nan

    dict_header = {
        'CDELT1': lon_bin,
        'NAXIS1': len(lon),
        'CRVAL1': lon.min(),
        'CRPIX1': 1,
        'CRPIX2': 1,
        'CUNIT1': "deg",
        'CTYPE1': "HG",
        'CDELT2': lat_bin,
        'NAXIS2': len(lat),
        'CRVAL2': lat.min(),
        'CUNIT2': "deg",
        'CTYPE2': "HG",
        'DATE_OBS': map.meta['date-obs']
    }

    header = dict_header
    transformed_map = sunpy.map.Map(newdata, header)
    transformed_map.name = map.name

    return transformed_map

def map_hg_to_hpc_rotate(map, epi_lon = 90, epi_lat = 0, xbin = 2.4, ybin = 2.4):
    """
    Transform raw data in HG' coordinates to HPC coordinates

    HG' = HG, except center at wave epicenter
    """

    #Origin grid, HG'
    lon_grid, lat_grid = sunpy.wcs.convert_pixel_to_data([map.shape[1], map.shape[0]],
                                                         [map.scale['x'], map.scale['y']],
                                                         [map.reference_pixel['x'], map.reference_pixel['y']],
                                                         [map.reference_coordinate['x'], map.reference_coordinate['y']])

    #Origin grid, HG' to HCC'
    #HCC' = HCC, except centered at wave epicenter
    x, y, z = sunpy.wcs.convert_hg_hcc(lon_grid, lat_grid, z=True)

    #Origin grid, HCC' to HCC''
    #Moves the wave epicenter to initial conditions
    #HCC'' = HCC, except assuming that HGLT_OBS = 0
    zpp, xpp, ypp = euler_zyz((z, x, y), (epi_lon, 90.-epi_lat, 0.))

    #Origin grid, HCC to HPC (arcsec)
    #xx, yy = sunpy.wcs.convert_hcc_hpc(current_wave_map.header, xpp, ypp)
    xx, yy = sunpy.wcs.convert_hcc_hpc(xpp, ypp)
                                       #dsun_meters=map.dsun)

    #Destination HPC grid
    hpcx_range = (np.nanmin(xx), np.nanmax(xx))
    hpcy_range = (np.nanmin(yy), np.nanmax(yy))

    hpcx = np.arange(hpcx_range[0], hpcx_range[1], xbin)
    hpcy = np.arange(hpcy_range[0], hpcy_range[1], ybin)
    newgrid_x, newgrid_y = np.meshgrid(hpcx, hpcy)

    #Coordinate positions (HPC) with corresponding map data
    points = np.vstack((xx.ravel(), yy.ravel())).T
    values = np.array(map).ravel()

    #2D interpolation from origin grid to destination grid
    newdata = griddata(points[zpp.ravel() >= 0], values[zpp.ravel() >= 0],
                       (newgrid_x, newgrid_y), method="linear")

    dict_header = {
        "CDELT1": xbin,
        "NAXIS1": len(hpcx),
        "CRVAL1": hpcx.min(),
        "CRPIX1": 1, #this makes hpcx.min() the center of the first bin
        "CUNIT1": "arcsec",
        "CTYPE1": "HPLN-TAN",
        "CDELT2": ybin,
        "NAXIS2": len(hpcy),
        "CRVAL2": hpcy.min(),
        "CRPIX2": 1, #this makes hpcy.min() the center of the first bin
        "CUNIT2": "arcsec",
        "CTYPE2": "HPLT-TAN",
        "HGLT_OBS": 0,
        "HGLN_OBS": 0,
        'DATE_OBS': map.meta['date-obs'],
    }

    header = sunpy.map.MapMeta(dict_header)

    transformed_map = sunpy.map.Map(newdata, header)
    transformed_map.name = map.name

    return transformed_map


def euler_zyz(xyz, angles):
    """
    Rotation with Euler angles defined in the ZYZ convention with left-handed
    positive sign convention.

    Parameters
    ----------
        xyz : tuple of ndarrays
            Input coordinates
        angles : tuple of scalars
            Angular rotations are applied in the following order
            * angles[2] is the angle CCW around Z axis (intrinsic rotation)
            * angles[1] is the angle CCW around Y axis (polar angle)
            * angles[0] is the angle CCW around Z axis (azimuth angle)

    Returns
    -------
        X, Y, Z : ndarray
            Output coordinates

    Notes
    -----
    angles = (phi, theta, psi) inverts angles = (-psi, -theta, -phi)

    References
    ----------
    https://en.wikipedia.org/wiki/Euler_angles#Matrix_orientation
        ("Left-handed positive sign convention", "ZYZ")

    Examples
    --------
    >>> wave2d.euler_zyz((np.array([1,0,0]),
                          np.array([0,1,0]),
                          np.array([0,0,1])),
                         (45,45,0))
    (array([ 0.5       , -0.70710678,  0.5       ]),
     array([ 0.5       ,  0.70710678,  0.5       ]),
     array([-0.70710678,  0.        ,  0.70710678]))
    """
    c1 = np.cos(np.deg2rad(angles[0]))
    s1 = np.sin(np.deg2rad(angles[0]))
    c2 = np.cos(np.deg2rad(angles[1]))
    s2 = np.sin(np.deg2rad(angles[1]))
    c3 = np.cos(np.deg2rad(angles[2]))
    s3 = np.sin(np.deg2rad(angles[2]))
    x = (c1*c2*c3-s1*s3)*xyz[0]+(-c3*s1-c1*c2*s3)*xyz[1]+(c1*s2)*xyz[2]
    y = (+c1*s3+c2*c3*s1)*xyz[0]+(c1*c3-c2*s1*s3)*xyz[1]+(s1*s2)*xyz[2]
    z = (-c3*s2)*xyz[0]+(s2*s3)*xyz[1]+(c2)*xyz[2]
    return x, y, z

def str2func(function):
    if isinstance(function, str):
        if function.lower() == 'gaussian':
            # p[0] is normalization
            # p[1] is mean (first raw moment)
            # p[2] is sigma (square root of second central moment)
            f = lambda p, x: p[0]/np.sqrt(2.*np.pi)/p[2]*np.exp(-((x-p[1])/p[2])**2/2.)
        else:
            raise ValueError
    return f

def fitfunc(x, y, function, initial, free=None, yerr=None, **kwargs):
    """Wrapper to scipy.optimize.leastsq to fit data to an arbitrary function."""
    f = str2func(function) if isinstance(function, str) else function
    if free is None: free = np.ones(np.shape(initial))
    if yerr is None: yerr = np.ones(np.shape(y))

    errfunc = lambda p, xp, yp, yerrp: (yp-f(p*free+initial*np.logical_not(free), xp))/yerrp

    return optimize.leastsq(errfunc, initial, args=(x, y, yerr), **kwargs)

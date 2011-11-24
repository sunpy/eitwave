"""
Simulates a wave
"""

from __future__ import absolute_import

__all__ = ["simulate", "simulate_raw", "transform", "add_noise"]

__authors__ = ["Albert Shih"]
__email__ = "albert.y.shih@nasa.gov"

import numpy as np
import sunpy

def prep_coeff(coeff, order=2):
    """
    Prepares polynomial coefficients out to a certain order, outputs as ndarray
    """
    new_coeff = np.zeros(order+1)
    if type(coeff) == list or type(coeff) == np.ndarray:
        size = min(len(coeff), len(new_coeff))
        new_coeff[0:size] = coeff[0:size]
    else:
        new_coeff[0] = coeff
    return new_coeff

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

def simulate_raw(params, verbose = False):
    """
    Simulate data in HG' coordinates
    
    HG' = HG, except center at wave epicenter
    """
    from sunpy.util import util
    import datetime
    from scipy.special import ndtr

    cadence = params["cadence"]
    direction = 180.+params["direction"]
    
    width_coeff = prep_coeff(params["width"])
    wave_thickness_coeff = prep_coeff(params["wave_thickness"])
    wave_normalization_coeff = prep_coeff(params["wave_normalization"])
    speed_coeff = prep_coeff(params["speed"])
    
    steps = params["max_steps"]
    
    lat_min = params["lat_min"]
    lat_max = params["lat_max"]
    lat_bin = params["lat_bin"]
    lon_min = params["lon_min"]
    lon_max = params["lon_max"]
    lon_bin = params["lon_bin"]

    #This roundabout approach recalculates lat_bin and lon_bin to produce
    #equally sized bins to exactly span the min/max ranges
    lat_num = int(round((lat_max-lat_min)/lat_bin))
    lat_edges, lat_bin = np.linspace(lat_min, lat_max, lat_num+1,
                                     retstep=True)
    lon_num = int(round((lon_max-lon_min)/lon_bin))
    lon_edges, lon_bin = np.linspace(lon_min, lon_max, lon_num+1,
                                     retstep=True)
    
    #Propagates from 90. down to lat_min, irrespective of lat_max
    p = np.poly1d([speed_coeff[2]/3., speed_coeff[1]/2., speed_coeff[0],
                   -(90.-lat_min)])
    
    #Will fail if wave does not propogate all the way to lat_min
    #duration = p.r[np.logical_and(p.r.real > 0, p.r.imag == 0)][0]
    
    #steps = int(duration/cadence)+1
    #if steps > params["max_steps"]:
    #    steps = params["max_steps"]
    
    #Maybe used np.poly1d() instead to do the polynomial calculation?
    time = np.arange(steps)*cadence
    time_powers = np.vstack((time**0, time**1, time**2))
    
    width = np.dot(width_coeff, time_powers).ravel()
    wave_thickness = np.dot(wave_thickness_coeff, time_powers).ravel()
    wave_normalization = np.dot(wave_normalization_coeff, time_powers).ravel()
    
    #Propagates from 90., irrespective of lat_max
    wave_peak = 90.-(p(time)+(90.-lat_min))
    
    out_of_bounds = np.logical_or(wave_peak < lat_min, wave_peak > lat_max)
    if out_of_bounds.any():
        steps = np.where(out_of_bounds)[0][0]
        
    wave_maps = []
    
    dict_header = {
        "cdelt1": lon_bin,
        "naxis1": lon_num,
        "crval1": lon_min,
        "crpix1": 0.5, #this makes lon_min the left edge of the first bin
        "cunit1": "deg",
        "ctype1": "HG",
        "cdelt2": lat_bin,
        "naxis2": lat_num,
        "crval2": lat_min,
        "crpix2": 0.5, #this makes lat_min the left edge of the first bin
        "cunit2": "deg",
        "ctype2": "HG",
        "hglt_obs": 0,
        "hgln_obs": 0,
        "rsun_obs": 963.879683,
        "rsun_ref": 696000000.0,
        "dsun_obs": 148940609626.98
    }
    
    header = sunpy.map.MapHeader(dict_header)
    
    if verbose:
        print("Simulating "+str(steps)+" raw maps")
    
    for istep in xrange(steps):
        #Gaussian profile in longitudinal direction
        #Does not take into account spherical geometry
        #(i.e., change in area element)
        if (wave_thickness[istep] <= 0):
            print("ERROR: wave thickness is non-physical!")
        z = (lat_edges-wave_peak[istep])/wave_thickness[istep]
        wave_1d = wave_normalization[istep]*(ndtr(np.roll(z, -1))-ndtr(z))[0:lat_num]
        wave_1d /= lat_bin
        
        wave_lon_min = direction-width[istep]/2
        wave_lon_max = direction+width[istep]/2

        if (width[istep]< 360.):
            #Do these need to be np.remainder() instead?
            wave_lon_min_mod = ((wave_lon_min+180.) % 360.)-180.
            wave_lon_max_mod = ((wave_lon_max+180.) % 360.)-180.
            
            index1 = np.arange(lon_num+1)[np.roll(lon_edges, -1) > min(wave_lon_min_mod, wave_lon_max_mod)][0]
            index2 = np.roll(np.arange(lon_num+1)[lon_edges < max(wave_lon_min_mod, wave_lon_max_mod)], 1)[0]
    
            wave_lon = np.zeros(lon_num)
            wave_lon[index1+1:index2] = 1.
            #Possible weirdness if index1 == index2
            wave_lon[index1] += (lon_edges[index1+1]-min(wave_lon_min_mod, wave_lon_max_mod))/lon_bin
            wave_lon[index2] += (max(wave_lon_min_mod, wave_lon_max_mod)-lon_edges[index2])/lon_bin
            
            if (wave_lon_min_mod > wave_lon_max_mod):
                wave_lon = 1.-wave_lon
        else:
            wave_lon = np.ones(lon_num)
        
        #Could be accomplished with np.dot() without casting as matrices?
        wave = np.mat(wave_1d).T*np.mat(wave_lon)
        
        wave_maps += [sunpy.map.BaseMap(wave, header)]
        wave_maps[istep].name = "Simulation"
        wave_maps[istep].date = util.anytim("2011-11-11")+datetime.timedelta(0, istep*cadence)
    
    return wave_maps

def transform(params, wave_maps, verbose = False):
    """
    Transform raw data in HG' coordinates to HPC coordinates
    
    HG' = HG, except center at wave epicenter
    """
    from scipy.interpolate import griddata
    
    hglt_obs = params["hglt_obs"]
    rotation = params["rotation"]
    
    epi_lat = params["epi_lat"]
    epi_lon = params["epi_lon"]
    
    hpcx_min = params["hpcx_min"]
    hpcx_max = params["hpcx_max"]
    hpcx_bin = params["hpcx_bin"]
    hpcy_min = params["hpcy_min"]
    hpcy_max = params["hpcy_max"]
    hpcy_bin = params["hpcy_bin"]
    
    hpcx_num = int(round((hpcx_max-hpcx_min)/hpcx_bin))
    hpcy_num = int(round((hpcy_max-hpcy_min)/hpcy_bin))
    
    wave_maps_transformed = []
    
    dict_header = {
        "cdelt1": hpcx_bin,
        "naxis1": hpcx_num,
        "crval1": hpcx_min,
        "crpix1": 0.5, #this makes hpcx_min the left edge of the first bin
        "cunit1": "arcsec",
        "ctype1": "HPC",
        "cdelt2": hpcy_bin,
        "naxis2": hpcy_num,
        "crval2": hpcy_min,
        "crpix2": 0.5, #this makes hpcy_min the left edge of the first bin
        "cunit2": "arcsec",
        "ctype2": "HPC",
        "hglt_obs": hglt_obs,
        "hgln_obs": 0,
        "rsun_obs": 963.879683,
        "rsun_ref": 696000000.0,
        "dsun_obs": 148940609626.98
    }
    
    header = sunpy.map.MapHeader(dict_header)
    
    start_date = wave_maps[0].date
    
    #Origin grid, HG'
    lon_grid, lat_grid = sunpy.wcs.convert_pixel_to_data(wave_maps[0].header)
    
    #Origin grid, HG' to HCC'
    #HCC' = HCC, except centered at wave epicenter
    x, y, z = sunpy.wcs.convert_hg_hcc_xyz(wave_maps[0].header,
                                           lon_grid, lat_grid)
    
    #Origin grid, HCC' to HCC''
    #Moves the wave epicenter to initial conditions
    #HCC'' = HCC, except assuming that HGLT_OBS = 0
    zxy_p = euler_zyz((z, x, y), (epi_lon, 90.-epi_lat, 0.))
    
    #Destination HPC grid
    hpcx_grid, hpcy_grid = sunpy.wcs.convert_pixel_to_data(dict_header)
    
    for current_wave_map in wave_maps:
        if verbose:
            print("Transforming map at "+str(current_wave_map.date))
        
        #Origin grid, HCC'' to HCC
        #Moves the observer to HGLT_OBS and adds rigid solar rotation
        td = current_wave_map.date-start_date
        total_seconds = (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6
        zpp, xpp, ypp = euler_zyz(zxy_p, (0., hglt_obs, total_seconds*rotation))
        
        #Origin grid, HCC to HPC (arcsec)
        xx, yy = sunpy.wcs.convert_hcc_hpc(current_wave_map.header, xpp, ypp)
        xx *= 3600
        yy *= 3600
        
        #Coordinate positions (HPC) with corresponding map data
        points = np.vstack((xx.ravel(), yy.ravel())).T
        values = np.array(current_wave_map).ravel()
        
        #2D interpolation from origin grid to destination grid
        grid = griddata(points[zpp.ravel() >= 0], values[zpp.ravel() >= 0],
                        (hpcx_grid, hpcy_grid), method="linear")
        
        transformed_wave_map = sunpy.map.BaseMap(grid, header)
        transformed_wave_map.name = current_wave_map.name
        transformed_wave_map.date = current_wave_map.date
        wave_maps_transformed += [transformed_wave_map]

    return wave_maps_transformed

def add_noise(params, wave_maps, verbose = False):
    """
    Adds simulated noise to a map
    """
    from scipy import stats
    
    noise_type = params["noise_type"]
    noise_scale = params["noise_scale"]
    noise_mean = params["noise_mean"]
    noise_sdev = params["noise_sdev"]
    
    if noise_type is None:
        if verbose:
            print("No added noise")
        return wave_maps
    else:
        wave_maps_noise = []
        for current_wave_map in wave_maps:
            wave = np.array(current_wave_map)
            
            if noise_type == "Normal":
                if verbose:
                    print("Adding normal noise to map at "+str(current_wave_map.date))
                noise = noise_scale*stats.norm.rvs(loc=noise_mean, scale=noise_sdev, size=wave.size).reshape(wave.shape)
            elif noise_type == "Poisson":
                if verbose:
                    print("Adding Poisson noise to map at "+str(current_wave_map.date))
                noise = noise_scale*stats.poisson.rvs(noise_mean, size=wave.size).reshape(wave.shape)
            else:
                if verbose:
                    print("Unknown noise requested to map at "+str(current_wave_map.date))
                noise = np.zeros_like(wave)
            
            new_wave_map = sunpy.map.BaseMap(wave+noise,
                                             current_wave_map.header)
            new_wave_map.name = current_wave_map.name
            new_wave_map.date = current_wave_map.date
            wave_maps_noise += [new_wave_map]
            
        return wave_maps_noise

def simulate(params, verbose = False):
    """
    Simulates wave in HPC coordinates with added noise
    """
    wave_maps_raw = simulate_raw(params, verbose)
    wave_maps_transformed = transform(params, wave_maps_raw, verbose)
    wave_maps_noise = add_noise(params, wave_maps_transformed, verbose)
    return wave_maps_noise

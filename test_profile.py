import numpy as np
import matplotlib.pyplot as plt

from sim import wave2d
from visualize import visualize

import util

m2deg = 360./(2*3.1415926*6.96e8)

params = {
    "cadence": 12., #seconds
    
    "hglt_obs": 0., #degrees
#    "rotation": 360./(27.*86400.), #degrees/s, rigid solar rotation
    "rotation": 0., #degrees/s, rigid solar rotation
    
    #Wave parameters that are initial conditions
    "direction": 25., #degrees, measured CCW from HG +latitude
    "epi_lat": 0., #degrees, HG latitude of wave epicenter
    "epi_lon": 0., #degrees, HG longitude of wave epicenter
    
    #Wave parameters that can evolve over time
    #The first element is constant in time
    #The second element (if present) is linear in time
    #The third element (if present) is quadratic in time
    #Be very careful of non-physical behavior
    "width": [90., 1.5], #degrees, full angle in azimuth, centered at 'direction'
    "wave_thickness": [6.0e6*m2deg,6.0e4*m2deg], #degrees, sigma of Gaussian profile in longitudinal direction
    "wave_normalization": [1.], #integrated value of the 1D Gaussian profile
    "speed": [9.33e5*m2deg, -1.495e3*m2deg], #degrees/s, make sure that wave propagates all the way to lat_min for polynomial speed
    
    #Random noise parameters
    "noise_type": None, #can be None, "Normal", or "Poisson"
    "noise_scale": 0.1,
    "noise_mean": 1.,
    "noise_sdev": 1.,
    
    #Structured noise parameters
    "struct_type": None, #can be None, "Arcs", or "Random"
    "struct_scale": 5.,
    "struct_num": 10,
    "struct_seed": 13092,
    
    "max_steps": 20,
    
    "clean_nans": True,
    
    #HG grid, probably would only want to change the bin sizes
    "lat_min": -90.,
    "lat_max": 90.,
    "lat_bin": 0.2,
    "lon_min": -180.,
    "lon_max": 180.,
    "lon_bin": 5.,
    
    #HPC grid, probably would only want to change the bin sizes
    "hpcx_min": -1228.8,
    "hpcx_max": 1228.8,
    "hpcx_bin": 2.4,
    "hpcy_min": -1228.8,
    "hpcy_max": 1228.8,
    "hpcy_bin": 2.4
}

#wave_maps = wave2d.simulate(params)
wave_maps = wave2d.simulate(params, verbose = True)

#To get simulated HG' maps (centered at wave epicenter):
wave_maps_raw = wave2d.simulate_raw(params)
#wave_maps_raw_noise = wave2d.add_noise(params, wave_maps_raw)

visualize(wave_maps)


new_wave_maps = []

for wave in wave_maps:
    print("Unraveling map at "+str(wave.date))
    new_wave_maps += [util.map_hpc_to_hg_rotate(wave,
                                                epi_lon = params["epi_lon"],
                                                epi_lat = params["epi_lat"],
                                                xbin = 5, ybin = 0.2)]

lat_min = params["lat_min"]
cadence = params["cadence"]
width_coeff = wave2d.wave2d.prep_coeff(params["width"])
wave_thickness_coeff = wave2d.wave2d.prep_coeff(params["wave_thickness"])
wave_normalization_coeff = wave2d.wave2d.prep_coeff(params["wave_normalization"])
speed_coeff = wave2d.wave2d.prep_coeff(params["speed"])
p = np.poly1d([speed_coeff[2]/3., speed_coeff[1]/2., speed_coeff[0],
               -(90.-lat_min)])
time = np.arange(len(wave_maps))*cadence
time_powers = np.vstack((time**0, time**1, time**2))

width = np.dot(width_coeff, time_powers).ravel()
wave_thickness = np.dot(wave_thickness_coeff, time_powers).ravel()
wave_normalization = np.dot(wave_normalization_coeff, time_powers).ravel()
wave_peak = 90.-(p(time)+(90.-lat_min))

n0 = wave_normalization*width.clip(0,360)/params["lon_bin"]
m0 = wave_peak
s0 = wave_thickness

n1 = []
m1 = []
s1 = []

for wave in wave_maps_raw:
    yy = np.sum(np.asarray(wave), axis=1)
    xx = np.arange(wave.yrange[0],wave.yrange[1],wave.scale['y'])+wave.scale['y']/2.
    p, success = util.fitfunc(xx, yy, 'gaussian', [1, 90, 1])
    if p[0] < 0.1:
        p, success = util.fitfunc(xx, yy, 'gaussian', [1, 75, 1])       
    n1 += [p[0]]
    m1 += [p[1]]
    s1 += [p[2]]

n2 = []
m2 = []
s2 = []

for wave in new_wave_maps:
    data = np.array(wave)
    data[np.isnan(data)] = 0.
    yy = np.sum(data, axis=1)
    xx = np.linspace(wave.yrange[0],wave.yrange[1],wave.shape[0])+wave.scale['y']/2.
    p, success = util.fitfunc(xx, yy, 'gaussian', [1, 90, 1])
    if p[0] < 0.1:
        p, success = util.fitfunc(xx, yy, 'gaussian', [1, 75, 1])       
    n2 += [p[0]]
    m2 += [p[1]]
    s2 += [p[2]]

plt.figure()
plt.plot(n0)
plt.plot(n1,'+')
plt.plot(n2,'x')

plt.figure()
plt.plot(m0)
plt.plot(m1,'+')
plt.plot(m2,'x')

plt.figure()
plt.plot(s0)
plt.plot(s1,'+')
plt.plot(s2,'x')

xx0 = np.linspace(params["lat_min"],params["lat_max"],num=10000,endpoint=True)

wave = wave_maps_raw[3]
data1 = np.array(wave)
yy1 = np.sum(data1, axis=1)
xx1 = np.arange(wave.yrange[0],wave.yrange[1],wave.scale['y'])+wave.scale['y']/2.

wave = new_wave_maps[3]
data2 = np.array(wave)
data2[np.isnan(data2)] = 0.
yy2 = np.sum(data2, axis=1)
xx2 = np.linspace(wave.yrange[0],wave.yrange[1],wave.shape[0])+wave.scale['y']/2.


plt.figure()
plt.plot(xx0,util.str2func('gaussian')([n0[3],m0[3],s0[3]],xx0))
plt.plot(xx1,yy1,'+')
plt.plot(xx2,yy2,'x')

"""
from matplotlib import colors

wave_maps_raw = wave2d.simulate_raw(params)
wave_maps_transformed = wave2d.transform(params, wave_maps_raw, verbose = True)

#First simulation slide
wave_maps_raw[19].show()
wave_maps_transformed[19].show()

#Second simulation slide
wave_maps[19].show(norm = colors.Normalize(0,1))
new_wave_maps[19].show(norm = colors.Normalize(0,1))
"""
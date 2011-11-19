from sim import wave2d
from visualize import visualize

params = {
    "cadence": 12., #seconds
    
    "direction": 25., #degrees, measured CCW from HPC +X
    "epi_lat": 30., #degrees, HG latitude of wave epicenter
    "epi_lon": 20., #degrees, HG longitude of wave epicenter
    
    #Parameters that evolve over time, be very careful of non-physical behavior
    "width": [90., 1.], #degrees, full angle in azimuth, centered at 'direction'
    "wave_thickness": [1., 0.02, -0.00005], #degrees, sigma of Gaussian profile in longitudinal direction
    "wave_normalization": [1.], #integrated value of the 1D Gaussian profile
    "speed": [0.2, 0.002], #degrees/s, make sure that wave propagates all the way to lat_min for polynomial speed
    
    "max_steps": 50,
    
    #HG grid, probably would only want to change the bin sizes
    "lat_min": -90.,
    "lat_max": 90.,
    "lat_bin": 1.,
    "lon_min": -180.,
    "lon_max": 180.,
    "lon_bin": 1.,
    
    #HPC grid, probably would only want to change the bin sizes
    "hpcx_min": -1000.,
    "hpcx_max": 1000.,
    "hpcx_bin": 10.,
    "hpcy_min": -1000.,
    "hpcy_max": 1000.,
    "hpcy_bin": 10.
}

wave_maps = wave2d.simulate(params)
visualize(wave_maps)

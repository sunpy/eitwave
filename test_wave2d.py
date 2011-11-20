from sim import wave2d
from visualize import visualize

params = {
    "cadence": 12., #seconds
    
    "hglt_obs": 0., #degrees
    "rotation": 360./(27.*86400.), #degrees/s, rigid solar rotation
    
    #Wave parameters that are initial conditions
    "direction": 25., #degrees, measured CCW from HG +longitude
    "epi_lat": 30., #degrees, HG latitude of wave epicenter
    "epi_lon": 20., #degrees, HG longitude of wave epicenter
    
    #Wave parameters that can evolve over time
    #The first element is constant in time
    #The second element (if present) is linear in time
    #The third element (if present) is quadratic in time
    #Be very careful of non-physical behavior
    "width": [90., 1.5], #degrees, full angle in azimuth, centered at 'direction'
    "wave_thickness": [0.5, 0.02, -0.00005], #degrees, sigma of Gaussian profile in longitudinal direction
    "wave_normalization": [1.], #integrated value of the 1D Gaussian profile
    "speed": [0.2, 0.002], #degrees/s, make sure that wave propagates all the way to lat_min for polynomial speed
    
    #Noise parameters
    "noise_type": "Poisson", #can be None, "Normal", or "Poisson"
    "noise_scale": 0.05,
    "noise_mean": 1.,
    "noise_sdev": 1.,
    
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
#wave_maps = wave2d.simulate(params, verbose = True)
visualize(wave_maps)

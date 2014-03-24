from __future__ import absolute_import

__authors__ = ["Albert Shih"]
__email__ = "albert.y.shih@nasa.gov"

def visualize(wave_maps, delay = 0.1, range = None):
    """
    Visualizes a list of SunPy Maps.  Requires matplotlib 1.1
    """
    import matplotlib.pyplot as plt
    from matplotlib import colors

    fig = plt.figure()
    
    axes = fig.add_subplot(111)
    axes.set_xlabel('X-position [' + wave_maps[0].units['x'] + ']')
    axes.set_ylabel('Y-position [' + wave_maps[0].units['y'] + ']')
    
    extent = wave_maps[0].xrange + wave_maps[0].yrange
    axes.set_title("%s %s" % (wave_maps[0].name, wave_maps[0].date))
    params = {
        "cmap": wave_maps[0].cmap,
        "norm": wave_maps[0].mpl_color_normalizer
    }
    if range != None:
        params["norm"] = colors.Normalize(range[0],range[1])
    img = axes.imshow(wave_maps[0], origin='lower', extent=extent, **params)
    fig.colorbar(img)
    fig.show()
    
    for current_wave_map in wave_maps[1:]:
        axes.set_title("%s %s" % (current_wave_map.name, current_wave_map.date))
        img.set_data(current_wave_map)
        plt.pause(delay)

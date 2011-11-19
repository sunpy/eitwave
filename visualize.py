def visualize(wave_maps,delay=0.1):
    import matplotlib.pyplot as plt

    fig = plt.figure()
    
    axes = fig.add_subplot(111)
    axes.set_xlabel('X-position [' + wave_maps[0].units['x'] + ']')
    axes.set_ylabel('Y-position [' + wave_maps[0].units['y'] + ']')
    
    extent = wave_maps[0].xrange + wave_maps[0].yrange
    axes.set_title("%s %s" % (wave_maps[0].name, wave_maps[0].date))
    params = {
        "cmap": wave_maps[0].cmap,
        "norm": wave_maps[0].norm()
    }
    im = axes.imshow(wave_maps[0], origin='lower', extent=extent, **params)
    fig.colorbar(im)
    fig.show()
    
    for current_wave_map in wave_maps[1:]:
        axes.set_title("%s %s" % (current_wave_map.name, current_wave_map.date))
        im.set_data(current_wave_map)
        plt.pause(delay)

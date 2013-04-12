
from visualize import visualize
import eitwaveutils
import copy
import os
from test_wave2d import test_wave2d
from sunpy.time import TimeRange, parse_time
from sunpy.net import hek


def main(source_data='.jp2',
         time_range=TimeRange('2011/10/01 09:45:00', '2011/10/01 10:15:59'),
         algorithm='hough', feed_directory='~/Data/eitwave/jp2/20111001/'):
    '''
    source_data { jp2 | fits | test }
    look for helioviewer JP2 files, FITS files, or load the test data
    respectively

    time_range : a TimeRange object
    time range that is searched for EUV waves

    feed_directory
    If set to a string, look in this directory for the jp2 files.  Assumes that
    all the JP2 files are relevant to the flare detection.

    algorithm: { 'hough' : 'phough' }
    algorithm used to find the wave
    '''

    if source_data == 'test':
        # where to store those data
        maps = test_wave2d()
    elif source_data == '.jp2':
        # where to store those data
        data_storage = "~/Data/eitwave/jp2/AGU/"

    # Query the HEK for flare information we need
    client = hek.HEKClient()
    hek_result = client.query(hek.attrs.Time(time_range.t1, time_range.t2),
                              hek.attrs.EventType('FL'))
                              #hek.attrs.FRM.Name == '')
    # no flares, no analysis possible
    if hek_result is None:
        return None

    # Flares!
    print('Number of flares found = ' + str(len(hek_result)))

    for flare in hek_result[10:11]:

        if feed_directory is None:
            print('Acquiring data for flare')
            filelist = eitwaveutils.acquire_data(data_storage, source_data,
                                                 flare)
        else:
            # Assumes that the necessary files are already present
            filelist = eitwaveutils.listdir_fullpath(feed_directory)

        # reduce the number of files to those that happen after the flare has
        # started
        files = []
        for f in filelist:
            fhv = f.split(os.sep)[-1]
            if eitwaveutils.hv_filename2datetime(fhv) > \
            parse_time(flare['event_starttime']):
                files.append(f)
        print('Number of files :' + str(len(files)))
        if len(files) == 0:
            print('No files found.  Returning.')
            return None

        # Define the transform parameters
        # params = eitwaveutils.params(flare='test')
        params = eitwaveutils.params(flare)

        # read in files and accumulate them
        maps = eitwaveutils.accumulate(files, accum=1, nsuper=4,
                                   verbose=True)

        # Unravel the maps
        new_maps = eitwaveutils.map_unravel(maps, params, verbose=True)

        #sometimes unravelling maps leads to slight variations in the unraveled
        #image dimensions.  check dimensions of maps and resample to dimensions
        #of first image in sequence if need be.
        new_maps = eitwaveutils.check_dims(new_maps)

        # calculate the differences
        diffs = eitwaveutils.map_diff(new_maps)

        #determine the threshold to apply to the difference maps.
        #diffs > diff_thresh will be True, otherwise False.
        threshold_maps = eitwaveutils.map_threshold(new_maps, factor=0.2)

        # transform difference maps into binary maps
        binary_maps = eitwaveutils.map_binary(diffs, threshold_maps)

        if algorithm == 'hough':
            # detection based on the hough transform
            detection = eitwaveutils.hough_detect(binary_maps, vote_thresh=10)
        elif algorithm == 'prob_hough':
            # detection based on the probabilistic hough transform.  Takes the
            # keywords of the probabilistic hough transform - see the documentation
            # of skimage.transform.probabilistic_hough (scikit-image.org)
            detection = eitwaveutils.prob_hough_detect(binary_maps, threshold=10)

        # Remove areas that are too small or that don't have enough detections
        detection = eitwaveutils.cleanup(detection,
                                         size_thresh=50,
                                         inv_thresh=8)

        #If there is anything left in 'detection', fit a function to the original
        #diffmaps in the region defined by 'detection'. Simplest case: fit a
        #Gaussian in the y-direction for some x or range of x.
        #eitwaveutils.fit_wavefront should probably take the arguments of fitfunc.
        #use 'detection' to guess starting fit parameters?

        #get just the positive elements of the difference map. Perform fitting on
        #these positive diffmaps.
        posdiffs = copy.deepcopy(diffs)
        for i in range(0, len(diffs)):
            temp = diffs[i] < 0
            posdiffs[i][temp] = 0

        #fit a function to the difference maps in the cases where there has been a
        #detection
        wavefront = eitwaveutils.fit_wavefront(posdiffs, detection)

        #strip out the velocity information from the wavefront fitting
        velocity = eitwaveutils.wavefront_velocity(wavefront[0])

        #strip out the position and width information from the wavefront fitting
        pos_width = eitwaveutils.wavefront_position_and_width(wavefront[0])

        visualize(detection)

    return maps, new_maps, diffs, threshold_maps, binary_maps
    #, detection,wavefront, velocity, pos_width

if __name__ == '__main__':
    main()

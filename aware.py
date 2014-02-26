
from visualize import visualize
import aware_utils
import copy
import os
import cPickle as pickle
from test_wave2d import test_wave2d
from sunpy.time import TimeRange, parse_time
from sunpy.net import hek
import os


def main(source_data='.jp2',
         time_range=TimeRange('2011/10/01 09:45:00', '2011/10/01 10:15:59'),
         algorithm='hough', feed_directory='~/Data/eitwave/jp2/20111001_jp2/',
         use_pickle=None,diff_type='running'):
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
    if feed_directory != None:
        feed_directory = os.path.expanduser(feed_directory)

    if source_data == 'test':
        # where to store those data
        maps = test_wave2d()
    elif source_data == '.jp2':
        # where to store those data
        data_storage = "~/Data/eitwave/jp2/20120607/"

    if not os.path.exists(os.path.expanduser(data_storage)):
            os.makedirs(os.path.expanduser(data_storage))

    # Query the HEK for flare information we need
    client = hek.HEKClient()
    hek_result = client.query(hek.attrs.Time(time_range.t1, time_range.t2),
                              hek.attrs.EventType('FL'),hek.attrs.OBS.ChannelID == '211')
    #hek.attrs.FRM.Name == '')
    if hek_result is None:
    # no flares, no analysis possible
        return None

    # Flares!
    print('Number of flares found = ' + str(len(hek_result)))

    for flare in hek_result[0:1]:

        if feed_directory is None:
            print('Acquiring data for flare')
            filelist = aware_utils.acquire_data(data_storage, source_data,
                                                 flare)
        else:
            # Assumes that the necessary files are already present

            filelist = aware_utils.listdir_fullpath(feed_directory,
                                                     filetype ='jp2')

        #filter to only grab the data files with the source_data extn in the directory
        files_tmp = []
        for f in filelist:
            if f.endswith(source_data):
                files_tmp.append(f)
            files = files_tmp

        # reduce the number of files to those that happen after the flare has
        # started
        files = []
        for f in files_tmp:
            fhv = f.split(os.sep)[-1]
            if aware_utils.hv_filename2datetime(fhv) > \
            parse_time(flare['event_starttime']):
                files.append(f)
        print('Number of files :' + str(len(files)))
        if len(files) == 0:
            print('No files found.  Returning.')
            return None

        # Define the transform parameters
        params = aware_utils.params(flare)

        # read in files and accumulate them
        if use_pickle != None:
            # load in a pickle file of the data
            pfile = open(feed_directory + use_pickle, 'rb')
            a = pickle.load(pfile)
            maps = a[0]
            new_maps = a[1]
            diffs = a[2]
            pfile.close()
        else:
            maps = aware_utils.accumulate(files[6:30], accum=1, nsuper=4,
                                   verbose=True)

            #temporary fix for exposure control and S/N changes
            long_maps = []
            for m in maps:
                if m.exposure_time > 2.0:
                    long_maps.append(m)
            maps=long_maps

            # Unravel the maps
            new_maps = aware_utils.map_unravel(maps, params, verbose=True)
            #return new_maps

            #sometimes unravelling maps leads to slight variations in the unraveled
            #image dimensions.  check dimensions of maps and resample to dimensions
            #of first image in sequence if need be.
            #new_maps[0].peek()
            new_maps = aware_utils.check_dims(new_maps)

            # calculate the differences
            if diff_type == 'base':
                diffs=aware_utils.map_basediff(new_maps)
            else:
                diffs = aware_utils.map_diff(new_maps)


        #generate persistence maps - currently bugged, so skip this step
        #persistence_maps = eitwaveutils.map_persistence(diffs)
        persistence_maps=[]

        #determine the threshold to apply to the difference maps.
        #diffs > diff_thresh will be True, otherwise False.
        threshold_maps = aware_utils.map_threshold(new_maps, factor=0.2)
        #return threshold_maps

        # transform difference maps into binary maps
        binary_maps = aware_utils.map_binary(diffs, threshold_maps)

        if algorithm == 'hough':
            # detection based on the hough transform
            detection = aware_utils.hough_detect(binary_maps, vote_thresh=10)
        elif algorithm == 'prob_hough':
            # detection based on the probabilistic hough transform.  Takes the
            # keywords of the probabilistic hough transform - see the documentation
            # of skimage.transform.probabilistic_hough (scikit-image.org)
            detection = aware_utils.prob_hough_detect(binary_maps, threshold=10)

        # Remove areas that are too small or that don't have enough detections
        detection = aware_utils.cleanup(detection,
                                         size_thresh=50,
                                         inv_thresh=8)

        detection_maps = copy.deepcopy(binary_maps)
        for i in range(0,len(detection)):
            detection_maps[i].data = detection[i]
        #If there is anything left in 'detection', fit a function to the original
        #diffmaps in the region defined by 'detection'. Simplest case: fit a
        #Gaussian in the y-direction for some x or range of x.
        #eitwaveutils.fit_wavefront should probably take the arguments of fitfunc.
        #use 'detection' to guess starting fit parameters.

        #get just the positive elements of the difference map. Perform fitting on
        #these positive diffmaps.
        posdiffs = copy.deepcopy(diffs)
        for i in range(0, len(diffs)):
            temp = diffs[i].data < 0
            posdiffs[i].data[temp] = 0

        #fit a function to the difference maps in the cases where there has been a
        #detection
        wavefront = aware_utils.fit_wavefront(posdiffs, detection)

        #transform the detected model wavefront back into heliocentric coordinates so it can be overlayed
        wavefront_hc = aware_utils.map_reravel(wavefront[1],params,verbose=True)

        #strip out the velocity information from the wavefront fitting
        velocity = aware_utils.wavefront_velocity(wavefront[0])

        #strip out the position and width information from the wavefront fitting
        pos_width = aware_utils.wavefront_position_and_width(wavefront[0])

        visualize(detection_maps)

    return maps, new_maps, diffs, threshold_maps, binary_maps, detection, wavefront, velocity, pos_width, persistence_maps, wavefront_hc,params


if __name__ == '__main__':
    main()


import copy
import os
import cPickle as pickle

from sunpy.time import TimeRange, parse_time
from sunpy.net import hek

from test_wave2d import test_wave2d
from visualize import visualize
import aware_utils


def main(source_data='.jp2',
         time_range=TimeRange('2011/10/01 09:45:00', '2011/10/01 10:15:59'),
         algorithm='hough', feed_directory='~/Data/eitwave/jp2/20111001_jp2/',
         use_pickle=None,diff_type='running',data_savedir=None):
    '''
    This is the main executable for the Automated EUV Wave Analysis and Reduction (AWARE)
    code. The procedure is as follows:
        - Query the HEK to find whether any flares were observed in SDO/AIA 211A during the input time range
        - If yes, then read in (or download) a sequence of solar images corresponding to the time range
        - Transform these images from Helioprojective Coordinates to Heliographic Coordinates,
          with the origin centered at the flare origin
        - Create a sequence of difference maps from these transformed maps
        - Use a threshold method to create a binary map from the the difference maps.
        - Apply the Hough Transform to the binary map to search for strong lines in the image
        - Use the results of the Hough Transform to detect whether an EUV wave is present
        - Fit an appropriate function (e.g. Gaussian) to the detected wavefront as a function of longitude
        - Record fit results and return data products
    
    Parameters
    ----------    
    source_data : string
        description of the type of data being input. Allowed is '.jp2', 'fits', or 'test'.
        will look for helioviewer JP2 files, FITS files, or load the test data respectively

    time_range : a TimeRange object
        time range within which to search for EUV waves

    feed_directory : string
        A directory containing data files to be analysed. If set, AWARE will assume data files are
        already download and will search in this directory instead. Assumes that all files in the
        directory with the appropriate extension (e.g. .jp2, .fits) are relevant to the flare detection.

    use_pickle : string
        BUGGED - currently not supported, always set to None

    diff_type : string
        The type of image differencing to use. Allowed values are 'running' or 'base'. Default is 'running'
        Will perform either running differencing or base differencing on the image sequence.

    data_savedir : string
        directory in which to save downloaded jp2 files from Helioviewer. If None, then AWARE will construct a directory
        based on the start time of the query.

    Returns
    -------

    Outputs a pickle file containing the following data products (in order):
        1) a list of maps modelling the detected wavefront, transformed back to original HPC coordinates  
    '''


    if feed_directory != None:
        feed_directory = os.path.expanduser(feed_directory)

    #Check which type of data is being analysed, and establish the directory to store downloaded files,
    #if appropriate
    if source_data == 'test':
        maps = test_wave2d()
    elif source_data == '.jp2' and data_savedir == None:
        data_savedir = '~/aware_data/jp2/' + time_range.start().strftime('%Y%m%d_%H%M')
    elif source_data == '.fits' and data_savedir == None:
        data_savedir = '~/aware_data/fits/' + time_range.start().strftime('%Y%m%d_%H%M')
            
    if not os.path.exists(os.path.expanduser(data_savedir)):
        os.makedirs(os.path.expanduser(data_savedir))

    # Query the HEK to see whether there were any flares during the time range specified
    # Concentrate on the AIA 211A channel as it has clearest observations of global waves
    client = hek.HEKClient()
    hek_result = client.query(hek.attrs.Time(time_range.t1, time_range.t2),
                              hek.attrs.EventType('FL'),hek.attrs.OBS.ChannelID == '211')
    if hek_result is None:
    # if no flares found, no analysis possible. Return
        print 'No flares found in HEK database during specified time range.'
        print 'No analysis possible. Returning.'
        return None

    # Otherwise, we have found at least one flare
    print('Number of flares found = ' + str(len(hek_result)))

    #assume the first result of the HEK query has the correct information
    for flare in hek_result[0:1]:

        if feed_directory is None:
            print('Acquiring data for flare')
            filelist = aware_utils.acquire_data(data_savedir, source_data,
                                                 flare)
        else:
            # Assumes that the necessary files are already present
            filelist = aware_utils.listdir_fullpath(feed_directory,
                                                     filetype = source_data)

        #filter to only grab the data files with the source_data extn in the directory
        #this looks like duplication of listdir_fullpath
        files_tmp = []
        for f in filelist:
            if f.endswith(source_data):
                files_tmp.append(f)
            files = files_tmp

        # reduce the number of files to those that happen after the flare has
        # started
        files = []
        if source_data == '.jp2':
            for f in files_tmp:
                fhv = f.split(os.sep)[-1]
                if aware_utils.hv_filename2datetime(fhv) > \
                parse_time(flare['event_starttime']):
                    files.append(f)
            print('Number of files :' + str(len(files)))
            if len(files) == 0:
                print('No files found.  Returning.')
                return None
        else:
            files = files_tmp

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
        #diffs > diff_thresh will be 1, otherwise 0.
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
        #diffmaps in the region defined by 'detection'. Currently fits a
        #Gaussian in the y-direction for each x
        #use 'detection' to guess starting fit parameters.

        #get just the positive elements of the difference map. Perform fitting on
        #these positive diffmaps.
        posdiffs = copy.deepcopy(diffs)
        for i in range(0, len(diffs)):
            temp = diffs[i].data < 0
            posdiffs[i].data[temp] = 0

        #fit a function to the difference maps in the cases where there has been a
        #detection
        fitparams, wavefront = aware_utils.fit_wavefront(posdiffs, detection)

        #transform the detected model wavefront back into heliocentric coordinates so it can be overlayed
        wavefront_hc = aware_utils.map_reravel(wavefront,params,verbose=True)

        #strip out the velocity information from the wavefront fitting
        velocity = aware_utils.wavefront_velocity(fitparams)

        #strip out the position and width information from the wavefront fitting
        pos_width = aware_utils.wavefront_position_and_width(fitparams)

        #now save products we have created in a pickle file for future reference
        #Will save output in ~/aware_results
        extn=time_range.start().strftime('%Y%m%d_%H%M')
        save_path=os.path.expanduser('~/aware_results/')
        save_file='aware_results_' + extn + '.pickle'
        
        if not os.path.exists(save_path):
            os.makedirs(save_path)
                          
        output=open(save_path + save_file,'wb')
        print 'Saving result products to: '+ save_path + save_file
                           
        pickle.dump(wavefront_hc,output)
        output.close()

        #visualize the model wavefront
        visualize(wavefront_hc)
        
    return maps, new_maps, diffs, threshold_maps, binary_maps, detection_maps, wavefront, velocity, pos_width, persistence_maps, wavefront_hc


if __name__ == '__main__':
    main()

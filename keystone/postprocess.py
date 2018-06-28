import yaml
import io
import os
import keystone
from astropy.utils.data import get_pkg_data_filename

def get_baselineRegion(region='W48', dirname=None):
    # Read YAML file (contains clipping parameters by region)
    filename = get_pkg_data_filename('./data/postprocess.yaml',
                                     package='keystone')

    with open(filename, 'r') as stream:
    	data_loaded = yaml.load(stream)
    
    clip = data_loaded['clip_params']
    # Figure out which clipping parameters to grab based on dirname
    if region in clip:
        if ('NH3_11' in dirname):
	    line = clip[region]['nh3_11']
            baselineRegion = [slice(1362+line['lowFreq'], 1462+line['lowFreq'], 1), slice(2634-line['highFreq'], 2734-line['highFreq'], 1)]
	    startChannel = 1362+line['lowFreq']
	    endChannel = 2734-line['highFreq']
        elif ('NH3_22' in dirname):
	    line = clip[region]['nh3_22']
            baselineRegion = [slice(1362+line['lowFreq'], 1662+line['lowFreq'], 1), slice(2434-line['highFreq'], 2734-line['highFreq'], 1)]
	    startChannel = 1362+line['lowFreq']
            endChannel = 2734-line['highFreq']
        else:
	    line = clip[region]['nh3_22']
            baselineRegion = [slice(1500+line['lowFreq'], 1800+line['lowFreq'], 1), slice(2093-line['highFreq'], 2597-line['highFreq'], 1)]
	    startChannel = 1500+line['lowFreq']
            endChannel = 2597-line['highFreq']
	return baselineRegion, startChannel, endChannel

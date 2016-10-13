import subprocess
import glob
import os
import warnings
from astropy.table import Table, join
import numpy as np

def updateLogs(output='ObservationLog.csv',release=None):
    if (release is None) or ('all' in release):
        command = "wget --no-check-certificate --output-document="+output+" 'https://docs.google.com/spreadsheet/ccc?key=1319WA-stpit5_Hxe9MbY673ELhxC4p6oKG2l9eyEUM8&output=csv'"
        # The returns from subprocess are the error codes from the OS
        # If 0 then it worked so we should return True
        return not subprocess.call(command,shell=True)
    if 'DR1' in release:
        from astropy.utils.data import get_pkg_data_filename
        filename = get_pkg_data_filename('data/ObservationLog_DR1.csv',
                                         package='KEYSTONE')
        command = 'cp '+filename+' ./ObservationLog.csv'
        return not subprocess.call(command,shell=True)
    else:
        warnings.warn('Updating logs failed for non-existent data release.')
        return False

def updateCatalog(output='RegionCatalog.csv',release=None):
    if (release is None) or ('all' in release):
        command = "wget --no-check-certificate --output-document="+output+" 'https://docs.google.com/spreadsheet/ccc?key=12Vs1VeDfq3S7LmgvBs5_bbNvo4ui7_eiFpsGfGehlOQ&output=csv'"
        return not subprocess.call(command,shell=True)
    if 'DR1' in release:
        from astropy.utils.data import get_pkg_data_filename
        filename = get_pkg_data_filename('data/RegionCatalog_DR1.csv',
                                         package='KEYSTONE')
        command = 'cp '+filename+' ./'+output
        return not subprocess.call(command,shell=True)
    else:
        warnings.warn('Updating logs failed for non-existent data release.')
        return False

def GenerateRegions(refresh=False,release='all'):

    if refresh:
        updateLogs(release=release)
        updateCatalog(release=release)
    if not os.access('ObservationLog.csv',os.R_OK):
        updateLogs(release = release)
    if not os.access('RegionCatalog.csv',os.R_OK):
        updateCatalog(release = release)
    obs = Table.read('ObservationLog.csv')
    cat = Table.read('RegionCatalog.csv')
    cat.rename_column('Region/Box name','BoxName')
    obs.rename_column('Source','BoxName')

# This takes out rows that are empty
# This needs to be done twice for some reason... 
    for thiscat in (obs,obs,cat,cat,cat):
        for idx, row in enumerate(thiscat):
            if not row['BoxName']:
                thiscat.remove_row(idx)
    cat.replace_column('VLSR',np.asarray(cat['VLSR'],dtype=np.float))
    joincat = join(obs,cat,keys='BoxName')
    groupcat = joincat.group_by('Region name')
    min_values = groupcat.groups.aggregate(np.min)
    max_values = groupcat.groups.aggregate(np.max)
    mean_values = groupcat.groups.aggregate(np.mean)
    vavg = 0.5*(min_values['VLSR'] + max_values['VLSR'])
    vrange = max_values['VLSR']- min_values['VLSR']
    mean_values['VAVG'] = vavg
    mean_values['VRANGE'] = vrange

    return(mean_values)


def parseLog(logfile='ObservationLog.csv'):
    """
    Ingests a CSV log file into an astropy table
    """
    try:
        from astropy.table import Table
    except:
        warnings.warn('KEYSTONE Pipeline requires astropy.  Try Anaconda!')
        return
    t = Table.read(logfile)
    return(t)

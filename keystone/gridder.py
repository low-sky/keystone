#from sdpy import makecube
import numpy as np
import glob
from astropy.io import fits
import astropy.wcs as wcs
import itertools
from scipy.special import j1
import pdb
import numpy.fft as fft
import astropy.utils.console as console
import astropy.units as u
import astropy.constants as con
import numpy.polynomial.legendre as legendre
import warnings
import os
import gbtpipe
from .utils import VlsrByCoord
from . import __version__
from .postprocess import *
import postprocess


def baselineSpectrum(spectrum, order=1, baselineIndex=()):
    x = np.linspace(-1, 1, len(spectrum))
    coeffs = legendre.legfit(x[baselineIndex], spectrum[baselineIndex], order)
    spectrum -= legendre.legval(x, coeffs)
    return(spectrum)


def freqShiftValue(freqIn, vshift, convention='RADIO'):
    cms = 299792458.
    if convention.upper() in 'OPTICAL':
        return freqIn / (1.0 + vshift / cms)
    if convention.upper() in 'TRUE':
        return freqIn * ((cms + vshift) / (cms - vshift))**0.5
    if convention.upper() in 'RADIO':
        return freqIn * (1.0 - vshift / cms)


def channelShift(x, ChanShift):
    # Shift a spectrum by a set number of channels.
    ftx = np.fft.fft(x)
    m = np.fft.fftfreq(len(x))
    phase = np.exp(2 * np.pi * m * 1j * ChanShift)
    x2 = np.real(np.fft.ifft(ftx * phase))
    return(x2)


def jincGrid(xpix, ypix, xdata, ydata, pixPerBeam=None):
    a = 1.55 / (3.0 / pixPerBeam)
    b = 2.52 / (3.0 / pixPerBeam)

    Rsup = 1.09 * pixPerBeam  # Support radius is ~1 FWHM (Leroy likes 1.09)
    dmin = 1e-4
    dx = (xdata - xpix)
    dy = (ydata - ypix)

    pia = np.pi / a
    b2 = 1. / (b**2)
    distance = np.sqrt(dx**2 + dy**2)

    ind  = (np.where(distance <= Rsup))
    d = distance[ind]
    wt = j1(d * pia) / \
        (d * pia) * \
        np.exp(-d**2 * b2)
#    wt[ind] = np.exp(-distance[ind]**2*b2)*\
#              np.sin(pia*distance[ind])/\
#              (pia*distance[ind])
    wt[(d < dmin)] = 0.5  # Peak of the jinc function is 0.5 not 1.0

    return(wt, ind)


def VframeInterpolator(scan):
    # Find cases where the scan number is 
    startidx = scan['PROCSEQN']!=np.roll(scan['PROCSEQN'],1)
    scanstarts = scan[startidx]
    indices = np.arange(len(scan))
    startindices = indices[startidx]
    scannum = scanstarts['PROCSEQN']
    vfs = scanstarts['VFRAME']

    odds = (scannum % 2) == 1
    evens = (scannum % 2) == 0
    
    coeff_odds,_,_,_ = np.linalg.lstsq(\
        np.c_[startindices[odds]*1.0,
              np.ones_like(startindices[odds])],
        vfs[odds])

    coeff_evens,_,_,_ = np.linalg.lstsq(\
        np.c_[startindices[evens]*1.0,
              np.ones_like(startindices[evens])],
        vfs[evens])

    vfit = np.zeros(len(scan))+np.nan

    for thisone, singlescan in enumerate(scan):
        startv = vfs[scannum == singlescan['PROCSEQN']]
        startt = startindices[scannum == singlescan['PROCSEQN']]
        if singlescan['PROCSEQN'] % 2 == 0:
            endv = coeff_odds[1] + coeff_odds[0] * (startt+94)
        if singlescan['PROCSEQN'] % 2 == 1:
            endv = coeff_evens[1] + coeff_evens[0] * (startt+94)

        endt = startt+94
        try:
            vfit[thisone] = (thisone - startt) * \
                (endv - startv)/94 + startv
        except:
            pass
    return(vfit)
    
def autoHeader(filelist, beamSize=0.0087, pixPerBeam=3.0):
    RAlist = []
    DEClist = []
    for thisfile in filelist:
        s = fits.getdata(thisfile)
        try:
            RAlist = RAlist + [s['CRVAL2']]
            DEClist = DEClist + [s['CRVAL3']]
        except:
            pdb.set_trace()

    longitude = np.array(list(itertools.chain(*RAlist)))
    latitude = np.array(list(itertools.chain(*DEClist)))
    longitude = longitude[longitude != 0]
    latitude = latitude[latitude != 0]
    minLon = np.nanmin(longitude)
    maxLon = np.nanmax(longitude)
    minLat = np.nanmin(latitude)
    maxLat = np.nanmax(latitude)

    naxis2 = np.ceil((maxLat - minLat) /
                     (beamSize / pixPerBeam) + 2 * pixPerBeam)
    crpix2 = naxis2 / 2
    cdelt2 = beamSize / pixPerBeam
    crval2 = (maxLat + minLat) / 2
    ctype2 = s[0]['CTYPE3'] + '--TAN'
    # Negative to go in the usual direction on sky:
    cdelt1 = -beamSize / pixPerBeam

    naxis1 = np.ceil((maxLon - minLon) /
                     (beamSize / pixPerBeam) *
                     np.cos(crval2 / 180 * np.pi) + 2 * pixPerBeam)
    crpix1 = naxis1 / 2
    crval1 = (minLon + maxLon) / 2
    ctype1 = s[0]['CTYPE2'] + '---TAN'
    outdict = {'CRVAL1': crval1, 'CRPIX1': crpix1,
               'CDELT1': cdelt1, 'NAXIS1': naxis1,
               'CTYPE1': ctype1, 'CRVAL2': crval2,
               'CRPIX2': crpix2, 'CDELT2': cdelt2,
               'NAXIS2': naxis2, 'CTYPE2': ctype2}

    return(outdict)


def addHeader_nonStd(hdr, beamSize, Data_Unit):
    if Data_Unit == 'Tmb':
        hdr['BUNIT'] = 'K'
    hdr['BMAJ'] = beamSize
    hdr['BMIN'] = beamSize
    hdr['BPA'] = 0.0
    hdr['TELESCOP'] = 'GBT'
    hdr['INSTRUME'] = 'KFPA'
    return(hdr)

def gridall(region='NGC7538', **kwargs):
    suffix = ['NH3_22', 'NH3_33', 'NH3_44', 'NH3_55',
              'C2S_2_1', 'CH3OH_10_9', 'CH3OH_12_11', 
              'H2O', 'HC5N_8_7', 'HC5N_9_8', 'HC7N_19_18',
              'HNCO_1_0']
    griddata(region = region, dirname = region + '_NH3_11',
             outdir = './images/', rebase=True,
             **kwargs)
    templatehdr = fits.getheader('./images/' + region + '_NH3_11_all.fits')
    for thisline in suffix:
        griddata(region = region, dirname = region + '_' + thisline,
                 outdir = './images/', rebase=True,
                 templateHeader=templatehdr,
                 **kwargs)

def getGainDict():
    gainDict = {('0', '0') : 0.979128705,
                ('1', '0') : 0.916095905,
                ('2', '0') : 0.875017,
                ('3', '0') : 0.784742095,
                ('4', '0') : 0.93435506,
                ('5', '0') : 0.742008405,
                ('6', '0') : 0.87569747,
                ('0', '1') : 0.94387634,
                ('1', '1') : 0.86831252,
                ('2', '1') : 0.87293043,
                ('3', '1') : 0.780018805,
                ('4', '1') : 0.8046946,
                ('5', '1') : 0.532584695,
                ('6', '1') : 0.97160044}
    return(gainDict)

def griddata(pixPerBeam=3.5,
             templateHeader=None,
             gridFunction=jincGrid,
             rootdir='/lustre/pipeline/scratch/KEYSTONE/',
             region='NGC7538',
             dirname='NGC7538_NH3_11',
             startChannel=None, endChannel=None,
             doBaseline=True,
             baselineRegion=None,
             blorder=1,
             Sessions=None,
             file_extension=None,
             rebase=True, 
             beamSize=None,
             OnlineDoppler=True,
             flagRMS=True,
             flagRipple=True,
             flagSpike=True,
             blankSpike=True,
             rmsThresh=1.5,
             plotTimeSeries=True,
             gainDict=None,
             filelist=[],
             outdir=None, **kwargs):

    # This uses a fixed set of gain dicts.  
    if gainDict is None:
        gainDict = getGainDict()
    if outdir is None:
        outdir = os.getcwd()
    
    blRegion,sChannel,eChannel = postprocess.get_baselineRegion(region=region, dirname=dirname)
    if baselineRegion is None:
        baselineRegion=blRegion
    if startChannel is None:
        startChannel=sChannel
    if endChannel is None:
        endChannel=eChannel

# If the filelist is not specified, go ahead and build it from groupings.
    if not filelist:
        if not Sessions:
            filelist = glob.glob(rootdir + '/' + region + '/' + dirname + '/*fits')
            if not file_extension:
                file_extension = '_all'
            history_message = 'Gridding of data using all sessions'
        else:
            filelist = []
            for scan_i in Sessions:
                    filelist.extend(glob.glob(rootdir + '/' + region +
                                              '/' + dirname + '/*_sess' +
                                              str(scan_i) + '.fits'))
            if isinstance(Sessions, list):
                if not file_extension:
                    file_extension = '_sess{0}-sess{1}'.format(Sessions[0],
                                                               Sessions[-1])
                if (Sessions[-1] + 1. - Sessions[0]) / len(Sessions) == 1.0:
                    history_message = 'Gridding of data using sessions' \
                        'between {0} and {1}'.format(Sessions[0], Sessions[-1])
                else:
                    history_message = 'Gridding of data using sessions: '
                    for scan_i in Sessions:
                        history_message += '{0}, '.format(scan_i)
            else:
                if not file_extension:
                    file_extension = '_sess{0}'.format(Sessions)
                history_message = 'Gridding of data using session'\
                    '{0}'.format(Sessions)

    if len(filelist) == 0:
        warnings.warn('There are no FITS files to process '
                      'in ' + rootdir + '/' + region + '/' + dirname)
        return
    # check that every file in the filelist is valid
    # If not then remove it and send warning message
    for file_i in filelist:
        try: 
           fits.open(file_i)
        except:
            warnings.warn('file {0} is corrupted'.format(file_i))
            filelist.remove(file_i)
            
    outdir = rootdir + '/images/'
    outname = dirname + file_extension
    gbtpipe.Gridding.griddata(filelist,
                              startChannel=startChannel,
                              endChannel=endChannel,
                              doBaseline=doBaseline,
                              baselineRegion=baselineRegion,
                              blorder=blorder, rebaseorder=3,
                              flagRMS=flagRMS, rmsThresh=rmsThresh,
                              outdir=outdir, outname=outname,
                              templateHeader=templateHeader,
                              VlsrByCoord=VlsrByCoord,
                              plotTimeSeries=plotTimeSeries,
                              blankSpike=blankSpike,
                              rebase=rebase, gainDict=gainDict,
                              flagSpike=flagSpike, **kwargs)
    # Convolve the beam size up by 10% in size
    #gbtpipe.Gridding.postConvolve(outdir+outname)




# # pull a test structure
#     hdulist = fits.open(filelist[0])
#     s = hdulist[1].data
# #    s = fits.getdata(filelist[0])
    
# # Constants block
#     sqrt2 = np.sqrt(2)
#     mad2rms = 1.4826
#     prefac = mad2rms / sqrt2
#     c = 299792458.
#     nu0 = s[0]['RESTFREQ']

#     Data_Unit = s[0]['TUNIT7']
#     if beamSize is None:
#         beamSize = 1.22 * (c / nu0 / 100.0) * 180 / np.pi  # in degrees
#     naxis3 = len(s[0]['DATA'][startChannel:endChannel])

# # Default behavior is to park the object velocity at
# # the center channel in the VRAD-LSR frame

#     crval3 = s[0]['RESTFREQ'] * (1 - s[0]['VELOCITY'] / c)
#     crpix3 = s[0]['CRPIX1'] - startChannel
#     ctype3 = s[0]['CTYPE1']
#     cdelt3 = s[0]['CDELT1']

#     w = wcs.WCS(naxis=3)

#     w.wcs.restfrq = nu0
#     w.wcs.radesys = s[0]['RADESYS']
#     w.wcs.equinox = s[0]['EQUINOX']
#     # We are forcing this conversion to make nice cubes.
#     w.wcs.specsys = 'LSRK'
#     w.wcs.ssysobs = 'TOPOCENT'

#     if templateHeader is None:
#         wcsdict = autoHeader(filelist, beamSize=beamSize,
#                              pixPerBeam=pixPerBeam)
#         w.wcs.crpix = [wcsdict['CRPIX1'], wcsdict['CRPIX2'], crpix3]
#         w.wcs.cdelt = np.array([wcsdict['CDELT1'], wcsdict['CDELT2'], cdelt3])
#         w.wcs.crval = [wcsdict['CRVAL1'], wcsdict['CRVAL2'], crval3]
#         w.wcs.ctype = [wcsdict['CTYPE1'], wcsdict['CTYPE2'], ctype3]
#         naxis2 = wcsdict['NAXIS2']
#         naxis1 = wcsdict['NAXIS1']
#     else:
#         w.wcs.crpix = [templateHeader['CRPIX1'],
#                        templateHeader['CRPIX2'], crpix3]
#         w.wcs.cdelt = np.array([templateHeader['CDELT1'],
#                                 templateHeader['CDELT2'], cdelt3])
#         w.wcs.crval = [templateHeader['CRVAL1'],
#                        templateHeader['CRVAL2'], crval3]
#         w.wcs.ctype = [templateHeader['CTYPE1'],
#                        templateHeader['CTYPE2'], ctype3]
#         naxis2 = templateHeader['NAXIS2']
#         naxis1 = templateHeader['NAXIS1']
#     outCube = np.zeros((int(naxis3), int(naxis2), int(naxis1)))
#     outWts = np.zeros((int(naxis2), int(naxis1)))

#     xmat, ymat = np.meshgrid(np.arange(naxis1), np.arange(naxis2),
#                              indexing='ij')
#     xmat = xmat.reshape(xmat.size)
#     ymat = ymat.reshape(ymat.size)
#     xmat = xmat.astype(np.int)
#     ymat = ymat.astype(np.int)

#     ctr = 0
#     for thisfile in filelist:
#         ctr += 1
#         s = fits.open(thisfile)
#         print("Now processing {0}".format(thisfile))
#         print("This is file {0} of {1}".format(ctr, len(filelist)))

#         nuindex = np.arange(len(s[1].data['DATA'][0]))

#         if not OnlineDoppler:
#             vframe = VframeInterpolator(s[1].data)
#         else:
#             vframe = s[1].data['VFRAME']

#         for idx, spectrum in enumerate(console.ProgressBar((s[1].data))):
#             # Generate Baseline regions
#             baselineIndex = np.concatenate([nuindex[ss]
#                                             for ss in baselineRegion])

#             specData = spectrum['DATA']
#             # baseline fit
#             if doBaseline & np.all(np.isfinite(specData)):
#                 specData = baselineSpectrum(specData, order=blorder,
#                                             baselineIndex=baselineIndex)

#             # This part takes the TOPOCENTRIC frequency that is at
#             # CRPIX1 (i.e., CRVAL1) and calculates the what frequency
#             # that would have in the LSRK frame with freqShiftValue.
#             # This then compares to the desired frequency CRVAL3.
                
#             DeltaNu = freqShiftValue(spectrum['CRVAL1'],
#                                      -vframe[idx]) - crval3
#             DeltaChan = DeltaNu / cdelt3
#             specData = channelShift(specData, -DeltaChan)
#             outslice = (specData)[startChannel:endChannel]
#             spectrum_wt = np.isfinite(outslice).astype(np.float)
#             outslice = np.nan_to_num(outslice)
#             xpoints, ypoints, zpoints = w.wcs_world2pix(spectrum['CRVAL2'],
#                                                         spectrum['CRVAL3'],
#                                                         spectrum['CRVAL1'], 0)
#             tsys = spectrum['TSYS']
#             if flagRMS:
#                 radiometer_rms = tsys / np.sqrt(np.abs(spectrum['CDELT1']) *
#                                                 spectrum['EXPOSURE'])
#                 scan_rms = prefac * np.median(np.abs(outslice[0:-2] -
#                                                         outslice[2:]))

#                 if scan_rms > 1.25 * radiometer_rms:
#                     tsys = 0 # Blank spectrum
#             if flagRipple:
#                 scan_rms = prefac * np.median(np.abs(outslice[0:-2] -
#                                                      outslice[2:]))
#                 ripple = prefac * sqrt2 * np.median(np.abs(outslice))

#                 if ripple > 2 * scan_rms:
#                     tsys = 0 # Blank spectrum
                
#             if (tsys > 10) and (xpoints > 0) and (xpoints < naxis1) \
#                     and (ypoints > 0) and (ypoints < naxis2):
#                 pixelWeight, Index = gridFunction(xmat, ymat,
#                                                   xpoints, ypoints,
#                                                   pixPerBeam=pixPerBeam)
#                 vector = np.outer(outslice * spectrum_wt,
#                                   pixelWeight / tsys**2)
#                 wts = pixelWeight / tsys**2
#                 outCube[:, ymat[Index], xmat[Index]] += vector
#                 outWts[ymat[Index], xmat[Index]] += wts
#         # Temporarily do a file write for every batch of scans.
#         outWtsTemp = np.copy(outWts)
#         outWtsTemp.shape = (1,) + outWtsTemp.shape
#         outCubeTemp = np.copy(outCube)
#         outCubeTemp /= outWtsTemp

#         hdr = fits.Header(w.to_header())
#         hdr = addHeader_nonStd(hdr, beamSize, Data_Unit)
#         #
#         hdu = fits.PrimaryHDU(outCubeTemp, header=hdr)
#         hdu.writeto(outdir + '/' + dirname + '.fits', clobber=True)

#     outWts.shape = (1,) + outWts.shape
#     outCube /= outWts

#     # Create basic fits header from WCS structure
#     hdr = fits.Header(w.to_header())
#     # Add non standard fits keyword
#     hdr = addHeader_nonStd(hdr, beamSize, Data_Unit)
#     # Adds history message
#     try:
#         hdr.add_history(history_message)
#     except UnboundLocalError:
#         pass
#     hdr.add_history('Using KEYSTONE pipeline version {0}'.format(__version__))
#     hdu = fits.PrimaryHDU(outCube, header=hdr)
#     hdu.writeto(outdir + '/' + dirname + '.fits', clobber=True)

#     w2 = w.dropaxis(2)
#     hdr2 = fits.Header(w2.to_header())
#     hdu2 = fits.PrimaryHDU(outWts, header=hdr2)
#     hdu2.writeto(outdir + '/' + dirname + '_wts.fits', clobber=True)

#     if rebase:

#         if 'NH3_11' in dirname:
#             baseline.rebaseline(outdir + '/' + dirname + '.fits',
#                                 windowFunction=baseline.ammoniaWindow,
#                                 line='oneone', **kwargs)

#         elif 'NH3_22' in dirname:
#             winfunc = baseline.ammoniaWindow
#             baseline.rebaseline(outdir + '/' + dirname + '.fits',
#                                 windowFunction=baseline.ammoniaWindow,
#                                 line='twotwo', **kwargs)

#         elif 'NH3_33' in dirname:
#             baseline.rebaseline(outdir + '/' + dirname + '.fits',
#                                 winfunc = baseline.ammoniaWindow,
#                                 line='threethree', **kwargs)
#         else:
#             baseline.rebaseline(outdir + '/' + dirname + '.fits',
#                                 windowFunction=baseline.tightWindow, 
#                                 **kwargs)

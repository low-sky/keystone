import numpy
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares as lsq
import numpy.polynomial.legendre as legendre
import pprocess
import time
import sys
import skimage
from gbtpipe.Baseline import mad1d, legendreLoss

def get_mask(spectra, mask_percent=0.4):
        """  
 Returns a mask of channels to be used for baseline fitting.
 Function calculates the standard deviation of a 31 pixel window
 centred on each pixel. A percentage of the pixels with the lowest 
 standard deviation for their window are then chosen for baseline fitting.    
   
 spectra = input spectra as numpy array
 mask_percent = percentage of pixels to select for baseline fitting    
        """ 
        spec_len = len(spectra)
        sample = int(spec_len*mask_percent)
        left = numpy.zeros(15)+numpy.std(spectra[0:31]) # For the first 15 entries, use first window
        right = numpy.zeros(15)+numpy.std(spectra[-31:]) # For the last 15 entries, use last window
        middle = numpy.std(skimage.util.view_as_windows(spectra, 31, 1),axis=1)
        stds = numpy.concatenate((left, middle, right))
        mask = numpy.arange(spec_len)[numpy.argsort(stds)[:sample]]
        #median_std = numpy.std(stds)*3
        #mask = numpy.where(numpy.array(stds)<median_std)[0]
        #mask = numpy.concatenate((mask, numpy.arange(-50, 50)))
        #if 1>0: #len(mask)<30
                #print 'yes'
                #mask = numpy.arange(len(spectra))[numpy.argsort(stds)[:500]]
                #mask=numpy.arange(len(spectra))[0::5]

        #plt.plot(range(len(spectra)), spectra)
        #plt.plot(numpy.arange(len(spectra))[mask], spectra[mask])
        #plt.show()
        return mask

def redchisqg(ydata,ymod,deg=2,sd=None):     
      """Returns the reduced chi-square error statistic for an arbitrary model,   
 chisq/nu, where nu is the number of degrees of freedom. If individual   
 standard deviations (array sd) are supplied, then the chi-square error   
 statistic is computed as the sum of squared errors divided by the standard   
 deviations. See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.  
   
 ydata,ymod,sd assumed to be Numpy arrays. deg integer.  
   
 Usage:  
 >>> chisq=redchisqg(ydata,ymod,n,sd)  
 where  
  ydata : data  
  ymod : model evaluated at the same x points as ydata  
  n : number of free parameters in the model  
  sd : uncertainties in ydata  
   
 Rodrigo Nemmen  
 http://goo.gl/8S1Oo"""  
      # Chi-square statistic  
      if sd==None:  
           chisq=numpy.sum((ydata-ymod)**2)  
      else:  
           chisq=numpy.sum( ((ydata-ymod)/sd)**2 )  
             
      # Number of degrees of freedom assuming 2 free parameters  
      nu=ydata.size-1-deg  
        
      return chisq/nu

def get_chi(blorder, ydata, xdata, blindex, noise):
    """
 Returns the best-fit Legendre polynomial to an input spectra, 
 along with that model's reduced chi-squared value.     
   
 blorder = order of polynomial to fit
 ydata = spectrum y values
 xdata = spectrum x values
 blindex = the indices of the spectra to include in baseline fitting
 noise = rms noise of spectrum"""
    opts = lsq(legendreLoss, np.zeros(blorder + 1), args=(ydata[blindex],
                                                          xdata[blindex],
                                                        noise), loss='arctan')
    ymod = legendre.legval(xdata, opts.x)
    chi = redchisqg(ydata,ymod,deg=blorder+1)
    return ymod, chi

def robustBaseline_chi(y, baselineIndex, blorder_max=3, noiserms=None):
    """  
 Returns a baseline subtracted spectrum, based on the best-fitting polynomial 
 for a range of polynomial orders.    
   
 y = input spectra
 baselineIndex = indices of spectra to include in baseline fitting
 blorder_max = largest order polynomial to fit (fit from blorder_max down to order of 1) 
 noiserms = rms noise of spectrum"""
    x = np.linspace(-1, 1, len(y))
    if noiserms is None:
        noiserms = mad1d((y - np.roll(y, -2))[baselineIndex])
    out = []
    for i in numpy.arange(1,blorder_max+1):
        b = get_chi(blorder=i, ydata=y, xdata=x, blindex=baselineIndex, noise=noiserms)
        out.append(b)
    out = numpy.array(out)
    find_low = numpy.where(out[:,1]==min(out[:,1]))
    low_model = out[find_low][0][0]
    #print low_model
    #plt.plot(range(len(y)), y)
    #plt.plot(numpy.arange(len(y))[baselineIndex], y[baselineIndex])  
    #plt.plot(range(len(x)), low_model, color='red')
    #plt.show()
    return y - low_model

def rebase(i, j, spectra, mask_percent=0.4, blorder_max=3):
        """  
 Parallelizable function to feed into pprocess. Returns a baseline-subtracted
 spectrum and its indices on the image plane.    
   
 i,j = x and y indices for spectrum location on image plane
 spectra = spectrum on which a baseline will be subtracted
 mask_percent = percentage of pixels to select for baseline fitting
 blorder_max = largest order polynomial to fit (fit from blorder_max down to order of 1) 
        """
        if (False in numpy.isnan(spectra)): #and (m/std > 10.):
                mask = get_mask(spectra, mask_percent=mask_percent)
                spectra = robustBaseline_chi(spectra, mask, blorder_max=blorder_max, noiserms=None)
        return i, j, spectra

def rebase_multi(filename, nproc=8, mask_percent=0.4, blorder_max=3):
        """  
 Returns a baseline-subtracted cube. Can be run with parallel processes.    
   
 filename = name of datacube to process (including its path)
 nproc = number of parallel processes desired
 mask_percent = percentage of pixels to select for baseline fitting
 blorder_max = largest order polynomial to fit (fit from blorder_max down to order of 1) 
        """
        cube = SpectralCube.read(filename)
        data = numpy.array(cube.unmasked_data[:,:,:])   

        queue = pprocess.Queue(limit=nproc)
        calc = queue.manage(pprocess.MakeParallel(rebase))
        tic = time.time()

        # create cube to store rebaselined data
        #cube_out = data.copy()
        cube_out = np.zeros(cube.shape) * np.nan
        shape = numpy.shape(data)
        pixels = shape[1] * shape[2]

        counter = 0
        for (i,j), value in numpy.ndenumerate(data[0]):
                calc(i,j, spectra=numpy.array(data[:,i,j]), mask_percent=mask_percent, blorder_max=blorder_max)

        for i, j, ss in queue:
                cube_out[:,i,j]=ss
                counter+=1
                print(str(counter) + ' of ' + str(pixels) + ' pixels completed \r'),
                sys.stdout.flush()
        print("\n %f s for parallel computation." % (time.time() - tic))
        
        cube_final = SpectralCube(data=cube_out, wcs=cube.wcs, header=cube.header)
        cube_final.write(filename[0:-5] + '_rebase_multi.fits', format='fits', overwrite=True)

# Usage examples:
#rebase_multi('/lustre/pipeline/scratch/KEYSTONE/images/W48_NH3_11_sess31-sess31.fits')
#rebase_multi('/lustre/pipeline/scratch/KEYSTONE/images/NGC2264_NH3_11_all.fits', mask_percent=0.2)
#rebase_multi('/lustre/pipeline/scratch/KEYSTONE/images/Rosette_NH3_11_all.fits', mask_percent=0.4)
#rebase_multi('/lustre/pipeline/scratch/KEYSTONE/images/NGC2264_NH3_22_all.fits', mask_percent=0.2, blorder_max=1)

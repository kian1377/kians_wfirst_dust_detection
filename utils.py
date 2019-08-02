"""
convenience utilities for image processing from

https://github.com/douglase/poppy_nulling/

"""
from numpy.lib.stride_tricks import as_strided as ast
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

def downsample(A, block= (2,2), subarr=False):
    '''
    downsample, a downsampling function
    Conserves flux.
    
    Take a 2D numpy array, A, break it into subarrays of shape bloc
    k and sum the counts in them, returning a new array of the summed blocks.
    if subarray=True, if the size of the last block spills out of the array the block 
    will be discarded, otherwise an error will be returned.
    
    striding approach from: http://stackoverflow.com/a/5078155/2142498
    
    as_strided() is not limited to the memory block of your array, so added  check of dimensions.
    http://scipy-lectures.github.io/advanced/advanced_numpy/index.html#stride-manipulation-label
    '''
    if (np.remainder(block[0],np.floor(block[0])) !=0) or (np.remainder(block[0],np.floor(block[0])) !=0) :
        raise ValueError("Block size must be integers")

    if (np.remainder(A.shape[0],block[0]) !=0):
        if subarr:
            A=A[0:block[0]*int(np.floor(A.shape[0]/block[0])),:]
        else:
            raise ValueError("not an integer number of blocks in first dimension")
            
    if (np.remainder(A.shape[1],block[1]) !=0):
        if subarr:
            A=A[:,0:block[1]*int(np.floor(A.shape[1]/block[1]))]
        else:
            raise ValueError("not an integer number of blocks in second dimension")
    #shape of new array:
        
    shape= (A.shape[0]/ block[0], A.shape[1]/ block[1])+ block

    #strides that fill the new array:
    strides= (block[0]* A.strides[0], block[1]* A.strides[1])+ A.strides
    #create an array of sub_arrays
    
    blocked= ast(A, shape= shape, strides= strides)
    #sum along the third and fourth axes, forming the rebinned, flux conserved array  
    #print(blocked)

    binned=blocked.sum(axis=(2,3))
    return binned

def add_poisson_noise(photons):
	'''takes a numpy array of  values and finds a 
	random number from a poisson distribution centered 
	on the number of photons in that bin'''
	from scipy.stats import poisson
	
	vpoisson_rvs=np.vectorize(poisson.rvs) 

	if str(type(input)) == "<class 'astropy.io.fits.hdu.hdulist.HDUList'>":
		noisy_array=vpoisson_rvs(photons[0].data)
		noisy = fits.HDUList(fits.PrimaryHDU(data=noisy_array,header=photons[0].header))
	else:
		noisy=vpoisson_rvs(photons)
	return noisy


def add_poisson_noise(photons):
	'''Takes a numpy array of values and finds a 
	random number from a poisson distribution centered 
	on the number of photons in that bin
    
    '''
	from scipy.stats import poisson
	
	vpoisson_rvs=np.vectorize(poisson.rvs) 

	if str(type(input)) == "<class 'astropy.io.fits.hdu.hdulist.HDUList'>":
		noisy_array=vpoisson_rvs(photons[0].data)
		noisy = fits.HDUList(fits.PrimaryHDU(data=noisy_array,header=photons[0].header))
	else:
		noisy=vpoisson_rvs(photons)
	return noisy


from skimage.transform import resize
import scipy.ndimage
import astropy.units as u

def apply_det(image,
              ps_det=0.0211*u.arcsec/u.pixel,#CGI_Imaging_Pixel_Scale
              ps_input=0.005*u.arcsec/u.pixel,#OS6 PSF simulations
              det_shape=(150,150)):
    try:
        unit=image.unit
    except Exception as err:
        print(err)
        print("setting units to 1")
        unit=1
    #zoomed= resize(image/image.sum(),shape,anti_aliasing=True,order=1) #(image, output_shape, order=1, mode='reflect', cval=0, clip=True, preserve_range=False, anti_aliasing=True, anti_aliasing_sigma=None)
    #return 
    zoom = (ps_input/ps_det).decompose().value
    resampled_image = scipy.ndimage.interpolation.zoom(image, zoom,
                                                                 output=image.dtype,
                                                                 order=3)
    
    lx, ly = resampled_image.shape
    print("zoom: "+str(zoom))
    # crop down to match size of detector:
    if det_shape == None:
        lx_w,ly_w = resampled_image.shape
    else:
        lx_w, ly_w = det_shape
    border_x = np.abs(lx - lx_w) // 2
    border_y = np.abs(ly - ly_w) // 2


    if zoom<1:
        new_im = np.zeros([lx_w, ly_w])
        new_im[border_x:border_x + resampled_image.shape[0],
                                        border_y:border_y + resampled_image.shape[1]] = resampled_image
    else:
        new_im=resampled_image[border_x:border_x + lx_w, border_y:border_y + ly_w]
    # print(image.sum())
    # print(new_im.sum())
    renormed_im = new_im*image.sum()/new_im.sum()
    #print(renormed_im.sum())
    return renormed_im #check this conserves flux well enough
    
def displ_scale(array,ps=1*u.arcsec,ax=None,cmap=plt.cm.viridis,**kwargs):
    nx = np.int(array.shape[1])
    ny= np.int(array.shape[0])
    halfpix = ps.to(u.arcsec).value*0.5
    extent = ([-(nx/2*ps).to(u.arcsec).value-halfpix,
                (nx/2*ps).to(u.arcsec).value+halfpix,
                -(ny/2*ps).to(u.arcsec).value-halfpix,
                +(ny/2*ps).to(u.arcsec).value+halfpix
              ])
    #plt.figure(figsize=[9,4])
    #plt.subplot(121)
    print(nx,ny,ps,halfpix)

    try:
        plt.imshow(array.decompose().value,interpolation='nearest',extent=extent,cmap=cmap,**kwargs)
        cm=plt.colorbar(format='%.0e')
        cm.set_label(array.decompose().unit)
    except:
        plt.imshow(array,interpolation='nearest',extent=extent,**kwargs)
        plt.colorbar(format='%.0e')
        
    #plt.grid()

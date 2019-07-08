
import astropy.units as u
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

try:
    import cv2
    USE_OPENCV=True
except:
    USE_OPENCV=False
    import scipy.ndimage

def rotate(img,angle):
        n_y,n_x=img.shape
        if USE_OPENCV:
                #print("using openCV for rotations" )
                return cv2.warpAffine(img,
                                  cv2.getRotationMatrix2D((n_y/2,n_x/2),angle,1),
                                  (n_y,n_x),
                                  flags=cv2.INTER_AREA)
        else:
                return scipy.ndimage.rotate(img,angle,reshape=False,order=3)

def fit_to_hlc(fitsfile,
               interpolating_function,
               xmas,
               HLC_plate_scale,
               n = 200,#zodi_fits[0].data.shape[0]
               thresh=10,
               display=False,
               core_mask_radius=0*u.arcsec,
               mask_radius=None):

    zodi_fits = fits.open(fitsfile)
    index = zodi_fits[0].data < zodi_fits[0].data.max()/thresh
    zodi = np.ma.masked_array(zodi_fits[0].data,index)
    zodi_pixscale = zodi_fits[0].header["PIXELSCL"]*u.arcsecond
    pixnum = np.int(zodi_fits[0].data.shape[0])
    xpix,ypix = np.meshgrid(np.arange(-pixnum/2,pixnum/2),np.arange(-pixnum/2,pixnum/2))
    x = (xpix+.5).flatten()*zodi_pixscale
    y = (ypix+.5).flatten()*zodi_pixscale
    xcenter = (pixnum/2+.5)*zodi_pixscale
    ycenter = (pixnum/2+.5)*zodi_pixscale
    index[(np.sqrt((x)**2 + (y)**2)<core_mask_radius).reshape([pixnum,pixnum])]=True
    
    if display:
        plt.imshow(index,cmap=plt.cm.bone)
        plt.title("Mask")
        #plt.colorbar()
        
    xpix_flat = xpix.flatten()

    ypix_flat = ypix.flatten()

    im = np.zeros([n,n])
    
    if display:
        plt.figure()
        plt.imshow(zodi)
    for i,zodi_val in enumerate(zodi.flatten()):

        if np.ma.is_masked(zodi_val):
            continue
        im += zodi_val*closest_monochrome_PSF(x[i],y[i],
                                              interpolating_function,
                                              xmas,
                                              HLC_plate_scale,
                                              mask_radius=mask_radius)


    return im.T,zodi
   

def closest_monochrome_PSF(x,y,
                           HLCinterpfun,
                           xmas,
                           HLC_plate_scale,
                           n=200,# size of arrays from Krist
                           display=False,
                           mask_radius=None):
    
    '''
    Returns the nearest PSF realization by interpolating across grid of input 
    
    '''
    
    r,theta = tilt_converter(x,y)
    nearest=(np.abs(xmas - r)).argmin()
    grid=np.meshgrid(range(n),range(n))
    mask=1
    if mask_radius is not None:
        if mask_radius > n:
            raise ValueError("Mask radius exceeds array size")
        xp,yp=((+n/2-x/HLC_plate_scale).decompose(),(y/HLC_plate_scale).decompose()+n/2)
        mask=np.sqrt((-grid[0]+yp)**2+(grid[1]-xp)**2)<=mask_radius
        #plt.imshow(mask)
        #plt.show()
        #plt.figure()
    
    r=r.to(u.milliarcsecond)
    
    pts = np.vstack([np.vstack([grid[0].flatten(),grid[1].flatten()]),
                   r.value*np.ones(len(grid[0].flatten()))]).T
    
    interpped = rotate(HLCinterpfun(pts).reshape(n,n).T,
                     -theta.to(u.deg).value)*mask
    #print(theta)
    return interpped


## Calculate the theta/r value -- This is annoying since the wavefront is in-xy space.
## but the wfirst CGI doesn't inherit OpticalSystem, so I don't see an easy way to redefine 
## input_wavefront() to use xy coordinates

# These functions should move to an external .py module once done debugging.

def tilt_converter(x,y):
    """
    converts from xy coordinates to angle and rotation angle

    Parameters
    ----------  
    x: float
     pixel position
    y: float
     pixel position

    
    Returns
    ----------  

    (r,theta):
            where r is in the units of x and y and theta is in degrees
    
    """
    return (np.sqrt(x**2+y**2), np.arctan2(-x,y))*u.deg

def cartesian_off_axis_psf_poppy(cgi,x,y,wavels,**kwargs):
    r,theta = tilt_converter(x,y)
    #print(x,y,r,theta)
    ifs_spc.options['source_offset_r'] = r # arcsec
    ifs_spc.options['source_offset_theta'] = theta # deg w.r.t. North
    ifs_psf = ifs_spc.calc_datacube(wavels, **kwargs)
    return ifs_psf

import scipy.ndimage
import scipy.interpolate
def cartesian_off_axis_psf_interpol(cgi,x,y,wavels,**kwargs):
    r,theta = tilt_converter(x,y)
    PSF = fits.open()
    return ifs_psf

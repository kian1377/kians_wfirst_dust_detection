from importlib import reload
import hlc_processing
reload(hlc_processing)

import astropy.units as u
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

import os
my_home_dir = os.getcwd()

def disk_through_hlc(fitsfile,
                     my_interp_fun,
                     xmas,
                     HLC_plate_scale_AS,
                     thresh=1000,
                     display=False,
                     x_extent=.5,y_extent=0.5,
                     load_existing=False,
                     localzodi=0):
    
    if fitsfile.find("89") !=-1:
        thresh=50000
    
    print("Running disk through HLC ...")
    zodi_fits = fits.open(fitsfile)
    fitsfile_parts = fitsfile.split("\\")
    outfits = my_home_dir+"\KianDebesModels_OS5\\"+fitsfile_parts[-2]+"\\" + fitsfile_parts[-1][:-5]+"_HLC_"+str(thresh)+'.fits'
    print("\nInput file: ")
    print(fitsfile)
    print("\nOutput file: ")
    print(outfits)

    print("\nThresh = " + str(thresh))
    print("local zodi: {}".format(localzodi))

    zodi = np.ma.masked_array(zodi_fits[0].data, zodi_fits[0].data < zodi_fits[0].data.max()/thresh)
    zodi_fits[0].data = zodi
    
    n = 128
    zodi_pixscale = zodi_fits[0].header["PIXELSCL"]*u.arcsecond

    pixnum = np.int(zodi_fits[0].data.shape[0])
    x,y = np.meshgrid(np.arange(-pixnum/2,pixnum/2),np.arange(-pixnum/2,pixnum/2))
    x = (x+.5).flatten()*zodi_pixscale

    y = (y+.5).flatten()*zodi_pixscale
    im = np.zeros([n,n])

    zodi.mask#[zodi<zodi.max()/100] = True
    
    if display:
        plt.imshow(zodi)
    
    try:
        if load_existing:
            print("\nload_existing keyword True.")
            im = fits.getdata(outfits)
        else:
            print("\nload_existing keyword False, throwing error to regenerate file.")
            raise ValueError("")

    except:
        for i,zodi_val in enumerate(zodi.flatten()):
            if np.ma.is_masked(zodi_val):
                continue
            im += (localzodi+zodi_val)*hlc_processing.closest_monochrome_PSF(x[i],y[i],
                                                                             my_interp_fun,
                                                                             xmas,
                                                                             HLC_plate_scale_AS,
                                                                             n=128)
    if display:
        halfpix = zodi_pixscale.to(u.arcsec).value*0.5
        extent = ([x.min().to(u.arcsec).value-halfpix,
                   x.max().to(u.arcsec).value+halfpix, 
                   y.min().to(u.arcsec).value-halfpix,
                   y.max().to(u.arcsec).value+halfpix])
        
        plt.figure(figsize=[9,4])
        
        plt.subplot(121)
        plt.title("Input Flux")
        plt.xlim([-x_extent,x_extent])
        plt.ylim([-y_extent,y_extent])
        plt.ylabel("$\prime\prime$")
        plt.xlabel("$\prime\prime$")
        plt.grid()
        plt.imshow(zodi,extent=extent,norm=LogNorm(zodi.data[zodi>0].min(),zodi.data.max()))#zodi.data.min(),zodi.data.max))
        axc = plt.colorbar()
        axc.set_label("Jy")
        
        plt.subplot(122)
        pixnum = np.int(im.shape[0])
        x,y = np.meshgrid(np.arange(-pixnum/2,pixnum/2),np.arange(-pixnum/2,pixnum/2))
        x = (x+.5).flatten()*HLC_plate_scale_AS
        y = (y+.5).flatten()*HLC_plate_scale_AS
        halfpix = .05*HLC_plate_scale_AS.value
        extent = ([x.min().to(u.arcsec).value-halfpix,
        x.max().to(u.arcsec).value+halfpix, 
        y.min().to(u.arcsec).value-halfpix,
        y.max().to(u.arcsec).value+halfpix])
        plt.imshow(im.T,extent=extent)
        plt.colorbar()
        plt.title("HLC Image of Included Flux")
        plt.xlim([-x_extent,x_extent])
        plt.ylim([-y_extent,y_extent])
        plt.grid()
        plt.xlabel("$\prime\prime$")
        plt.tight_layout()
        
        plt.show()
        
    return im,zodi_fits,outfits

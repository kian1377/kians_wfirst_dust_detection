WFIRST-CGI Off-axis PSFs in HLC narrow FOV imaging mode
----------------------------------------------------------------------

Hanying Zhou
Jet Propulsion Laboratory, California Institute of Technology
January 7, 2019
Cleared for public distribution.



Realistic off-axis PSF images in HLC narrow FOV imaging mode (filter Band 1) are provided for science modeling purposes. 
 
The simulation is based on a recent HLC design for the (pre-) phase B pupil, and uses near-latest (phase A) OTA misalignment + OTA+CGI surface aberrations + (phase B) polarization effects.  



----------------------------------------------------------------------
This file package contains five *.fits files, one png file, and one .txt file:

20180718_hlc_nfov_PSFs_1Dcube.fits
20180718_hlc_nfov_PSFs_1Dcube_info.fits
20180718_hlc_nfov_PSFs_2Dcube.fits
20180718_hlc_nfov_PSFs_2Dcube_info.fits
DH_mask0.fits
offaxis_psfs_2D_grid.png
readme.txt

----------------------------------------------------------------------
Description of PSF generation:

The model consists of a full Fresnel diffraction model for high-accuracy contrast truth and an economical compact model for realistic wavefront control (dark hole digging), mimicking flight-constrained EFC performance.  

PSFs are taken post dark hole digging of a system*  with CBE (Phase A) OTA misalignments + OTA/CGI surface aberrations + (Phase B) polarization effects. The dark hole is static, and no jitter was added. No detector noise was added.

Post EFC control, off-axis sources are introduced as wavefront tilts, one at a time, and each resulting PSF image is obtained. A coarse 2D off-axis grid and a finer 1D off-axis scan of PSFs are simulated.


HLC design: Dwight Moody 20180925-251 (narrow fov imager 1, 575nm @ 10% bw), for 20180718 (pre-) Phase B pupil. It includes an occulter mask and a Lyot stop. No field stop is used in this simulation, to provide full field information for users.

The dark hole region is from 2.5~9 lambda/D, for D = 2.3631 m.

Detector plane sampling: 0.2 lam/D per pixel 



* Most details about the model, the system condition, and the modeling tool, are described in:
1) Hanying Zhou, et al., "High Accuracy Coronagraph Flight Model for 
WFIRST-CGI Raw Contrast Sensitivity Analysis", Proc. SPIE, 10698-92 (2018)

2) John Krist, et al.,“WFIRST coronagraph optical modeling,” Proc. SPIE,10400-4 (2017)

3) John Krist, "PROPER: an optical modeling program for IDL," Proc. SPIE,66750P (2007) 
Also:  http://proper-library.sourceforge.net


----------------------------------------------------------------------
Details on the HLC coarse 2D and fine 1D PSF cubes and their off-axis info.

Column variables in xx_info.fits files:

1. ith PSF (as in xx_cube.fits file)
2. source x-offset in lam/D
3. source y-offset in lam/D
4. core throughput, in %
   Fraction of the power inside the PSF core (those >=0.5*PSF peak) 
   relative to the power incident on the telescope primary	
5. core area, in arcsec^2
   Pixels whose value exceed half max of PSF peak	


Data in xx_2Dcube.fits file:
2D PSF of 128x128 size, total 429 
each offset is sequentially defined in the companion xx_info.fits
resolutions:
0.5 lam/D for 2.5- 5lam/D
1   lam/D for the rest of 0-9 lam/D
see also the grid.png file

Data in xx_1Dcube.fits file:
2D PSF of 128x128 size, total 85 
each offset is sequentially defined in the companion xx_info.fits 
resolutions:
3mas (or 0.06 lam/D) for <4 lam/D 
0.5 lam/D for 4-9 lam/D
0.25lam/D for 9-11 lam/D

DH_mask0.fits
This file can be used as visualization aid (use 'contour' in Matlab)
to superimpose the mask on to a PSF image to see the boundary of dark hole area 
or to mask out (multiply directly with a PSF image) light outside the dark hole region

----------------------------------------------------------------------

contact hanying.zhou@jpl.caltech.edu for any questions or suggestions


(c) 2019 California Institute of Technology. Government sponsorship acknowledged.
pro run_zodipic
  !PATH=!PATH+':'+Expand_Path('+~/github/zodipic')
  !PATH=!PATH+':'+Expand_Path('+~/github/IDLAstro')
    !PATH=!PATH+':'+Expand_Path('+~/github/cmsvlib')

  pixel_size=3              ;mas 
  wavel=0.5750                   ;um
  ;radin = .1
  distance=10
  alpha=1.5
  zodipic, fnu, pixel_size, wavel, ring=0, bands=0, blob=0, pixnum=256, distance=distance, alpha=alpha 
  
  mkhdr, h, fnu
  
  sxaddpar, h, 'FLUXUNIT', 'Jy', 'Flux units of each pixel'
  sxaddpar, h, 'PIXELSCL', pixel_size/1000, 'pixel dimension, arcsec'
  sxaddpar, h, 'dist', distance, 'distance to source, parsecs'
  sxaddpar, h, 'wavel', wavel, 'wavelength, microns'
  sxaddpar, h, 'alpha',alpha, 'radial power law of dust distribution'

  writefits, "zodipic_"+STRING(distance, FORMAT='(I2.2)')+"pc"+STRING(pixel_size, FORMAT='(I2.2)')+"mas.fits",$
             fnu, h, /CHECKSUM ;;Syntax - WRITEFITS, filename, data,[ header, /APPEND,]

                                ; Run higheralpha
  alpha=2.0
  sxaddpar, h, 'alpha', alpha, 'radial power law of dust distribution'
  zodipic, fnu, pixel_size, wavel, ring=0, bands=0, blob=0, pixnum=256, distance=distance, alpha=alpha 
 
  writefits, "zodipic_"+STRING(distance, FORMAT='(I2.2)')+"pc"+STRING(pixel_size, FORMAT='(I2.2)')+"mas"+STRING(alpha, FORMAT='(I2.2)')+"alpha.fits",$
             fnu, h, /CHECKSUM ;;Syntax - WRITEFITS, filename, data,[ header, /APPEND,]
  
   alpha=1.0
  sxaddpar, h, 'alpha', alpha, 'radial power law of dust distribution'
  zodipic, fnu, pixel_size, wavel, ring=0, bands=0, blob=0, pixnum=256, distance=distance, alpha=alpha 
 
  writefits, "zodipic_"+STRING(distance, FORMAT='(I2.2)')+"pc"+STRING(pixel_size, FORMAT='(I2.2)')+"mas"+STRING(alpha, FORMAT='(I2.2)')+"alpha.fits",$
             fnu, h, /CHECKSUM ;;Syntax - WRITEFITS, filename, data,[ header, /APPEND,]
stop
                                ;clear inside IWA
  radin=0.15*distance
  zodipic, fnu, pixel_size, wavel, ring=0,bands=0, blob=0, pixnum=256, radin=radin,distance=distance             ;,;starname='Epsilon_Eridani', zodis=1,radin=radin

  mkhdr, h, fnu
  
  sxaddpar, h, 'FLUXUNIT', 'Jy', 'Flux units of each pixel'
  sxaddpar, h, 'PIXELSCL', pixel_size/1000, 'pixel dimension, arcsec'
  sxaddpar, h, 'dist', distance, 'distance to source, parsecs'
  sxaddpar, h, 'wavel', wavel, 'wavelength, microns'
  sxaddpar, h, 'radin', radin, 'inner radius cutoff, AU'

  writefits, "zodipic_"+STRING(distance, FORMAT='(I2.2)')+"pc"+STRING(pixel_size, FORMAT='(I2.2)')+"mas_cutoff.fits", fnu, h, /CHECKSUM ;;Syntax - WRITEFITS, filename, data,[ header, /APPEND,]

stop

;eps eri
      dist=3.2                      ;pc, for converting arcsec to AU
  radin = 0.75*1.25/2
  radout=1.25*dist/5
  zodis=723; LBTI
  inclination=17.3              ;Macgregor 2017 outer belt gaussian fit
  magAO_pix=8
    zodipic, fnu, magAO_pix, 1,inclination=inclination,radin=radin,radout=radout,zodis=zodis, pixnum=1024,starname='Epsilon_Eridani';,/ADDSTAR;,dustsize=3
  mkhdr, h, fnu

  sxaddpar, h, 'FLUXUNIT', 'Jy', 'Flux units of each pixel'
  sxaddpar, h, 'PIXELSCL', magAO_pix, 'pixel dimension, mas'


    writefits, "epseri.fits", fnu, h, /CHECKSUM ;;Syntax - WRITEFITS, filename, data,[ header, /APPEND,]

  ;print,fnu
end

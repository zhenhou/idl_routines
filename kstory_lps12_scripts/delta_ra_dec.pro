;;;
; NAME: delta_ra_dec.pro
; PURPOSE:
;   General script for today
;
; NOTES:
;
; MODIFICATION HISTORY:
;  07/20/2012: (KTS) Created
;;;

;;; Get field properties for table 1
pro delta_ra_dec
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
adir = '/home/kstory/lps12/masks/apod/'

dra    = fltarr(20)
ddec   = fltarr(20)
thresh = 0.5
area   = fltarr(20)
aarea  = fltarr(20) ; apod mask only


for i=0, 19 do begin
    if i eq 2 then i=3

    mask = get_lps12_mask(i)

    nx = n_elements(mask[*,0])
    ny = n_elements(mask[0,*])
    radec0 = [f[i].ra0, f[i].dec0]
    radec0[0] = 180.
    xc = round(0.5*nx)
    yc = round(0.5*ny)
    xind = findgen(nx)
    yind = findgen(ny)

; dRA
    tmp = mask[*,yc]
    wh=where(xind lt 0.8*xc)
    xmin = wh[where_closest(tmp[wh],thresh*max(tmp))]
    wh=where(xind gt 1.2*xc)
    xmax = wh[where_closest(tmp[wh],thresh*max(tmp))]

    pix2ang_proj5, [nx,ny], radec0, 1., ra, dec, $
                   xpix=[xmin,xmax], ypix=[yc,yc]

    dra[i] = abs(ra[1]-ra[0])

; dDEC
    tmp = mask[xc,*]
    wh=where(yind lt 0.8*yc)
    ymin = wh[where_closest(tmp[wh],thresh*max(tmp))]
    wh=where(yind gt 1.2*yc)
    ymax = wh[where_closest(tmp[wh],thresh*max(tmp))]

    pix2ang_proj5, [nx,ny], radec0, 1., ra, dec, $
                   xpix=[xc,xc], ypix=[ymin,ymax]

    ddec[i] = abs(dec[1]-dec[0])

; area
    restore, adir+'apod_'+f[i].name+'_60_0.0500_30.sav'
    area[i] = total(mask) / (60^2.)
    aarea[i] = total(apod) / (60^2.) ; apod only
     ;print,f[i].name
    print, sigfig(f[i].ra0,3),' & ',sigfig(f[i].dec0,3),' & ',$
      sigfig(dra[i],3),' & ',sigfig(ddec[i],3), ' & ',$
      sigfig(area[i],4)
      

endfor

print, "total area = ", total(area)
print, "are masked by point sources = ", total(aarea) - total(area)
stop
END

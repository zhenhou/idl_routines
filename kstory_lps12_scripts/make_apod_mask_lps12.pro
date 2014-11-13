;;;
; NAME: make_apod_mask_lps12
; PURPOSE:
;   Explore how to make apodization masks
;
; INPUTS:
;   overlap,      4 element array in deg, [+x, -x, +y, -y], see get_map_overlap
;
; NOTES:
; 1) Now using the lpf-ds IDF's from Ryan, Jan3 2012
; 2) Need to check hanning choices
;
; MODIFICATION HISTORY:
;  02/10/2012: (KTS) Created from /home/cr/code/spt/coverage/mk_mask_run4_2009.pro
;  02/21/2012: (KTS) Add info from /home/rkeisler/ps09/make_apod.pro
;  03/06/2012: (KTS) Deal with overlap
;  03/12/2012: (KTS) Add option to use coadded weight masks instead of
;                    exp_nmiss.  Make weights the default.
;  03/21/2012: (KTS) Finalize procedure, re-name
;  03/26/2012: (KTS) use get_lps12_fieldinfo instead of get_lps12_map_npix
;  04/19/2012: (KTS) Change default coadd dir to maps/20120420/coadds/
;  06/28/2012: (KTS) Change default coadd dir to maps/20120620/coadds/
;;;


PRO make_apod_mask_lps12, field_idx, $
                          overlap=overlap, $
                          use_nmiss=use_nmiss, $
                          coadd_dir=coadd_dir, $
                          taper_arcmin=taper_arcmin, $
                          threshold=threshold, $
                          fwhm_smooth_arcmin=fwhm_smooth_arcmin, $
                          ret_apod=ret_apod, $
                          stopit=stopit
compile_opt IDL2, HIDDEN

; Output dir
savdir = '/home/kstory/lps12/masks/apod/'

; coadd dir (only if using weight map)
if ~keyword_set(use_nmiss) then begin
    if ~keyword_set(coadd_dir) then coadd_dir = '/data/kstory/projects/lps12/maps/20120620/coadds/'
endif

; Read command line arguments
if n_elements(taper_arcmin) eq 0 then taper_arcmin=60.
staper = sigfig(taper_arcmin,2)

if n_elements(threshold) eq 0 then threshold=0.05
sthresh = sigfig(threshold,3)

if n_elements(fwhm_smooth_arcmin) eq 0 then fwhm_smooth_arcmin=30.
ssmooth = sigfig(fwhm_smooth_arcmin,2)


;-------------------------
; Setup
;-------------------------

; Get field information
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]

field_name = fst.name
field_dir_name = fst.dir_name

print, 'MAKE_APOD_MASK_LPS12, field = ', field_name


;-------------------------
; Define the initial apod mask
;-------------------------
; There are two ways to define the apod mask:
; 1) exp(-nmiss)
; 2) weight map (DEFAULT)

; Option 2: exp(-nmiss)
if keyword_set(use_nmiss) then begin
    print, 'MAKE_APOD_MASK_LPS12, USE EXP(-NMISS) METHOD'

    restore, '/data/kstory/projects/lps12/masks/nmiss_sav/nmiss_'+field_name+'.sav' ; get nmiss
    apod = exp(-nmiss)

; Option 1: use the weight mask
endif else begin
    print, 'MAKE_APOD_MASK_LPS12, USE WEIGHT_MAP METHOD'

    coadd_name = coadd_dir+'coadd_'+field_name+'_50mJy.fits'
    d = expand_fits_struct(read_spt_fits(coadd_name))
    apod = d.weight.map
    apod /=max(apod)

endelse


;-------------------------
; We want to control the overlap between fields.  Use
; get_map_overlap() to define how much overlap in deg to allow.
;-------------------------

; get ra and dec for each pixel
npixels = (get_lps12_fieldinfo(field_idx)).npix
radec0  = [ fst.ra0, fst.dec0 ]
reso_arcmin = 1.0

pix2ang_proj5, npixels, radec0, reso_arcmin, ra, dec

; Deal with the fact that ra wraps:

; fields between 0h and 6h
if radec0[0] lt 180 then begin
    wh_neg = where(ra gt 180, nwh)
    if (nwh gt 0) then ra[wh_neg] = ra[wh_neg] - 360.
endif 

; fields between 21h and 24h
if radec0[0] gt 180 then begin
    wh_neg = where(ra lt 180, nwh)
    if (nwh gt 0) then ra[wh_neg] = ra[wh_neg] + 360.
endif 


dx = fst.dx
dy = fst.dy

; allow some overlap in regions
if ~keyword_set(overlap) then begin 
    overlap = get_lps12_map_overlap(field_idx)
endif

wh_outside = [-1]

ra_max=-1000 & ra_min=-1000 & dec_max=-1000 & dec_min=-1000

; ra max
if (overlap[0] ge 0) then begin
    ra_max  = radec0[0] + ( dx/2. + overlap[0])
    wh_good = where( ra lt ra_max , complement=wh_bad)
    wh_outside = union(wh_outside, wh_bad)
    print, 'apply ra_max = ', ra_max
endif

; ra min
if (overlap[1] ge 0) then begin
    ra_min  = radec0[0] - ( dx/2. + overlap[1])
    wh_good = where( ra gt ra_min , complement=wh_bad)
    wh_outside = union(wh_outside, wh_bad)
    print, 'apply ra_min = ', ra_min
endif

; dec max (most neg)
if (overlap[2] ge 0) then begin
    dec_max = radec0[1] - ( dy/2. + overlap[2])
    wh_good = where( abs(dec) lt abs(dec_max) , complement=wh_bad)
    wh_outside = union(wh_outside, wh_bad)
    print, 'apply dec_max = ', abs(dec_max)
endif

; dec min (least neg)
if (overlap[3] ge 0) then begin
    dec_min = radec0[1] + ( dy/2. + overlap[3])
    wh_good = where( abs(dec) gt abs(dec_min) , complement=wh_bad)
    wh_outside = union(wh_outside, wh_bad)
    print, 'apply dec_min = ', abs(dec_min)
endif

; clean up
if n_elements(wh_outside) gt 1 then begin
    wh_outside = wh_outside[1:*]
endif

; Set all regions outisde of the field region to zero
if (wh_outside[0] ne -1) then apod[wh_outside] = 0.


; Set the boarders to zero
npixzero = 20 ; number of pixels to zero at the edge
apod[0:npixzero-1,*]=0.0
apod[*,0:npixzero-1]=0.0
apod[*,n_elements(apod[0,*])-npixzero-1-1:n_elements(apod[0,*])-1]=0
apod[n_elements(apod[*,0])-npixzero-1-1:n_elements(apod[*,0])-1,*]=0

print, "Call hanningize_mask()"
apod = hanningize_mask(apod,taper_arcmin,threshold,fwhm_smooth_arcmin)


; --------------------------
; Print the total of the 'nominal area all equals 1' mask
ra_max0  = radec0[0] + ( dx/2.)
ra_min0  = radec0[0] - ( dx/2.)
dec_max0 = radec0[1] + ( dy/2.)
dec_min0 = radec0[1] - ( dy/2.)
wh_in_region0 = where( (ra gt ra_min0) and (ra lt ra_max0) and $
                      (dec gt dec_min0) and (dec lt dec_max0) , $
                    complement=wh_outside0)

apod0 = 0*apod
apod0[wh_in_region0] = 1
print, "ovelap = ", overlap
print, "max region = ", ra_max, ra_min, dec_max, dec_min
print, "Total of nominal area = ", total(apod0)
print, "Percent for this apod mask = ", total(apod) / total(apod0)
print, "nominal map region: ", ra_max0, ra_min0, dec_max0, dec_min0
; --------------------------

; Save the mask
;savename = savdir+'apod_test_'+field_name+'.sav'
savename=savdir+'apod_'+field_name+'_'+staper+'_'+sthresh+'_'+ssmooth+'.sav'
print, "MAKE_APOD_MASK_LPS12, output mask it file:"
print, "   "+savename
save,apod,filename=savename

; return the apodization mask
if keyword_set(ret_apod) then begin
    print, 'returning apod '
    ret_apod=apod
endif
if keyword_set(stopit) then stop
end




;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Reference code from RK, CR
;
;;;;;;;;;;;;;;;;;;;;;;;;;;

; ;-------------------------
; ; Get data from .sav files, from CR
; ;-------------------------
; if 0 then begin
; ; coadd directory:
;     dir = '/data11/cr/run4_2009/'+field_dir_name+'/coadds/'
    
; ; Restore the coadd
;     ifile=dir+'coadd_'+field_dir_name+'_150_nowt.sav'
;     restore,ifile
    
; ; Get the weight mask
;     wmask = map.weight
;     wmask/= max(wmask)
;     apod = wmask
; endif

; ;-------------------------
; ; Get data from .fits files, from RK
; ;-------------------------
; if 1 then begin
;     sband = '150'
;     coadd_name = coadd_dir+'coadd_'+field_name+'_50mJy.fits'
; ; ra5h30dec-55 had two different scans with different coverage.  use
; ; the more conservative coverage area.
;     if field eq 'ra5h30dec-55' then begin
;         coaddname='pre_29Apr2008_coadd_ra5h30dec-55_150.fits'
;     endif
;     d=read_spt_fits(coaddname)
;     apod=d.weight.map
; endif else begin
;     restore,'1.18/coadd_'+field+'.sav'
;     apod = exp(-nmiss)
; endelse


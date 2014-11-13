;;;
; NAME: make_apod_mask_0305
; PURPOSE:
;   Explore how to make apodization masks
;
; NOTES:
; 1) Now using the lpf-ds IDF's from Ryan, Jan3 2012
; 2) Need to check hanning choices
;
; MODIFICATION HISTORY:
;  02/10/2012: (KTS) Created from /home/cr/code/spt/coverage/mk_mask_run4_2009.pro
;  02/21/2012: (KTS) Add info from /home/rkeisler/ps09/make_apod.pro
;  03/06/2012: (KTS) add allow_overlap
;;;


pro make_apod_mask_0305, field_idx, $
                         allow_overlap=allow_overlap, $
                         taper_arcmin=taper_arcmin, $
                         threshold=threshold, $
                         fwhm_smooth_arcmin=fwhm_smooth_arcmin, $
                         stopit=stopit

; Read command line arguments
if n_elements(taper_arcmin) eq 0 then taper_arcmin=60.
staper = sigfig(taper_arcmin,2)

if n_elements(threshold) eq 0 then threshold=0.05
sthresh = sigfig(threshold,3)

if n_elements(fwhm_smooth_arcmin) eq 0 then fwhm_smooth_arcmin=40.
ssmooth = sigfig(fwhm_smooth_arcmin,2)


;-------------------------
; Setup
;-------------------------

; Get field information
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]

field_name = fst.name
field_dir_name = fst.dir_name

; get nmiss array:
restore, '/data/kstory/projects/lps12/masks/nmiss_sav/nmiss_'+field_name+'.sav'

;;; Coadds: Obsolete?
; coadd_dir = '/data/kstory/projects/lps12/maps/20120224/coadds/'
; coadd_name = coadd_dir+'coadd_'+field_name+'_50mJy.fits'
; ; ra5h30dec-55 had two different scans with different coverage.  use
; ; the more conservative coverage area.
; if field_idx eq 0 then begin
;     coadd_name = coadd_dir+'pre_29Apr2008_coadd_'+field_name+'_50mJy.fits'
; endif    

;-------------------------
; Read the coadd, get the weight map
; Obsolete?
;-------------------------
;d = read_spt_fits(coadd_name)
;apod = d.weight.map
;apod /= max(apod)


;-------------------------
; We want to ensure that there is no overlap between fields, so that
; no two apodization masks are non-zero in the same pixels on the
; sky.  Do this the most naive way; cut off any pixels that are
; outside of the pre-defined boundaries for each field.

; get ra and dec for each pixel
npixels = get_lps12_map_npix(field_idx)
radec0  = [ fst.ra0, fst.dec0 ]
reso_arcmin = 1.0

pix2ang_proj5, npixels, radec0, reso_arcmin, ra, dec

dx = fst.dx
dy = fst.dy

ra_min  = radec0[0] - dx/2.
ra_max  = radec0[0] + dx/2.
dec_min = radec0[1] - dy/2.
dec_max = radec0[1] + dy/2.
wh_in_region = where( (ra gt ra_min) and (ra lt ra_max) and $
                      (dec gt dec_min) and (dec lt dec_max) , $
                    complement=wh_outside)

apod = exp(-nmiss)

; Set all regions outisde of the field region to zero
if ~keyword_set(allow_overlap) then apod[wh_outside] = 0.

; Set the boarders to zero
npixzero = 20 ; number of pixels to zero at the edge
apod[0:npixzero-1,*]=0.0
apod[*,0:npixzero-1]=0.0
apod[*,n_elements(apod[0,*])-npixzero-1-1:n_elements(apod[0,*])-1]=0
apod[n_elements(apod[*,0])-npixzero-1-1:n_elements(apod[*,0])-1,*]=0

print, "Call hanningize_mask()"
apod = hanningize_mask(apod,taper_arcmin,threshold,fwhm_smooth_arcmin)

; Save the mask
ktsdir = '/home/kstory/lps12/masks/apod_0305/'
savename=ktsdir+'apod_'+field_name+'_'+staper+'_'+sthresh+'_'+ssmooth+'.sav'
if keyword_set(allow_overlap) then savename=ktsdir+'apod_allow_overlap_'+field_name+'_'+staper+'_'+sthresh+'_'+ssmooth+'.sav'
print, "MAKE_APOD_MASK_0305, output mask it file:"
print, "   "+savename
save,apod,filename=savename

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


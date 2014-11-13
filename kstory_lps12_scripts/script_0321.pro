;;;
; NAME: script_0320
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) sss, Study L-R leakage into 2010, 2011 maps
;
; MODIFICATION HISTORY:
;  03/20/2012: (KTS) Created
;;;


;...................................................................
; Script
PRO mk_kern, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name
print, 'Analyzing field ', field_idx, ', name = ', field_name

; inputs
reso_arcmin = 1.0
maxell = 3000
maskfile = '/home/kstory/lps12/masks/masks_50mJy/mask_'+field_name+'.sav'
restore, maskfile

kernel=coupling_kernel(mask, reso_arcmin, maxell=maxell, /changevar, $
                       interp=1000, oversamp=8, /cheby, ellkern=ellkern, $
                       curlyw=curlyw)

stop
END


;...................................................................
; make coadds for all fields fixed from az-wrap issue
PRO run_coadd_maps
compile_opt IDL2, HIDDEN

fields = [12, 9, 1, 2, 7, 15, 18]
nfields = n_elements(fields)

for ii=0, nfields-1 do begin

    idx = fields[ii]
    coadd_maps_0305, idx

endfor
END


;...................................................................
; make all apod masks
PRO run_apod_mask
compile_opt IDL2, HIDDEN

for ii=0, 19 do begin

    make_apod_mask_lps12, ii, $
      coadd_dir = '/data/kstory/projects/lps12/maps/20120305/coadds/'

endfor
END


;...................................................................
; make all ptsrc masks
PRO run_ptsrc_mask
compile_opt IDL2, HIDDEN

for ii=0, 19 do begin

    make_ptsrc_mask_lps12, ii, $
      config_dir = '/home/kstory/lps12/ptsrc_lists/', $
      savdir = '/home/kstory/lps12/masks/ptsrc/'

endfor
END



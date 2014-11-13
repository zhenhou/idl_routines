;;;
; NAME: script_0326
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  03/26/2012: (KTS) Created
;;;

;...................................................................
; make tweight jackknife runlists
PRO test_kern
compile_opt IDL2, HIDDEN

npix_pad = 4320
reso_arcmin = 1.0
maxell = 3000
maskfile = '/home/kstory/lps12/masks/masks_50mJy/mask_ra5h30dec-55_2008.sav'
restore, maskfile

; pad mask to full size
mask_padded = dblarr(npix_pad, npix_pad)
nx = (size(mask))[1]
ny = (size(mask))[2]
mask_padded[0:nx-1, 0:ny-1] = mask

kernel=coupling_kernel(mask_padded, reso_arcmin, maxell=maxell, /changevar, $
                       interp=1000, oversamp=8, /cheby, ellkern=ellkern, $
                       curlyw=curlyw)
END

;...................................................................
; make tweight jackknife runlists
PRO make_kerns
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()

for ii=0, 19 do begin
    fst = f[ii]
    maskfile = '/home/kstory/lps12/masks/masks_50mJy/mask_'+fst.name+'.sav'
    savfile = '/data/kstory/projects/lps12/masks/masks_50mJy/kern_'+fst.name+'.sav'

    print, "Make coupling Kernel for field ", ii, ', '+fst.name
    make_coupling_kernel, ii, maskfile=maskfile, savfile=savfile
endfor

END

;...................................................................
; run all jacks
PRO run_lr_jack
compile_opt IDL2, HIDDEN
for ii=3, 19 do lps12_jack, ii, 'lr', savdir='/home/kstory/lps12/jacks/'
END

PRO run_12_jack
compile_opt IDL2, HIDDEN
for ii=0, 19 do lps12_jack, ii, '12', savdir='/home/kstory/lps12/jacks/'
END

PRO run_azrms_jack
compile_opt IDL2, HIDDEN
for ii=0, 19 do lps12_jack, ii, 'azrms', savdir='/home/kstory/lps12/jacks/'
END

PRO run_tweight_jack
compile_opt IDL2, HIDDEN
for ii=0, 19 do lps12_jack, ii, 'tweight', savdir='/home/kstory/lps12/jacks/'
END

PRO run_moon_jack
compile_opt IDL2, HIDDEN
for ii=0, 19 do lps12_jack, ii, 'moon', savdir='/home/kstory/lps12/jacks/'
END

PRO run_sun_jack
compile_opt IDL2, HIDDEN
list = [4,6,11,17,18,19]
for idx=0, n_elements(list)-1 do begin
    ii = list[idx]
    lps12_jack, ii, 'sun', savdir='/home/kstory/lps12/jacks/'
endfor
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                 ;
; Things for ra3h30dec-42.5 field ;
;                                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;...................................................................
; make tweight jackknife runlists
PRO runlist_17
compile_opt IDL2, HIDDEN

; make runlist
make_runlist_for_lps12_single_field, 17

END


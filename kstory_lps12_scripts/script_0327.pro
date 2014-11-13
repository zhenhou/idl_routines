;;;
; NAME: script_0327
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
save, kernel, filename='/home/kstory/lps12/masks/masks_50mJy/test_kern_ra5h30dec-55_2008.sav'
END

;...................................................................
; make tweight jackknife runlists
PRO make_kerns_2to10
compile_opt IDL2, HIDDEN
f = lps12_fieldstruct()

for ii=2, 10 do begin
    fst = f[ii]
    maskfile = '/home/kstory/lps12/masks/masks_50mJy/mask_'+fst.name+'.sav'
    savfile = '/data/kstory/projects/lps12/masks/masks_50mJy/kern_'+fst.name+'.sav'

    print, "Make coupling Kernel for field ", ii, ', '+fst.name
    make_coupling_kernel, ii, maskfile=maskfile, savfile=savfile
endfor
END

PRO make_kerns_11to19
compile_opt IDL2, HIDDEN
f = lps12_fieldstruct()

for ii=11, 19 do begin
    fst = f[ii]
    maskfile = '/home/kstory/lps12/masks/masks_50mJy/mask_'+fst.name+'.sav'
    savfile = '/data/kstory/projects/lps12/masks/masks_50mJy/kern_'+fst.name+'.sav'

    print, "Make coupling Kernel for field ", ii, ', '+fst.name
    make_coupling_kernel, ii, maskfile=maskfile, savfile=savfile
endfor
END

PRO make_kerns_0
compile_opt IDL2, HIDDEN
f = lps12_fieldstruct()

ii = 0
fst = f[ii]
maskfile = '/home/kstory/lps12/masks/masks_50mJy/mask_'+fst.name+'.sav'
savfile = '/data/kstory/projects/lps12/masks/masks_50mJy/kern_'+fst.name+'.sav'

print, "Make coupling Kernel for field ", ii, ', '+fst.name
make_coupling_kernel, ii, maskfile=maskfile, savfile=savfile
END

;...................................................................
; run all jacks

; local6
PRO run_jacks_local6_0327
compile_opt IDL2, HIDDEN
run_12_jack
run_tweight_jack
END

;; local7
PRO run_jacks_local7_0327
compile_opt IDL2, HIDDEN
run_azrms_jack
run_moon_jack
run_sun_jack
END


PRO run_lr_jack
compile_opt IDL2, HIDDEN
for ii=4, 19 do lps12_jack, ii, 'lr', savdir='/home/kstory/lps12/jacks/'
END

PRO run_12_jack
compile_opt IDL2, HIDDEN
;for ii=0, 19 do lps12_jack, ii, '12', savdir='/home/kstory/lps12/jacks/'
for ii=1, 19 do lps12_jack, ii, '12', savdir='/home/kstory/lps12/jacks/'
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


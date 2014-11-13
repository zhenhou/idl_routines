;;;
; NAME: script_0401
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) test_distr; prob of getting 0.0012 from 120 draws from uniform distribution
;
; MODIFICATION HISTORY:
;  04/01/2012: (KTS) Created
;;;

;...................................................................
; 
PRO sss
compile_opt IDL2, HIDDEN

;coupling kernels
restore, '/data17/rkeisler/ps09/end_ra5h30dec-55_1.24_kweight.sav'
mask_rk = kernel
mask_rk = mask

mask_dir = '/data/kstory/projects/lps12/masks/masks_50mJy/'
mask_file = mask_dir+'mask_ra5h30dec-55_2008.sav'
restore, mask_file
mask_ks = mask

; kern_file = kern_dir+'kern_ra5h30dec-55_2008.sav'
; restore, kern_file
; kern_ks = kernel
; ellkern_ks = ellkern

stop
END

;...................................................................
; make new coupling kernel
PRO coupling_kernel_field0_0403
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()
fst = (f[0])
field_name = fst.name

make_coupling_kernel, 0, $
  maskfile='/home/kstory/lps12/masks/masks_50mJy/mask_'+field_name+'.sav', $
  savfile='/data/kstory/projects/lps12/masks/masks_50mJy/kern_'+field_name+'.sav'

END

;...................................................................
; test coupling kernel with different pad length
PRO coupling_kernel_test_field0_0403
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()
fst = (f[0])
field_name = fst.name

make_coupling_kernel, 0, $
  maskfile='/home/kstory/lps12/masks/masks_50mJy/mask_'+field_name+'.sav', $
  savfile='/data/kstory/projects/lps12/masks/masks_50mJy/kern_smallPad_'+field_name+'.sav', $
  test_pad=2160

END


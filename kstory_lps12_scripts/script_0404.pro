;;;
; NAME: script_0404
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) coupling kernels
;
; MODIFICATION HISTORY:
;  04/04/2012: (KTS) Created
;;;

;...................................................................
; 
PRO sss
compile_opt IDL2, HIDDEN

;coupling kernels
;restore, '/data17/rkeisler/ps09/end_ra5h30dec-55_1.24_kweight.sav'
restore, '/data/kstory/projects/lps12/scripts/sav_files/mask_rk_ra5h30dec-55.sav'
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
; test coupling kernel with different pad length
PRO coupling_kernel_rk
compile_opt IDL2, HIDDEN

maskfile = '/data/kstory/projects/lps12/scripts/sav_files/mask_rk_ra5h30dec-55.sav'

f = lps12_fieldstruct()
fst = (f[0])
field_name = fst.name

make_coupling_kernel, 0, $
  maskfile=maskfile, $
  savfile='/data/kstory/projects/lps12/masks/masks_50mJy/kern_rk_'+field_name+'.sav', $
  test_pad=2160

END

;...................................................................
; test coupling kernels RK to KS
PRO kern_check
compile_opt IDL2, HIDDEN

; re-make RK kernel
maskfile1 = '/data/kstory/projects/lps12/scripts/sav_files/mask_rk_ra5h30dec-55.sav'
make_coupling_kernel, 0, $
  maskfile=maskfile1, $
  savfile='/data/kstory/projects/lps12/masks/masks_50mJy/kern_rk3_ra5h30dec-55.sav', $
  test_pad=2160

; re-make KS kernel
maskfile2 = '/data/kstory/projects/lps12/masks/masks_50mJy/mask_ra5h30dec-55_2008.sav'
make_coupling_kernel, 0, $
  maskfile=maskfile2, $
  savfile='/data/kstory/projects/lps12/masks/masks_50mJy/kern_ks3_ra5h30dec-55.sav', $
  test_pad=2160

END


;...................................................................
; test coupling kernels RK to KS
PRO plot_kerns
compile_opt IDL2, HIDDEN

figdir = '/data/kstory/projects/lps12/figs/'

; re-make RK kernel
maskfile1 = '/data/kstory/projects/lps12/scripts/sav_files/mask_rk_ra5h30dec-55.sav'
restore, maskfile1
mask_rk = mask
tv_spt_map, mask_rk[0:959, 0:959], /norms, title="RK mask: ra5h30dec-55"
err = tvread(/png, filename=figdir+'mask_rk_0404', /nodialog)

; re-make KS kernel
maskfile2 = '/data/kstory/projects/lps12/masks/masks_50mJy/mask_ra5h30dec-55_2008.sav'
restore, maskfile2
mask_ks = mask
tv_spt_map, mask_ks[0:959, 0:959], /norms, title="KTS mask: ra5h30dec-55"
err = tvread(/png, filename=figdir+'mask_ks_0404', /nodialog)

tv_spt_map, mask_rk[0:959,0:959] - mask_ks[0:959,0:959], title="RK - KTS, ra5h30dec-55", /norms
err = tvread(/png, filename=figdir+'mask_diff_0404', /nodialog)

stop
END


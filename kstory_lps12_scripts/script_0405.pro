;;;
; NAME: script_0405
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) remake_kerns_0405, remake all coupling kernels with correct
; padding, renormalized pad size, etc.
;
; MODIFICATION HISTORY:
;  04/05/2012: (KTS) Created
;;;

;...................................................................
; 
PRO remake_kerns_0405
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()

for ii=0, 19 do begin

    fst = f[ii]
    field_name = fst.name

    ; re-make KS kernel
    maskfile = '/data/kstory/projects/lps12/masks/masks_50mJy/mask_'+field_name+'.sav'
    make_coupling_kernel, ii, maskfile=maskfile, $
      savfile='/data/kstory/projects/lps12/masks/masks_50mJy/kern_'+field_name+'.sav'

endfor    

END

;...................................................................
; test coupling kernels RK to KS
PRO kern_check
compile_opt IDL2, HIDDEN

; re-make RK kernel
maskfile1 = '/data/kstory/projects/lps12/scripts/sav_files/mask_rk_ra5h30dec-55.sav'
make_coupling_kernel, 0, $
  maskfile=maskfile1, $
  savfile='/data/kstory/projects/lps12/masks/masks_50mJy/kern_rk5_ra5h30dec-55.sav', $
  test_pad=2160

; re-make KS kernel
maskfile2 = '/data/kstory/projects/lps12/masks/masks_50mJy/mask_ra5h30dec-55_2008.sav'
make_coupling_kernel, 0, $
  maskfile=maskfile2, $
  savfile='/data/kstory/projects/lps12/masks/masks_50mJy/kern_ks5_ra5h30dec-55.sav', $
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


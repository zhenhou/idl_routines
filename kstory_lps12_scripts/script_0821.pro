;;;
; NAME: script_0815
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Pipeline test2
;;;

FUNCTION get_theory
;;; get dl_th
; readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_uK2
; dl_uK2_2 = dl_uK2[50:3998]
; l_vec_2  = l_vec[50:3998]
; dl_th_2 = dl_uK2_2 # wf2use
theoryspectrumfiles=[ ['/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat',$
                       '/home/cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt',$
                       '/home/cr/cmb_models/foreground_sim09_150.txt']]

; get foregrounds
readcol,'/home/cr/cmb_models/foreground_sim09_150.txt', f_ell, f_dl
readcol,'~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt', ksz_ell, ksz_dl
readcol, '/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat', sim_ell, sim_dl

wh = indgen(15000)
dl_tot = sim_dl[wh] + f_dl[wh] + ksz_dl[wh]

RETURN, dl_tot
END


;;;;;;;;;;;;;;;;;;;;;;;
;;; 2) Pipeline test2, spec1
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_spec1
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_08p2_spec1_kweight_calib.sav'

; get pipe2
restore, file
calib = 1.
;calib = 0.825
;calib = 0.76
calib = 0.85
dl_p2 = spectrum*(1d12)^2. * (calib)
cov_p2 = cov

; get final combined
restore, edir+'run_08/combined_spectrum_20120717_174249_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

; readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_uK2
; dl_uK2_2 = dl_uK2[50:3998]
; l_vec_2  = l_vec[50:3998]
; dl_th_2 = dl_uK2_2 # wf2use

;----------------------------
;Get the theory spectrum
dl_th = get_theory()
x_th = indgen(3001) + 500
dl_th_plot = dl_th[x_th]

dl_th_2 = dl_th[50:3998]
dl_th_2 = reform(dl_th_2 # wf2use)

; Plotting
vec=indgen(47) + 9

window, 1
plot, l[vec], dl_p2[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl
oplot, l[vec], dl_p2[vec], color=!red, linestyle=0
;oplot, l[vec], dl_all[vec]

oplot, l[vec], dl_th_2[vec], color=!green, linestyle=0
oplot, x_th, dl_th_plot, color=!blue, linestyle=2

window, 2
plot, l[vec], (dl_th_2/dl_p2)[vec]
stop
END



;;;;;;;;;;;;;;;;;;;;;;;
;;; 2) Pipeline test2, spec2
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_spec2
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_08p2_spec2_kweight_calib.sav'

; get pipe2
restore, file
calib = 1.
;calib = 0.825
;calib = 0.76
;calib = 0.85
dl_p2 = spectrum*(1d12)^2. * (calib)
cov_p2 = cov

; get spec1
restore, edir+'end_ra3h30dec-60_08p2_spec1_kweight_calib.sav'
dl_p2_spec1 = spectrum*(1d12)^2. * (calib)

; get final combined
restore, edir+'run_08/combined_spectrum_20120717_174249_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

; readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_uK2
; dl_uK2_2 = dl_uK2[50:3998]
; l_vec_2  = l_vec[50:3998]
; dl_th_2 = dl_uK2_2 # wf2use

;----------------------------
;Get the theory spectrum
dl_th = get_theory()
x_th = indgen(3001) + 500
;x_th = x_th - 10 ; shift back 10
dl_th_plot = dl_th[x_th]

dl_th_2 = dl_th[50:3998]
dl_th_2 = reform(dl_th_2 # wf2use)

; Plotting
vec=indgen(47) + 9

window, 1
plot, l[vec], dl_p2[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl
oplot, l[vec], dl_p2[vec], color=!red, linestyle=0
;oplot, l[vec], dl_all[vec]

;plot, l[vec], dl_p2_spec1[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl
oplot, l[vec], dl_p2_spec1[vec], color=!blue, linestyle=1

; oplot, l[vec], dl_th_2[vec], color=!green, linestyle=0
; oplot, x_th, dl_th_plot, color=!purple, linestyle=2

window, 2
plot, l[vec], (dl_th_2/dl_p2)[vec]
stop
END

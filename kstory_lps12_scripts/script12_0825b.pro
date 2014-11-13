;;;
; NAME: script12_0825b
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Pipeline test2
;;;

FUNCTION get_theory, l_tot=l_tot, dl_tot=dl_tot
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

vec = indgen(15000)
wh_f = where(f_ell eq sim_ell[0])
wh_ksz = where(ksz_ell eq sim_ell[0])

l_tot = sim_ell
dl_tot = sim_dl[vec] + f_dl[vec+wh_f] + ksz_dl[vec+wh_ksz]
stop
RETURN, 1
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
;calib = 0.85
dl_p1 = spectrum*(1d12)^2. * (calib)
cov_p1 = cov

; get final combined
restore, edir+'run_08/combined_spectrum_20120717_174249_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

;----------------------------
;Get the theory spectrum 1
dl_th1 = get_theory()
dl_th1b = dl_th1[50:3998]
dl_th1b = reform(dl_th1b # wf2use)


;----------------------------
;Get the theory spectrum 2
readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec2,dl_th2
dl_th2b = dl_th2[50:3998]
l_vec2b  = l_vec2[50:3998]

dl_th2b = shift(dl_th2b, 10) ; shift by 10:
dl_th2b = dl_th2b # wf2use



;----------------------------
; Plotting
;----------------------------

vec=indgen(47) + 9

window, 1
plot, l[vec], dl_p1[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl
oplot, l[vec], dl_p1[vec], color=!red, linestyle=0
;oplot, l[vec], dl_all[vec]

oplot, l[vec], dl_th1b[vec], color=!green, linestyle=0
;oplot, x_th, dl_th_plot, color=!blue, linestyle=2

window, 2
plot, l[vec], (dl_th1b/dl_p1)[vec]
stop
END





;;;;;;;;;;;;;;;;;;;;;;;
;
; 2) Pipeline test2, spec1
;
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

; get final combined
restore, edir+'run_08/combined_spectrum_20120717_174249_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

;----------------------------
;Get the theory spectrum 1
dl_th1 = get_theory()
dl_th1b = dl_th1[50:3998]
dl_th1b = reform(dl_th1b # wf2use)


;----------------------------
;Get the theory spectrum 2
; readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec2,dl_th2
; dl_th2b = dl_th2[50:3998]
; l_vec2b  = l_vec2[50:3998]

; dl_th2b = shift(dl_th2b, 10) ; shift by 10:
; dl_th2b = dl_th2b # wf2use



;----------------------------
; Plotting
;----------------------------

vec=indgen(47) + 9

window, 1
plot, l[vec], dl_p2[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl
oplot, l[vec], dl_p2[vec], color=!red, linestyle=0
;oplot, l[vec], dl_all[vec]

oplot, l[vec], dl_th1b[vec], color=!green, linestyle=0
;oplot, x_th, dl_th_plot, color=!blue, linestyle=2

window, 2
plot, l[vec], (dl_th1b/dl_p2)[vec]
stop
END






;;;;;;;;;;;;;;;;;;;;;;;
;
; 2) Pipeline test2, all
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_all
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'

; get spec1
restore, edir+'end_ra3h30dec-60_08p2_spec1_kweight_calib.sav'
calib = 1.
dl_p1 = spectrum*(1d12)^2. * (calib)
cov_p1 = cov

; get spec2
restore, edir+'end_ra3h30dec-60_08p2_spec2_kweight_calib.sav'
calib = 1.
dl_p2 = spectrum*(1d12)^2. * (calib)
cov_p2 = cov

; get final combined
restore, edir+'run_08/combined_spectrum_20120717_174249_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

;----------------------------
;Get the theory spectrum 1
dl_th1 = get_theory()
dl_th1b = dl_th1[50:3998]
dl_th1b = reform(dl_th1b # wf2use)

;----------------------------
;Get the theory spectrum 2
readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec2,dl_th2
dl_th2b = dl_th2[50:3998]
l_vec2b  = l_vec2[50:3998]
dl_th2b = shift(dl_th2b, -10) ; shift by 10:
dl_th2b = dl_th2b # wf2use



;----------------------------
; Plotting
;----------------------------

vec=indgen(47) + 9

window, 1
plot, l[vec], dl_p1[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl

; spec1
oplot, l[vec], dl_p1[vec], color=!red, linestyle=0,thick=2

; spec2
oplot, l[vec], dl_p2[vec], color=!darkgreen, linestyle=0,thick=2

;oplot, l[vec], dl_all[vec]

oplot, l[vec], dl_th1b[vec], color=!purple, linestyle=2,thick=2
oplot, l[vec], dl_th2b[vec], color=!blue, linestyle=2,thick=2

legend,['spec1','spec2','th1','th2'],color=[!red,!darkgreen,!purple,!blue],linestyle=[0,0,2,2],$
  pos=[1800,4000]

ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_all_0825'
err=tvread(/png,/nodialog,filename=ff)
stop
END











;;;;;;;;;;;;;;;;;;;;;;;
;
; 2) Pipeline test2, spec2, OBSOLETE
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_spec2_BAK
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

readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec1,dl_th1
dl_th1b = dl_th1[50:3998]
l_vec1b  = l_vec1[50:3998]
dl_th1b = dl_th1b # wf2use

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

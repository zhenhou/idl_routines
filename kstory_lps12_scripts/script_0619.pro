;;;
; NAME: script_0619
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Plots
;
; MODIFICATION HISTORY:
;  06/19/2012: (KTS) Created
;;;




;;;;;;;;;;;;;;;;;;;;;;
; Plot error bars
;;;;;;;;;;;;;;;;;;;;;;

PRO plot_ps

; plotting setup
if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
whplot = indgen(45) + 10

;run_07
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_144542_kweight.sav'
dl_07 = dl_all
diag07 = diag_nobeam

;run_05
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120614_192137_kweight.sav' ; Should be final
dl_05  = dl_all
diag05 = diag_nobeam

; run_07; 0809 only
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120619_195218_kweight_0809.sav'
dl_0809     = dl_all
diag_0809 = diag_nobeam
dd_e2e = dblarr(58)
for i=0, 57 do dd_e2e[i] = sqrt(cov_all_nobeam_e2e[i,i])

restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_205039_kweight_0809.sav' ; from comb2_0809
;dl2_0809     = dl_all_init
diag2_0809 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_rk = diag_nobeam

; cl_theory
;readcol,'/home/kstory/lps12/theory_spectrum/Cls_100sims_ave_from_combined_map.txt',l_vec,cl_uK2
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)

;----------------
; Plots
;----------------
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

; ; ps plot
; window, 0
; plot,l[whplot],dl_07[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
;   xtitle=xtitle,ytitle=ytitle,chars=chars, title='Full 2500 deg^2, run_07'
; oplot,l_vec,dl_th,color=!red
; errplot,l[whplot],(dl_07-diag07)[whplot]*1d12,(dl_07+diag07)[whplot]*1d12
; ; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0619')

; ; compare 2500 against K11
; window, 1
; plot, l[whplot], ((err_rk )/diag07)[whplot], xtitle=xtitle, ytitle='diag_k11/diag07', title='Error bars, k11 v.s. run_07'
; oplot, [0,10000], [sqrt(2544./790),sqrt(2544./790)], linestyle=1
; ;err=tvread(/png,/nodialog,filename=fdir+'err_k11v07_0619')

; compare 0809 against K11
window, 2
plot, l[whplot], ((err_rk )/diag_0809)[whplot], xtitle=xtitle, ytitle='diag_k11/diag_0809', title='Error bars, comb_0809'
oplot, [0,10000], [1,1], linestyle=1
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v0809_0619')

;window, 3
;plot, l[whplot], ((err_rk )/diag2_0809)[whplot], xtitle=xtitle, ytitle='diag_k11/comb_0809', title='Error bars, comb_0809'

; window, 3
; plot, l[whplot], (dl_05/dl_07)[whplot], xtitle=xtitle, ytitle='dl_05 / dl_07'

; window, 4
; plot, l[whplot], (diag05/diag07)[whplot], xtitle=xtitle, ytitle='dl_05 / dl_07'


stop
END



;;;;;;;;;;;;;;;;;;;;
; Investigate ra3h30dec-60 errors
;;;;;;;;;;;;;;;;;;;;
PRO err_4
; get RK
; restore, '/data/rkeisler/ps09/ra3h30dec-60_diag_k11_for_kyle.sav'
; rk_diag_meas = diag_meas
; rk_diag_sv   = diag_sample
; rk_diag      = diag
restore, '/home/rkeisler/ps09/1.37/end_ra3h30dec-60_1.37_kweight.sav'
rk_diag_meas = dblarr(58)
rk_diag_sv   = dblarr(58)
rk_diag      = dblarr(58)
rk_dl        = spectrum
for i=0, 57 do begin
    rk_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    rk_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    rk_diag[i]      = sqrt(cov[i,0,i,0])
endfor

; get KS
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07_kweight.sav'
ks_diag_meas = dblarr(58)
ks_diag_sv   = dblarr(58)
ks_diag      = dblarr(58)
ks_dl        = spectrum
for i=0, 57 do begin
    ks_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks_diag[i]      = sqrt(cov[i,0,i,0])
endfor

; Try this:
; ks_diag_meas *= 0.825
; ks_diag_sv *= sqrt(0.825)
; ;ks_diag *= sqrt(0.825)
; ks_diag = sqrt(ks_diag_meas^2. + ks_diag_sv^2.)

; run_07; 0809 only
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120619_195218_kweight_0809.sav'
diag_0809 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
k11_diag = diag_nobeam

;---------------------------
; PLOTS
;---------------------------
whplot = indgen(45) + 10
fdir = '/home/kstory/public_html/notebook/spt_lps12/'
l = banddef-25.

window, 0, xsize=800, ysize=500
plot, l[whplot], (rk_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, ra3h30dec-60', $ ; tot
  ytitle='err_rk / err_ks', xtitle='ell'
oplot, l[whplot], (rk_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (rk_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
oplot, l[whplot], (k11_diag/diag_0809)[whplot], color=!red, linestyle=3       ; 0809
;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0619')

stop
END


PRO err_check
; get KS 07
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07_kweight.sav'
ks_diag_sv_07   = dblarr(58)
for i=0, 57 do begin
    ks_diag_sv_07[i]   = sqrt(sample_cov[i,0,i,0])
endfor

; get KS 05
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_05_kweight.sav'
ks_diag_sv_05   = dblarr(58)
for i=0, 57 do begin
    ks_diag_sv_05[i]   = sqrt(sample_cov[i,0,i,0])
endfor

; PLOTS
whplot = indgen(45) + 10
l = banddef-25.
window, 0, xsize=800, ysize=500
plot, l[whplot], (ks_diag_sv_05/ks_diag_sv_07)[whplot], yr=[.7,1.3],/yst,title='err_SV, 05 v 07'
stop
END


;;;;;;;;;;;;;;;;
; Compare kweights between K11 and lps12, ra3h30dec-60
;;;;;;;;;;;;;;;;
PRO err_check2
; get KS 07
;restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07_kweight.sav'
weight_ks = shift(get_lps12_kweight(4),2160,2160)

restore, '/home/rkeisler/ps09/1.37/end_ra3h30dec-60_1.37_kweight.sav'
weight_rk_orig = shift(weight_2D,1080,1080)

weight_rk = congrid(weight_rk_orig, 4320,4320)

xx = weight_ks / weight_rk
stop
END



;;;;;;;;;;;;;;;;;;;;
; Plot chisq table
;;;;;;;;;;;;;;;;;;;;

PRO chisq_table
arr = [ $
0.1429, 0.5913, 0.5523, 0.2829, 0.8622, $
0.9355, 0.2248, 0.6230, 0.7867, 0.0621, $
0.6158, 0.5241, 0.4908, 0.0668, 0.7850, $
0.6743, 0.7568, 0.6104, 0.1309, 0.4049, $
0.2484, 0.2599, 0.3400, 0.5380, 0.4459, 0.1867, $
0.2887, 0.7898, 0.7422, 0.0519, 0.8678, $
0.0398, 0.7270, 0.2814, 0.4490, 0.3418, 0.5987, $
0.9234, 0.1911, 0.3487, 0.1353, 0.9806, $
0.0326, 0.4760, 0.4829, 0.7768, 0.5356, $
0.5201, 0.1510, 0.0528, 0.4078, 0.3358, $
0.0212, 0.7528, 0.4710, 0.6239, 0.6395, $
0.0631, 0.1233, 0.1247, 0.9935, 0.6293, 0.5808, $
0.4802, 0.9686, 0.5769, 0.2053, 0.3233, $
0.7478, 0.8244, 0.6691, 0.9150, 0.3582, $
0.5175, 0.5726, 0.9999, 0.8742, 0.1721, $
0.3367, 0.6083, 0.0166, 0.2008, 0.3649, $
0.1122, 0.4916, 0.8304, 0.2184, 0.1433, $
0.4490, 0.0620, 0.0776, 0.1332, 0.6112, 0.6210, $
0.9811, 0.3201, 0.8993, 0.2339, 0.1405, 0.2465, $
0.4685, 0.6034, 0.6627, 0.3953, 0.7867, 0.9786 ]

hh = histogram(arr, nbins=23, locations=xx, min=-0.05, max=1.05)
plot, xx+0.025, hh, psym=10, xr=[-0.05,1.05],/xst,$
  xtitle='PTE', ytitle='num. jacks', title='PTE distribution for all individual-field jacks'

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
err=tvread(/png,/nodialog,filename=fdir+'pte_run07_0619')
stop
END




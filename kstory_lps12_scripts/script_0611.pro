;;;
; NAME: script_0611
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) end2end_no_kweight
; 2) end2end_calib
; 3) plot error bars
;
; MODIFICATION HISTORY:
;  06/11/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;
; end2end
;;;;;;;;;;;;;;;;;;;;;;

PRO end2end_no_kweight
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END


PRO end2end_calib
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume, /use_kweight, calib=0.680

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END


;;;;;;;;;;;;;;;;;;;;;;
; Plot error bars
;;;;;;;;;;;;;;;;;;;;;;

PRO plot_ps_run_05

; plotting setup
if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
whplot = indgen(45) + 10

;run_05
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120610_190210_kweight.sav' ; calib applied to cov
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120611_175439_kweight.sav' ; no calib on cov
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120611_225834_kweight.sav' ; new cov_sv
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120611_231531_kweight.sav' ; fix by-year
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120611_235136_kweight.sav' ; fix beam errors
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_022240_kweight.sav' ; to be sure.
dl_05  = dl_all
diag05 = diag_nobeam
diag05_e2e = dblarr(58)
cov_e2e = cov_all_nobeam_data + cov_all_nobeam_mc
for ii=0, 57 do diag05_e2e[ii] = sqrt(cov_e2e[ii, ii])

; run_05; 0809 only
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_024104_kweight_0809.sav'; best guess
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_025500_kweight_0809.sav' ; no combined_calib applied to cov_data
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_040328_kweight_0809.sav' ; corrected cov_sv
dl_0809     = dl_all
diag05_0809 = diag_nobeam

; run_04
restore, '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav'
diag04 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_rk = diag_nobeam


; cl_theory
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)

;----------------
; Plots
;----------------
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

; ps plot
window, 0
plot,l[whplot],dl_05[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars, title='Full 2500 deg^2, run_05'
oplot,l_vec,dl_th,color=!red
errplot,l[whplot],(dl_05-diag05)[whplot]*1d12,(dl_05+diag05)[whplot]*1d12
; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0611a')

window, 1
plot,l[whplot],dl_all[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars, title='2008, 2009 data, run_05'
oplot,l_vec,dl_th,color=!red
errplot,l[whplot],(dl_all-diag_nobeam)[whplot]*1d12,(dl_all+diag_nobeam)[whplot]*1d12
; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0611a')


; ;----------------
; ; make err plot comparing run 04 to 05

; window, 1
; plot, l[whplot], ((diag04 - diag05)/diag05)[whplot], xtitle=xtitle, ytitle='(diag04 - diag05)/diag05', $
;   title='Error bars, 04 vs 05'
; ;err=tvread(/png,/nodialog,filename=fdir+'err_04v05_0611a')

; ; compare against K11
window, 2
plot, l[whplot], ((err_rk )/diag05)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05', title='Error bars, K11 vs 05'
oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v05_0611a')

; window, 3
; plot, l[whplot], ((err_rk )/diag05_e2e)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05', title='Error bars, K11 vs 05_e2e'
; oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
; ;err=tvread(/png,/nodialog,filename=fdir+'err_k11v05e2e_0611a')

; window, 4
; plot, l[whplot], (diag05/diag05_e2e)[whplot], xtitle=xtitle, ytitle='diag05/diag05_e2e', title='Error bars, 05_sim vs 05_e2e'
; oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
; ;err=tvread(/png,/nodialog,filename=fdir+'err_05simv05e2e_0611a')

window, 5
plot, l[whplot], (err_rk/diag05_0809)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05_0809', title='Error bars, K11 vs 05_0809'
stop
END

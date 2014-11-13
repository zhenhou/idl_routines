;;;
; NAME: script_0613
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Plots
; 2) processing
;
; MODIFICATION HISTORY:
;  06/13/2012: (KTS) Created
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

;run_05
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_022240_kweight.sav' ; to be sure.
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_205636_kweight.sav' fix bug in cov_data*calib
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_213028_kweight.sav'
dl_05  = dl_all
diag05 = diag_nobeam
diag05_e2e = dblarr(58)
cov_e2e = cov_all_nobeam_data + cov_all_nobeam_mc
for ii=0, 57 do diag05_e2e[ii] = sqrt(cov_e2e[ii, ii])

;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_202229_kweight_calib.sav' ; e2e_calib
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_213316_kweight_calib.sav'
dl05_calib = dl_all
diag05_calib = diag_nobeam

; run_05; 0809 only
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_024104_kweight_0809.sav'; best guess
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_025500_kweight_0809.sav' ; no combined_calib applied to cov_data
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_040328_kweight_0809.sav' ; corrected cov_sv
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_222432_kweight_0809.sav' ; use cov_sv
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_222844_kweight_0809.sav' ; use cov from e2e
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

; window, 6
; plot,l[whplot],dl05_calib[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
;   xtitle=xtitle,ytitle=ytitle,chars=chars, title='Full 2500 deg^2, run_05'
; oplot,l_vec,dl_th,color=!red
; errplot,l[whplot],(dl05_calib-diag05_calib)[whplot]*1d12,(dl05_calib+diag05_calib)[whplot]*1d12
; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0611a')

; compare calib to debugged code
window, 1
plot, l[whplot], (diag05/diag05_calib)[whplot], xtitle=xtitle, ytitle='diag05/diag05_calib', $
  title='Error bars, 05 v 05_calib'
;err=tvread(/png,/nodialog,filename=fdir+'err_04v05_0611a')

; ; compare against K11
window, 2
plot, l[whplot], ((err_rk )/diag05)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05', title='Error bars, K11 vs 05'
oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v05_0613a')

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












;;;;;;;;;;;;;;;;;;;;;;
; End2End, run_06
;;;;;;;;;;;;;;;;;;;;;;

PRO end2end_run_06
tbegin = systime(0, /seconds)
;for idx=0, 19 do begin
for idx=11, 19 do begin
    print, 'run end2end_06, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_06, idx, /resume, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END




;;;;;;;;;;;;;;;;;;;;;;
; Jackknives, run_06
;;;;;;;;;;;;;;;;;;;;;;

PRO jacks_1
; ; azrms ; NEED jackgoodfiles
; for idx=0, 19 do begin
;     print, 'jack azrms, field ', idx
;     t0 = systime(0, /seconds)
;     lps12_jack, idx, 'azrms', run='06'
;     print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
; endfor

; tweight
for idx=3, 19 do begin
    print, 'jack tweight, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'tweight', run='06', savdir='/data18/kstory/lps12/jacks/run_06/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; 12
for idx=0, 19 do begin
    print, 'jack 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12', run='06', savdir='/data18/kstory/lps12/jacks/run_06/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


PRO jacks_2
; lr
;for idx=0, 19 do begin
for idx=11, 19 do begin
    print, 'jack lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr', run='06', savdir='/data18/kstory/lps12/jacks/run_06/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; moon
for idx=0, 19 do begin
    print, 'jack moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon', run='06', savdir='/data18/kstory/lps12/jacks/run_06/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; sun
idx_list = [4, 6, 11, 17, 18, 19]
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun', run='06', savdir='/data18/kstory/lps12/jacks/run_06/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


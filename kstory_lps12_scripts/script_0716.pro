;;;
; NAME: script_0716
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Test SV averaging
;
; MODIFICATION HISTORY:
;  07/16/2012: (KTS) Created
;;;

pro plot_err
; plotting setup
if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
whplot = indgen(45) + 10

; run_07
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_144542_kweight.sav'
dl_07 = dl_all
diag07 = diag_nobeam

; run_08
restore, '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120716_210919_kweight.sav'
dl_08 = dl_all
diag08 = diag_nobeam

; run_07; 0809 only
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_205039_kweight_0809.sav' ; from comb_0809
dl0_0809     = dl_all
diag0_0809 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
dl_rk = dl_all
err_rk = diag_nobeam

; Plots
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

window, 2
plot, l[whplot], ((err_rk )/diag08)[whplot], thick=3, chars=1.8,charthick=2,$
  xtitle=xtitle, ytitle='err_k11/err_08', title='Error bars, K11 v.s. 2500'
oplot, [0,10000], [sqrt(2500./790.),sqrt(2500./790.)], linestyle=1, thick=3

window, 3
plot, l[whplot], ((err_rk )/diag07)[whplot], thick=3, chars=1.8,charthick=2,$
  xtitle=xtitle, ytitle='err_k11/err_07', title='Error bars, K11 v.s. 2500'
oplot, [0,10000], [sqrt(2500./790.),sqrt(2500./790.)], linestyle=1, thick=3
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v07_0617')


stop
END


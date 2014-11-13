;;;
; NAME: script_0610
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) lr jack sim
; 2) make expected_dls for all sims (still need lr)
; 3) ps plots
;
; MODIFICATION HISTORY:
;  06/10/2012: (KTS) Created
;;;

PRO end2end_azrms95
tbegin = systime(0, /seconds)
;for idx=0, 19 do begin
for idx=9, 19 do begin
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume, /use_kweight, cut_name='azrms_95'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'

END


PRO jack_sims_lr ;; THIS FAILED - need to look into this and re-run.
; jack lr
for idx=0, 19 do begin
    print, 'jack sims lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


;--------------------
; calculate expected dls for run_05
PRO read_jksims
;read_jk_sim, 'lr', '05', calib_jk=0.680
read_jk_sim, '12', '05', calib_jk=0.680
read_jk_sim, 'azrms', '05', calib_jk=0.680
read_jk_sim, 'azrms_95', '05', calib_jk=0.680
read_jk_sim, 'azrms_90', '05', calib_jk=0.680
read_jk_sim, 'tweight', '05', calib_jk=0.680
read_jk_sim, 'moon', '05', calib_jk=0.680
read_jk_sim, 'sun', '05', calib_jk=0.680
END


PRO plot_ps_run_05
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120610_190210_kweight.sav'

readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)

if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
whplot = indgen(45) + 10

;----------------
window, 0
plot,l[whplot],dl_all[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars, title='Full 2500 deg^2, run_05'
oplot,l_vec,dl_th,color=!red
errplot,l[whplot],(dl_all-diag_nobeam)[whplot]*1d12,(dl_all+diag_nobeam)[whplot]*1d12

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
;err=tvread(/png,/nodialog,filename=fdir+'ps_0610')

;----------------
; make err plot comparing run 04 to 05
diag05 = diag_nobeam / 0.680
restore, '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav'
diag04 = diag_nobeam

window, 1
plot, l[whplot], ((diag04 - diag05)/diag05)[whplot], xtitle=xtitle, ytitle='(diag04 - diag05)/diag05', $
  title='Error bars, 04 vs 05'
;err=tvread(/png,/nodialog,filename=fdir+'err_04v05_0610')

;----------------
; compare against K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_rk = diag_nobeam

window, 2
;plot, l[whplot], ((err_rk - diag05)/diag05)[whplot], xtitle=xtitle, ytitle='(diag_k11 - diag05)/diag05', $
plot, l[whplot], ((err_rk )/diag05)[whplot], xtitle=xtitle, ytitle='(diag_k11 - diag05)/diag05', $
  title='Error bars, K11 vs 05'
oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v05_0610')


stop
END

;;;
; NAME: script_0605
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) processing new sims
;
; MODIFICATION HISTORY:
;  06/05/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Processing new sims
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; coadd for jacks
; twod_tfs
; twod_kweights

; TFS
PRO make_tf_prep05 ; Done
for idx=0, 19 do begin
    print, 'make tf_prep, field ', idx
    t0 = systime(0, /seconds)

    tf_prep, idx, '05'

    print, 'That tf_prep took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO make_tfs ; Done
for idx=0, 19 do begin
    print, 'make TF, field ', idx
    t0 = systime(0, /seconds)

    make_twod_tfs, idx, '05', /dosave

    print, 'That TF took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


; kweights
PRO make_kweights ; Done
for idx=0, 19 do begin
    print, 'make KWEIGHT, field ', idx
    t0 = systime(0, /seconds)

    make_twod_kweights_lps12, idx, /save_files

    print, 'That KWEIGHT took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


; end2end
PRO end2end_all ; Done
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Jacks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;--------------------------
; Jackknives for tonight
PRO night_jacks_1
; azrms
for idx=0, 19 do begin ; (0606 - running)
    print, 'jack azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; azrms_95
for idx=0, 19 do begin
    print, 'jack azrms_95, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_95'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; azrms_90
for idx=0, 19 do begin
    print, 'jack azrms_90, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_90'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; lr
for idx=0, 19 do begin
    print, 'jack lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

END



;;; night 2
PRO night_jacks_2
; 12
for idx=0, 19 do begin
    print, 'jack 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; tweight
for idx=0, 19 do begin
    print, 'jack tweight, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'tweight'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; moon
for idx=0, 19 do begin
    print, 'jack moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; sun
idx_list = [4, 6, 11, 17, 18, 19]
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END






;;;--------------------------
; Sim Jackknives for tonight
PRO night_sim_jacks

; coadd azrms
for idx=0, 19 do begin
    print, 'coadd jack sims azrms, field ', idx
    coadd_lps12_jack_sims, idx, 'azrms'
endfor
; azrms sim jack
for idx=0, 19 do begin
    print, 'jack sims azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


; coadd azrms_95
for idx=0, 19 do begin
    print, 'coadd jack sims azrms_95, field ', idx
    coadd_lps12_jack_sims, idx, 'azrms_95'
endfor
; azrms_95 sim jack
for idx=0, 19 do begin
    print, 'jack sims azrms_95, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_95', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


; coadd azrms_90
for idx=0, 19 do begin
    print, 'coadd jack sims azrms_90, field ', idx
    coadd_lps12_jack_sims, idx, 'azrms_90'
endfor
; azrms_90 sim jack
for idx=0, 19 do begin
    print, 'jack sims azrms_90, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_90', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


; coadd moon
for idx=0, 19 do begin
    print, 'coadd jack sims moon, field ', idx
    coadd_lps12_jack_sims, idx, 'moon'
endfor
; jack moon
for idx=0, 19 do begin
    print, 'jack sims moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


; coadd sun
idx_list = [4, 6, 11, 17, 18, 19]

for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'coadd jack sims sun, field ', idx
    coadd_lps12_jack_sims, idx, 'sun'
endfor
; jack sun
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sims sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

END
;;;--------------------------





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot sv error bars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_sv_err

restore, '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav'
restore, '/home/kstory/lps12/end2end/sv_err_04_kweight.sav' ; gets sv_err

v1 = indgen(50)+8
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

window, 1, xsize=700, ysize=500
plot, l[v1], sqrt(cov_all_nobeam[v1,v1]), /ylog, xtitle='ell', ytitle='dDl [K^2]', title='Simulated SV errors'
oplot, l[v1], sqrt(cov_all_nobeam_tot[v1,v1]), color=!red                                                    
oplot, l[v1], sqrt(cov_all_nobeam_data[v1,v1]), color=!skyblue
oplot, l[v1], sqrt(cov_all_nobeam_mc[v1,v1]), color=!lavender
oplot, l[v1], sv_err[v1], color=!green
oplot, l[v1], sqrt(cov_all_nobeam[v1,v1])
legend, ['Sim_total_err', 'original_err', 'data_err', 'mc_err', 'sv_err'], linestyle=[0,0,0,0,0], $
  colors=[!white, !red, !skyblue, !lavender, !green], pos=[2000, 5.e-11]

window, 2
plot, l[v1], sqrt(cov_all_nobeam[v1,v1]) / sqrt(cov_all_nobeam_tot[v1,v1]) - 1., xtitle='ell', ytitle='(err_SVsim - err_old) / err_old', title='Fraction change in error bars'

stop
END


PRO check_tf
a
END


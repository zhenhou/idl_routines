;;;
; NAME: script_0626
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Test SV averaging
;
; MODIFICATION HISTORY:
;  06/26/2012: (KTS) Created
;;;


PRO err_5
; get KS
restore, '/home/kstory/lps12/end2end/end_ra21hdec-50_07_kweight.sav'
ks07_diag_meas = dblarr(58)
ks07_diag_sv   = dblarr(58)
ks07_diag      = dblarr(58)
ks07_dl        = spectrum
ks07_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks07_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks07_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks07_diag[i]      = sqrt(cov[i,0,i,0])
    ;ks07_diag[i]      = sqrt( ks07_diag_meas[i]^2. + ks07_diag_sv[i]^2.)
endfor
ks07_tf = transfer
ks07_ellkern = ellkern

; get RK
restore, '/home/rkeisler/ps09/1.37/end_ra21hdec-50_1.37_kweight.sav'
rk_diag_meas = dblarr(58)
rk_diag_sv   = dblarr(58)
rk_diag      = dblarr(58)
rk_dl        = spectrum
rk_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    rk_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    rk_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    rk_diag[i]      = sqrt(cov[i,0,i,0])
    ;rk_diag[i]      = sqrt( rk_diag_meas[i]^2. + rk_diag_sv[i]^2.)
endfor
rk_tf = transfer
rk_ellkern = ellkern

; get KS run_07p5 (flatsky sims)
restore, '/home/kstory/lps12/end2end/end_ra21hdec-50_07p5_flatsky_kweight.sav'
ks07p5_diag_meas = dblarr(58)
ks07p5_diag_sv   = dblarr(58)
ks07p5_diag      = dblarr(58)
ks07p5_dl        = spectrum
ks07p5_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks07p5_diag_meas[i] = sqrt(meas_cov[i,0,i,0]) / 1d12
    ks07p5_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ;ks07p5_diag[i]      = sqrt(cov[i,0,i,0])
    ks07p5_diag[i]      = sqrt( ks07p5_diag_meas[i]^2. + ks07p5_diag_sv[i]^2.)
endfor

; get KS run_07p6 (use flatsky sims w k11 filtering)
restore, '/home/kstory/lps12/end2end/end_ra21hdec-50_07p6_flatsky_k11filtering_kweight.sav'
ks07p6_diag_meas = dblarr(58)
ks07p6_diag_sv   = dblarr(58)
ks07p6_diag      = dblarr(58)
ks07p6_dl        = spectrum
ks07p6_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks07p6_diag_meas[i] = sqrt(meas_cov[i,0,i,0]) / 1d12
    ks07p6_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks07p6_diag[i]      = sqrt( ks07p6_diag_meas[i]^2. + ks07p6_diag_sv[i]^2.)
endfor

; get KS run_07p7 (use flatsky sims w k11 filtering, k11 kweight)
restore, '/home/kstory/lps12/end2end/end_ra21hdec-50_07p7_flatsky_k11filtering_k11kweight_kweight.sav'
ks07p7_diag_meas = dblarr(58)
ks07p7_diag_sv   = dblarr(58)
ks07p7_diag      = dblarr(58)
ks07p7_dl        = spectrum
ks07p7_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks07p7_diag_meas[i] = sqrt(meas_cov[i,0,i,0]) / 1d12
    ks07p7_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks07p7_diag[i]      = sqrt( ks07p7_diag_meas[i]^2. + ks07p7_diag_sv[i]^2.)
endfor

; get KS run_07p8 (use flatsky sims w k11 filtering, k11 kweight)
restore, '/home/kstory/lps12/end2end/end_ra21hdec-50_07p8_flatsky_k11filtering_nokweight.sav'
ks07p8_diag_meas = dblarr(58)
ks07p8_diag_sv   = dblarr(58)
ks07p8_diag      = dblarr(58)
ks07p8_dl        = spectrum
ks07p8_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks07p8_diag_meas[i] = sqrt(meas_cov[i,0,i,0]) / 1d12
    ks07p8_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks07p8_diag[i]      = sqrt( ks07p8_diag_meas[i]^2. + ks07p8_diag_sv[i]^2.)
endfor

; get KS run_07p9 (use flatsky sims w k11 filtering, k11 kweight)
restore, '/home/kstory/lps12/end2end/end_ra21hdec-50_07p9_RKsims_kweight.sav'
ks07p9_diag_meas = dblarr(58)
ks07p9_diag_sv   = dblarr(58)
ks07p9_diag      = dblarr(58)
ks07p9_dl        = spectrum
ks07p9_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks07p9_diag_meas[i] = sqrt(meas_cov[i,0,i,0]) / 1d12
    ks07p9_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks07p9_diag[i]      = sqrt( ks07p9_diag_meas[i]^2. + ks07p9_diag_sv[i]^2.)
endfor

vec = indgen(27)+9
mm07 = mean((rk_diag_sv/ks07_diag_sv)[vec])
mm07p5 = mean((rk_diag_sv/ks07p5_diag_sv)[vec])
mm07p6 = mean((rk_diag_sv/ks07p6_diag_sv)[vec])
mm07p7 = mean((rk_diag_sv/ks07p7_diag_sv)[vec])
mm07p8 = mean((rk_diag_sv/ks07p8_diag_sv)[vec])
mm07p9 = mean((rk_diag_sv/ks07p9_diag_sv)[vec])


; Get 0809 numbers
; run_07; 0809 only
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120619_195218_kweight_0809.sav'
dl_0809     = dl_all
diag_0809 = diag_nobeam
dd_e2e = dblarr(58)
for i=0, 57 do dd_e2e[i] = sqrt(cov_all_nobeam_e2e[i,i])

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_rk = diag_nobeam



;---------------------------
; PLOTS
;---------------------------
whplot = indgen(45) + 10
fdir = '/home/kstory/public_html/notebook/spt_lps12/'
l = banddef-25.

;;;;;;;;;;;;;;;;
; combined plot
window, 10, xsize=800,ysize=500
plot, l[vec], (rk_diag_sv/ks07_diag_sv)[vec], yr=[.75,1.1],/yst,title='SV Errors, ra21hdec-50', $ ; tot
  ytitle='err_rk / err_07pX', xtitle='ell'
oplot, l[vec], (rk_diag_sv/ks07_diag_sv)[vec]
oplot, l[vec], (rk_diag_sv/ks07p5_diag_sv)[vec], color=!green
oplot, l[vec], (rk_diag_sv/ks07p6_diag_sv)[vec], color=!red
oplot, l[vec], (rk_diag_sv/ks07p7_diag_sv)[vec], color=!skyblue
oplot, l[vec], (rk_diag_sv/ks07p8_diag_sv)[vec], color=!purple
;oplot, l[vec], (rk_diag_sv/ks07p9_diag_sv)[vec], color=!orange
oplot, l[vec], (err_rk/diag_0809)[vec], color=!yellow, linestyle=1, thick=2

oplot, [0,10000], [mm07,mm07]
oplot, [0,10000], [mm07p5,mm07p5], color=!green
oplot, [0,10000], [mm07p6,mm07p6], color=!red
oplot, [0,10000], [mm07p7,mm07p7], color=!skyblue
oplot, [0,10000], [mm07p8,mm07p8], color=!purple
legend, ['07', '07p5', '07p6', '07p7', '07p8', '0809'], color=[!white, !green, !red, !skyblue, !purple, !yellow], linestyle=[0,0,0,0,0,1], pos=[550, 0.87]
;oplot, [0,10000], [mm07p9,mm07p9], color=!orange
;legend, ['07', '07p5', '07p6', '07p7', '07p8', '07p9'], color=[!white, !green, !red, !skyblue, !purple, !orange], linestyle=0, pos=[550, 0.87]
err=tvread(/png,/nodialog,filename=fdir+'err_07pX_0627a')



; window, 1
; plot, l[whplot], (rk_diag/ks07_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, ks07, ra21hdec-50', $ ; tot
;   ytitle='err_rk / err_07', xtitle='ell'
; ;oplot, l[whplot], (rk_diag_meas/ks07_diag_meas)[whplot], color=!skyblue         ; meas
; oplot, l[whplot], (rk_diag_sv/ks07_diag_sv)[whplot], color=!green               ; sv
; legend, ['Tot', 'SV', 'NV'], color=[!white, !green, !skyblue], linestyle=[0,0,0], pos=[600, 1.2]
; ;err=tvread(/png,/nodialog,filename=fdir+'err_rk_vs_07_0626')

; window, 5
; plot, l[whplot], (rk_diag/ks07p5_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, ks07p5, ra21hdec-50', $ ; tot
;   ytitle='err_rk / err_07p5', xtitle='ell'
; ;oplot, l[whplot], (rk_diag_meas/ks07p5_diag_meas)[whplot], color=!skyblue         ; meas
; oplot, l[whplot], (rk_diag_sv/ks07p5_diag_sv)[whplot], color=!green               ; sv
; legend, ['Tot', 'SV', 'NV'], color=[!white, !green, !skyblue], linestyle=[0,0,0], pos=[600, 1.2]
; ;err=tvread(/png,/nodialog,filename=fdir+'err_rk_vs_07p5_0626')

; window, 6
; plot, l[whplot], (rk_diag/ks07p6_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, ks07p6, ra21hdec-50', $ ; tot
;   ytitle='err_rk / err_07p6', xtitle='ell'
; ;oplot, l[whplot], (rk_diag_meas/ks07p6_diag_meas)[whplot], color=!skyblue         ; meas
; oplot, l[whplot], (rk_diag_sv/ks07p6_diag_sv)[whplot], color=!green               ; sv
; legend, ['Tot', 'SV', 'NV'], color=[!white, !green, !skyblue], linestyle=[0,0,0], pos=[600, 1.2]
; ;err=tvread(/png,/nodialog,filename=fdir+'err_rk_vs_07p6_0626')

; window, 7
; plot, l[whplot], (rk_diag/ks07p7_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, ks07p7, ra21hdec-50', $ ; tot
;   ytitle='err_rk / err_07p7', xtitle='ell'
; ;oplot, l[whplot], (rk_diag_meas/ks07p7_diag_meas)[whplot], color=!skyblue         ; meas
; oplot, l[whplot], (rk_diag_sv/ks07p7_diag_sv)[whplot], color=!green               ; sv
; legend, ['Tot', 'SV', 'NV'], color=[!white, !green, !skyblue], linestyle=[0,0,0], pos=[600, 1.2]
; ;err=tvread(/png,/nodialog,filename=fdir+'err_rk_vs_07p7_0626')

; window, 8
; plot, l[whplot], (rk_diag/ks07p8_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, ks07p8, ra21hdec-50', $ ; tot
;   ytitle='err_rk / err_07p8', xtitle='ell'
; ;oplot, l[whplot], (rk_diag_meas/ks07p8_diag_meas)[whplot], color=!skyblue         ; meas
; oplot, l[whplot], (rk_diag_sv/ks07p8_diag_sv)[whplot], color=!green               ; sv
; legend, ['Tot', 'SV', 'NV'], color=[!white, !green, !skyblue], linestyle=[0,0,0], pos=[600, 1.2]
; ;err=tvread(/png,/nodialog,filename=fdir+'err_rk_vs_07p8_0626')

stop
END


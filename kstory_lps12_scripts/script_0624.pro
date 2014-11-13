;;;
; NAME: script_0624
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Test SV averaging
;
; MODIFICATION HISTORY:
;  06/22/2012: (KTS) Created
;;;


PRO err_5
; get KS
restore, '/home/kstory/lps12/end2end/end_ra21hdec-50_07_kweight.sav'
ks_diag_meas = dblarr(58)
ks_diag_sv   = dblarr(58)
ks_diag      = dblarr(58)
ks_dl        = spectrum
ks_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks_diag[i]      = sqrt(cov[i,0,i,0])
    ;ks_diag[i]      = sqrt( ks_diag_meas[i]^2. + ks_diag_sv[i]^2.)
endfor
ks_tf = transfer
ks_ellkern = ellkern

; get KS run_07p5 (use end2end_lowell)
restore, '/home/kstory/lps12/end2end/end_ra21hdec-50_07p5_kweight.sav'
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



;---------------------------
; PLOTS
;---------------------------
whplot = indgen(45) + 10
fdir = '/home/kstory/public_html/notebook/spt_lps12/'
l = banddef-25.

window, 1
plot, l[whplot], (ks07p5_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, 7p5, ra21hdec-50', $ ; tot
  ytitle='err_07p5 / err_07', xtitle='ell'
oplot, l[whplot], (ks07p5_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (ks07p5_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
legend, ['Tot', 'SV', 'NV'], color=[!white, !green, !skyblue], linestyle=[0,0,0], pos=[600, 1.2]
;err=tvread(/png,/nodialog,filename=fdir+'err_07p5_vs_07_0624')

window, 2
plot, l[whplot], (rk_diag/ks07p5_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, 7p5, ra21hdec-50', $ ; tot
  ytitle='err_rk / err_07p5', xtitle='ell'
oplot, l[whplot], (rk_diag_meas/ks07p5_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (rk_diag_sv/ks07p5_diag_sv)[whplot], color=!green               ; sv
legend, ['Tot', 'SV', 'NV'], color=[!white, !green, !skyblue], linestyle=[0,0,0], pos=[600, 1.2]
err=tvread(/png,/nodialog,filename=fdir+'err_rk_vs_07p5_0624')

window, 3
plot, l[whplot], (rk_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, 7p5, ra21hdec-50', $ ; tot
  ytitle='err_rk / err_07', xtitle='ell'
oplot, l[whplot], (rk_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (rk_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
legend, ['Tot', 'SV', 'NV'], color=[!white, !green, !skyblue], linestyle=[0,0,0], pos=[600, 1.2]
;err=tvread(/png,/nodialog,filename=fdir+'err_rk_vs_07_0624')

; window, 4
; plot, l[whplot], (ks07p3_diag/ks07p5_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, 7p3, ra21hdec-50', $ ; tot
;   ytitle='err_2160 / err_4320', xtitle='ell'
; oplot, l[whplot], (ks07p3_diag_meas/ks07p5_diag_meas)[whplot], color=!skyblue         ; meas
; oplot, l[whplot], (ks07p3_diag_sv/ks07p5_diag_sv)[whplot], color=!green               ; sv
; ;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0624')

stop
END

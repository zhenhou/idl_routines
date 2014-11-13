;;;
; NAME: script_0622
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Test SV averaging
;
; MODIFICATION HISTORY:
;  06/22/2012: (KTS) Created
;;;




;;;;;;;;;;;;;;;;;;;;;;
; Test SV averaging
;;;;;;;;;;;;;;;;;;;;;;

PRO test_sv

nscan=100
nsamp=10000
seed=3
tod = randomn(seed,nscan,nsamp)




stop
END




;;;;;;;;;;;;;;;;;;;;;;
; Compare sv error in different fields
;;;;;;;;;;;;;;;;;;;;;;

PRO sv_err, idx, stopit=stopit
f = lps12_fieldstruct()

frac=[0.75, 0.8743, 0, 0.9191, 0.8670, 0.8829]

; KS
restore, '/home/kstory/lps12/end2end/end_'+f[idx].name+'_07_kweight.sav'
ks_diag_meas = dblarr(58)
ks_diag_sv   = dblarr(58)
ks_diag      = dblarr(58)
ks_dl        = spectrum
for i=0, 57 do begin
    ks_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks_diag[i]      = sqrt(cov[i,0,i,0])
    ;ks_diag[i]      = sqrt( ks_diag_meas[i]^2. + ks_diag_sv[i]^2.)
endfor
ks_area = total(mask_padded / (60.*60.)) ; in deg^2

; get RK
restore, '/home/rkeisler/ps09/1.37/end_'+f[idx].dir_name+'_1.37_kweight.sav'
rk_diag_meas = dblarr(58)
rk_diag_sv   = dblarr(58)
rk_diag      = dblarr(58)
rk_dl        = spectrum
for i=0, 57 do begin
    rk_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    rk_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    rk_diag[i]      = sqrt(cov[i,0,i,0])
endfor
rk_area = total(mask / (60.*60.)) ; in deg^2

;---------------------------
; PLOTS
;---------------------------
whplot = indgen(45) + 10
fdir = '/home/kstory/public_html/notebook/spt_lps12/'
l = banddef-25.

window, idx
plot, l[whplot], (rk_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, '+f[idx].name, $ ; tot
  ytitle='err_rk / err_ks', xtitle='ell'
oplot, l[whplot], (rk_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (rk_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0622')

print, '*** Field '+f[idx].name
print, 'Mean, Median SV                     = ', mean((rk_diag_sv/ks_diag_sv)[whplot]), median((rk_diag_sv/ks_diag_sv)[whplot])
print, 'sqrt(ks_area * frac[idx] / rk_area) = ', sqrt(ks_area * frac[idx] / rk_area)
print, 'sqrt(cut pix frac)                  = ', sqrt(frac[idx])
print, 'sqrt(ks_area / rk_area)             = ', sqrt(ks_area / rk_area)



if keyword_set(stopit) then stop
END

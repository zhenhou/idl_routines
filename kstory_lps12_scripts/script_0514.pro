;;;
; NAME: script_0514
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  05/14/2012: (KTS) Created
;;;

;;; Processing jackknives

PRO coadd_jack_sims_lr ; done
for idx=0, 19 do begin
    print, 'coadd jack sims lr, field ', idx
    coadd_lps12_jack_sims, idx, 'lr'
endfor
END

PRO coadd_jack_sims_tweight ; done
for idx=0, 19 do begin
    print, 'coadd jack sims tweight, field ', idx
    coadd_lps12_jack_sims, idx, 'tweight'
endfor
END

PRO coadd_jack_sims_azrms ; done
for idx=0, 19 do begin
    print, 'coadd jack sims azrms, field ', idx
    coadd_lps12_jack_sims, idx, 'azrms'
endfor
END


;;; jack sims
PRO jack_sims_lr ; Done
for idx=0, 19 do begin
    if idx eq 6 then idx = 7 ; already done with 6
    print, 'jack sims lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


PRO jack_sims_12 ; Done
for idx=0, 19 do begin
    print, 'jack sims 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_tweight ; running
for idx=0, 19 do begin
    print, 'jack sims tweight, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'tweight', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

;;;--------------------------
; Jackknives for tonight
PRO night_sim_jacks
; azrms sim jack
for idx=0, 19 do begin
    print, 'jack sims azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms', /sim
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


;;; jacks
PRO jack_azrms
for idx=0, 19 do begin
    print, 'jack azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sun
idx_list = [4, 6, 11, 17, 18, 19]

for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_moon
for idx=0, 19 do begin
    print, 'jack moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


;;; Make 'expected_jack_dls.sav'
PRO make_expected_jack_dls

jack_list = ['lr', '12', 'azrms', 'sun', 'moon']

jackdir = '/home/kstory/lps12/jacks/sims/'

i0=1
i1=5
banddef = (1+findgen(6))*500.
nbin=i1-i0+1

expected_dls=fltarr(nbin,20)
save, expected_dls, filename = jackdir+'expected_jack_dls_lr.sav'
save, expected_dls, filename = jackdir+'expected_jack_dls_12.sav'
save, expected_dls, filename = jackdir+'expected_jack_dls_azrms.sav'
save, expected_dls, filename = jackdir+'expected_jack_dls_tweight.sav'
save, expected_dls, filename = jackdir+'expected_jack_dls_sun.sav'
save, expected_dls, filename = jackdir+'expected_jack_dls_moon.sav'

END

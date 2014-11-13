;;;
; NAME: script_0513
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  05/13/2012: (KTS) Created
;;;

PRO coadd_jack_sims_12
for idx=0, 19 do begin
    print, 'coadd jack sims 12, field ', idx
    coadd_lps12_jack_sims, idx, '12'
endfor
END

PRO coadd_jack_sims_lr
for idx=0, 19 do begin
    print, 'coadd jack sims lr, field ', idx
    coadd_lps12_jack_sims, idx, 'lr'
endfor
END


PRO jack_sims_12
for idx=0, 19 do begin
    print, 'jack sims 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

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
for idx=0, 19 do begin
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


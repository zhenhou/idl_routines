;;;
; NAME: script_0804
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) jackknives, run_08
;;;


;;;;;;;;;;;;;;;;;;;;;;
; Run jobs
;;;;;;;;;;;;;;;;;;;;;;
PRO jacks1
;jack_lr
;jack_12
jack_moon
END

PRO jacks2
jack_azrms
jack_tweight
jack_sun
END

PRO jacks3
jack_sims_lr
jack_sims_12
jack_sims_moon
END

PRO jacks4
jack_sims_azrms
jack_sims_tweight
jack_sims_sun
END


;;;;;;;;;;;;;;;;;;;;;;
; Data jacks
;;;;;;;;;;;;;;;;;;;;;;

PRO jack_lr
for idx=0, 19 do begin
    print, 'jack lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr', run='08'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_12
for idx=0, 19 do begin
    print, 'jack 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12', run='08'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_tweight
for idx=0, 19 do begin
    print, 'jack tweight, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'tweight', run='08'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_moon
for idx=0, 19 do begin
    print, 'jack moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon', run='08'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sun
idx_list = [4, 6, 11, 17, 18, 19]
for ii=0, 5 do begin
    idx = idx_list[ii]
    print, 'jack sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun', run='08'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_azrms
;for idx=0, 19 do begin
for idx=1, 19 do begin
    print, 'jack azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms', run='08'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


;;;;;;;;;;;;;;;;;;;;;;
; Sim jacks
;;;;;;;;;;;;;;;;;;;;;;

PRO jack_sims_lr
for idx=0, 19 do begin
     print, 'coadd jack sims lr, field ', idx
    coadd_lps12_jack_sims, idx, 'lr', run='08'
endfor
for idx=0, 19 do begin
    print, 'jack sims lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr', run='08', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_12
for idx=0, 19 do begin
     print, 'coadd jack sims 12, field ', idx
    coadd_lps12_jack_sims, idx, '12', run='08'
endfor
for idx=0, 19 do begin
    print, 'jack sims 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12', run='08', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_tweight
for idx=0, 19 do begin
     print, 'coadd jack sims tweight, field ', idx
    coadd_lps12_jack_sims, idx, 'tweight', run='08'
endfor
for idx=0, 19 do begin
    print, 'jack sims tweight, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'tweight', run='08', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_moon
for idx=0, 19 do begin
     print, 'coadd jack sims moon, field ', idx
    coadd_lps12_jack_sims, idx, 'moon', run='08'
endfor
for idx=0, 19 do begin
    print, 'jack sims moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon', run='08', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_sun
idx_list = [4, 6, 11, 17, 18, 19]
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
     print, 'coadd jack sims sun, field ', idx
    coadd_lps12_jack_sims, idx, 'sun', run='08'
endfor
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sims sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun', run='08', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_azrms
for idx=0, 19 do begin
     print, 'coadd jack sims azrms, field ', idx
    coadd_lps12_jack_sims, idx, 'azrms', run='08'
endfor
for idx=0, 19 do begin
    print, 'jack sims azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms', run='08', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


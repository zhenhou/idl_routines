;;;
; NAME: script_0609
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Processing sim jacks
;
; MODIFICATION HISTORY:
;  06/09/2012: (KTS) Created
;;;

PRO jack_sims_moon
; coadd moon
for idx=0, 19 do begin
    print, 'coadd jack sims moon, field ', idx
    coadd_lps12_jack_sims, idx, 'moon', out_dir='/data18/kstory/lps12/sims/'
endfor
; jack moon
for idx=0, 19 do begin
    print, 'jack sims moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon', /sim, dir_sim='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


PRO jack_sims_sun
; coadd sun
idx_list = [4, 6, 11, 17, 18, 19]

for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'coadd jack sims sun, field ', idx
    coadd_lps12_jack_sims, idx, 'sun', out_dir='/data18/kstory/lps12/sims/'
endfor
; jack sun
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sims sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun', /sim, dir_sim='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

END


PRO jack_sims_lr ;; THIS FAILED - need to look into this and re-run.
; jack lr
for idx=0, 19 do begin
    print, 'jack sims lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr', /sim, dir_sim='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


;;;--------------------------


PRO jack_sims_12andtweight
; coadd 12
for idx=0, 19 do begin
    print, 'coadd jack sims 12, field ', idx
    coadd_lps12_jack_sims, idx, '12', out_dir='/data18/kstory/lps12/sims/'
endfor
; jack 12
for idx=0, 19 do begin
    print, 'jack sims 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12', /sim, dir_sim='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


; coadd tweight
for idx=0, 19 do begin
    print, 'coadd jack sims tweight, field ', idx
    coadd_lps12_jack_sims, idx, 'tweight', out_dir='/data18/kstory/lps12/sims/'
endfor
; jack tweight
for idx=0, 19 do begin
    print, 'jack sims tweight, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'tweight', /sim, dir_sim='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END



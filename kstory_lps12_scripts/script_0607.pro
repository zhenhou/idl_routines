;;;
; NAME: script_0607
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Processing
;
; MODIFICATION HISTORY:
;  06/07/2012: (KTS) Created
;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Re-run end2end run_05
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO end2end_run_05_0606 ; Leave this to run fields 8, 9, 10
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    if idx eq 0 then idx++
    if idx eq 1 then idx++
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume, /use_kweight, out_dir='/data18/kstory/lps12/end2end/'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END

PRO end2end_9to10 ; Job 1, part b
tbegin = systime(0, /seconds)
for idx=9, 10 do begin
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume, /use_kweight, out_dir='/data18/kstory/lps12/end2end/'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END

; Job 2
PRO end2end_11to13
tbegin = systime(0, /seconds)
for idx=11, 13 do begin
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume, /use_kweight, out_dir='/data18/kstory/lps12/end2end/'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END

; Job 3
PRO end2end_14to16
tbegin = systime(0, /seconds)
for idx=14, 16 do begin
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume, /use_kweight, out_dir='/data18/kstory/lps12/end2end/'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END

; Job 4
PRO end2end_17to19
tbegin = systime(0, /seconds)
for idx=17, 19 do begin
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume, /use_kweight, out_dir='/data18/kstory/lps12/end2end/'

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
; for idx=0, 19 do begin
;     print, 'jack azrms, field ', idx
;     t0 = systime(0, /seconds)
;     lps12_jack, idx, 'azrms'
;     print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
; endfor

; azrms_95
for idx=0, 19 do begin
    print, 'jack azrms_95, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_95', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; azrms_90
for idx=0, 19 do begin
    print, 'jack azrms_90, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_90', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; lr
for idx=0, 19 do begin
    print, 'jack lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

END



;;; night 2
PRO night_jacks_2
; 12
for idx=18, 19 do begin ; only last 2 files left to run
    print, 'jack 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; tweight
for idx=0, 19 do begin
    print, 'jack tweight, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'tweight', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; moon
for idx=0, 19 do begin
    print, 'jack moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; sun
idx_list = [4, 6, 11, 17, 18, 19]
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

;;;--------------------------
; Sim Jackknives for tonight
PRO night_sim_jacks

; ; coadd azrms
; for idx=0, 19 do begin
;     print, 'coadd jack sims azrms, field ', idx
;     coadd_lps12_jack_sims, idx, 'azrms'
; endfor
; ; azrms sim jack
; for idx=0, 19 do begin
;     print, 'jack sims azrms, field ', idx
;     t0 = systime(0, /seconds)
;     lps12_jack, idx, 'azrms', /sim
;     print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
; endfor


; ; coadd azrms_95
; for idx=0, 19 do begin
;     print, 'coadd jack sims azrms_95, field ', idx
;     coadd_lps12_jack_sims, idx, 'azrms_95', out_dir='/data18/kstory/lps12/sims/'
; endfor
; azrms_95 sim jack
for idx=0, 19 do begin
    print, 'jack sims azrms_95, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_95', /sim, dir_sim='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


; coadd azrms_90
for idx=0, 19 do begin
    print, 'coadd jack sims azrms_90, field ', idx
    coadd_lps12_jack_sims, idx, 'azrms_90', out_dir='/data18/kstory/lps12/sims/'
endfor
; azrms_90 sim jack
for idx=0, 19 do begin
    print, 'jack sims azrms_90, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_90', /sim, dir_sim='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


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


; ; coadd sun
; idx_list = [4, 6, 11, 17, 18, 19]

; for ii=0, n_elements(idx_list)-1 do begin
;     idx = idx_list[ii]
;     print, 'coadd jack sims sun, field ', idx
;     coadd_lps12_jack_sims, idx, 'sun', out_dir='/data18/kstory/lps12/sims/'
; endfor
; ; jack sun
; for ii=0, n_elements(idx_list)-1 do begin
;     idx = idx_list[ii]
;     print, 'jack sims sun, field ', idx
;     t0 = systime(0, /seconds)
;     lps12_jack, idx, 'sun', /sim, dir_sim='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
;     print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
; endfor

END
;;;--------------------------


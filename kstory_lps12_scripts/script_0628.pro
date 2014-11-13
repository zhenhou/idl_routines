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

;;;;;;;;;;;;;;;;;;;;;;
; End2End, run_08
;;;;;;;;;;;;;;;;;;;;;;

PRO end2end_run_08_a
tbegin = systime(0, /seconds)
for idx=0, 5 do begin
    print, 'run end2end_08, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_08, idx, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END

PRO end2end_run_08_b
tbegin = systime(0, /seconds)
for idx=6, 12 do begin
    print, 'run end2end_08, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_08, idx, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END

PRO end2end_run_08_c
tbegin = systime(0, /seconds)
for idx=13, 19 do begin
    print, 'run end2end_08, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_08, idx, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END


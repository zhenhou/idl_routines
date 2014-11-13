;;;
; NAME: script_0731
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Run pipetest2, run_08p2
; 2) SV-NV test; Run try_end2end_08, lr option
;;;

;;; 1) Run pipetest2, run_08p2
PRO end2end_run08p2
idx = 4
tbegin = systime(0, /seconds)
print, 'run end2end_08b1, field ', idx
t0 = systime(0, /seconds)

try_end2end_08p2, idx, /use_kweight, /resume

tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'

END




;;;;;;;;;;;;;;;;;;;;;;;
; 2) SV-NV test
;;;;;;;;;;;;;;;;;;;;;;;

; Run 0 to 6
PRO end2end_run08p2_0to6
tbegin = systime(0, /seconds)
t0 = systime(0, /seconds)
for idx=0, 6 do begin
    ti1 = systime(0, /seconds)
    try_end2end_08, idx, /use_kweight, /resume, /lr
    ti2 = systime(0, /seconds)
    print, 'Field ', idx, ' took ', (ti2-ti1)/60., 'minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END


; Run 7 to 13
PRO end2end_run08p2_7to13
tbegin = systime(0, /seconds)
t0 = systime(0, /seconds)
for idx=7, 13 do begin
    ti1 = systime(0, /seconds)
    try_end2end_08, idx, /use_kweight, /resume, /lr
    ti2 = systime(0, /seconds)
    print, 'Field ', idx, ' took ', (ti2-ti1)/60., 'minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END


; Run 14 to 19
PRO end2end_run08p2_14to19
tbegin = systime(0, /seconds)
t0 = systime(0, /seconds)
for idx=14, 19 do begin
    ti1 = systime(0, /seconds)
    try_end2end_08, idx, /use_kweight, /resume, /lr
    ti2 = systime(0, /seconds)
    print, 'Field ', idx, ' took ', (ti2-ti1)/60., 'minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END



;;;
; NAME: script_0615
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Make azrms jackgoodfiles, lists.
; 2) processing; end2end, run_07
; 3) sim jacks
;
; MODIFICATION HISTORY:
;  06/15/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;
; Make azrms jackgoodfiles
;;;;;;;;;;;;;;;;;;;;;;

PRO make_azrms_goodfiles
f = lps12_fieldstruct()

for idx=0, 19 do begin
    print, 'Make lists for field ', idx, ': ' + f[idx].name
    make_azrms_jack_list, idx
    make_azrms_jack_list, idx, /preaz
endfor
END




;;;;;;;;;;;;;;;;;;;;;;
; End2End, run_07
;;;;;;;;;;;;;;;;;;;;;;

PRO end2end_run_07
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    print, 'run end2end_07, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_07, idx, /resume, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END



;;;;;;;;;;;;;;;;;;;;;;
; Jack_sims
;;;;;;;;;;;;;;;;;;;;;;

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

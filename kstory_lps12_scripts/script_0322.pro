;;;
; NAME: script_0322
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) all of this work was re-done on 0323
;
; MODIFICATION HISTORY:
;  03/22/2012: (KTS) Created
;;;

;...................................................................
; make tweight jackknife runlists
PRO run_make_tweight_jack_list
compile_opt IDL2, HIDDEN

;for ii=0, nfields-1 do begin
for ii=1, 19 do begin

    make_tweight_jack_list, ii

endfor
END

;...................................................................
; run tweight jackknife
PRO run_tweight_jack
compile_opt IDL2, HIDDEN

;for ii=0, nfields-1 do begin
for ii=1, 19 do begin

    lps12_jack, ii, 'tweight', savdir='/home/kstory/lps12/jacks/'

endfor
END

;...................................................................
; run moon jackknife
PRO run_moon_jack
compile_opt IDL2, HIDDEN

;for ii=0, nfields-1 do begin
for ii=1, 19 do begin

    lps12_jack, ii, 'moon', savdir='/home/kstory/lps12/jacks/'

endfor
END

;...................................................................
; run moon jackknife
PRO run_sun_jack
compile_opt IDL2, HIDDEN

restore, '/data/kstory/projects/lps12/masks/masks_50mJy/sun_weights.sav'
idx = where(p_used gt 0.25)
nfields = n_elements(idx)

for ii=0, nfields-1 do begin

    field_idx = idx[ii]
    lps12_jack, field_idx, 'sun', savdir='/home/kstory/lps12/jacks/'

endfor
END


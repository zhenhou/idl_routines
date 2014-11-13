;;;
; NAME: script_0329
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  03/26/2012: (KTS) Created
;;;

;...................................................................
; run all jacks
PRO run_12_jack
compile_opt IDL2, HIDDEN
;for ii=0, 19 do lps12_jack, ii, '12', savdir='/home/kstory/lps12/jacks/'
for ii=7, 19 do lps12_jack, ii, '12', savdir='/home/kstory/lps12/jacks/', $
  mapdir='/scratch/kstory/20120305/'
END

PRO run_azrms_jack
compile_opt IDL2, HIDDEN
for ii=14, 19 do lps12_jack, ii, 'azrms', savdir='/home/kstory/lps12/jacks/', $
  mapdir='/scratch/kstory/20120305/'
END

PRO run_tweight_jack
compile_opt IDL2, HIDDEN
for ii=0, 19 do lps12_jack, ii, 'tweight', savdir='/home/kstory/lps12/jacks/', $
  mapdir='/scratch/kstory/20120305/'
END

PRO run_moon_jack
compile_opt IDL2, HIDDEN
for ii=0, 19 do lps12_jack, ii, 'moon', savdir='/home/kstory/lps12/jacks/', $
  mapdir='/scratch/kstory/20120305/'
END

PRO run_sun_jack
compile_opt IDL2, HIDDEN
list = [4,6,11,17,18,19]
for idx=0, n_elements(list)-1 do begin
    ii = list[idx]
    lps12_jack, ii, 'sun', savdir='/home/kstory/lps12/jacks/', $
      mapdir='/scratch/kstory/20120305/'
endfor
END


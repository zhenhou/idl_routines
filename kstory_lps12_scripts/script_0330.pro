;;;
; NAME: script_0330
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
; run all jacks for ra3h30dec-42.5 field, main node
PRO run_17_jack
compile_opt IDL2, HIDDEN
lps12_jack, 17, 'lr', savdir='/home/kstory/lps12/jacks/'
lps12_jack, 17, '12', savdir='/home/kstory/lps12/jacks/'
lps12_jack, 17, 'azrms', savdir='/home/kstory/lps12/jacks/'
lps12_jack, 17, 'tweight', savdir='/home/kstory/lps12/jacks/'
lps12_jack, 17, 'sun', savdir='/home/kstory/lps12/jacks/'
lps12_jack, 17, 'moon', savdir='/home/kstory/lps12/jacks/'
END



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
;for ii=0, 19 do lps12_jack, ii, 'tweight', savdir='/home/kstory/lps12/jacks/', $
for ii=17, 17 do lps12_jack, ii, 'tweight', savdir='/home/kstory/lps12/jacks/', $
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


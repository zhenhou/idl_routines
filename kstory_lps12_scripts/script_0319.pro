;;;
; NAME: script_0319
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) sss, Study L-R leakage into 2010, 2011 maps
;
; MODIFICATION HISTORY:
;  03/19/2012: (KTS) Created
;;;


;...................................................................
; run lr jack on all fields
PRO run_mk_az_quality_0319
compile_opt IDL2, HIDDEN
;for ii=0, 19 do begin
for ii=2, 19 do begin
    mk_az_quality, ii
endfor
END


;...................................................................
; run lr jack on all fields
PRO run_lr_jack_0319
compile_opt IDL2, HIDDEN
for ii=0, 19 do begin
    lps12_jack, ii, 'lr', savdir = '/home/kstory/lps12/jacks/0319/'
endfor
END


;...................................................................
; check azrms
PRO sss, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()

fst = field_arr_[field_idx]
field_name = fst.name

print, 'Analyzing field ', field_idx, ', name = ', field_name

savdir  = '/data/kstory/projects/lps12/jacks/az_quality/'
savname = savdir + 'ground_rms_'+field_name+'.sav'
restore, '/data/kstory/projects/lps12/jacks/az_quality/goodfiles_'+field_name+'_az_lowrms.txt'
restore, '/data/kstory/projects/lps12/jacks/az_quality/goodfiles_'+field_name+'_az_highrms.txt'

END

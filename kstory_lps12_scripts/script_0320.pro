;;;
; NAME: script_0320
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) sss, Study L-R leakage into 2010, 2011 maps
;
; MODIFICATION HISTORY:
;  03/20/2012: (KTS) Created
;;;


;...................................................................
; run lr jack on all fields
PRO run_mk_az_quality_0320
compile_opt IDL2, HIDDEN
;for ii=0, 19 do begin
for ii=4, 19 do begin
    mk_az_quality, ii
endfor
END


;...................................................................
; run lr jack on all fields
PRO run_lr_jack_2to10
compile_opt IDL2, HIDDEN
for ii=2, 10 do begin
    lps12_jack, ii, 'lr', savdir = '/home/kstory/lps12/jacks/0319/'
endfor
END

PRO run_lr_jack_11to19
compile_opt IDL2, HIDDEN
for ii=11, 19 do begin
    lps12_jack, ii, 'lr', savdir = '/home/kstory/lps12/jacks/0319/'
endfor
END


;...................................................................
; run azrms jack on all fields
PRO run_azrms_jack
compile_opt IDL2, HIDDEN
for ii=0, 19 do begin
    lps12_jack, ii, 'azrms', savdir = '/home/kstory/lps12/jacks/0319/'
endfor
END


;...................................................................
; Script
PRO sss, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()

fst = field_arr_[field_idx]
field_name = fst.name

print, 'Analyzing field ', field_idx, ', name = ', field_name
END


pro ss1
field = 'ra23h30dec-55_2010'
jfile  = '/data/kstory/projects/lps12/jacks/az_quality/goodfiles_'+field+'_az_lowrms.txt'
jfile2 = '/data/kstory/projects/lps12/jacks/az_quality/goodfiles_'+field+'_az_highrms.txt'
readcol,jfile,g1,format='a'
readcol,jfile2,g2,format='a'
stop
end

;;;
; NAME: script_0208
; PURPOSE:
;   Make IDL sav files for 1hz bolos
;
; NOTES:
;
; MODIFICATION HISTORY:
;  02/08/2012: (KTS) Created
;;;

; Note, no fits files for 2 yet
pro get_1hz_0208_0to5
compile_opt IDL2, HIDDEN
field_idx_list = [0,1,3,4,5]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** get_1hz_bolo_info for field, ',field_arr_[field_idx_list[ii]].name
    get_1hz_bolo_info, ii
endfor
end


; No fits files for anything else... :(

;;;
; NAME: script_0901
; PURPOSE:
;   
;
; NOTES:
;
; MODIFICATION HISTORY:
;  09/01/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro script_0901
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'
;exclude_list = [1,2,5,11]
;field_idx_list = setdifference(lindgen(n_elements(all_fields)), exclude_list)
field_idx_list = lindgen(n_elements(all_fields))

for ii=0, n_elements(field_idx_list)-1 do begin
    idx = field_idx_list[ii]
    field_name = all_fields[idx]
    print, "************* make_final_runlists, "+field_name
    make_final_runlists, field_name
endfor
stop
end


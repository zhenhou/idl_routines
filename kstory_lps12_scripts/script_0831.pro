;;;
; NAME: script_0831
; PURPOSE:
;   Make pdf maps from final runlists
;
; NOTES:
;
; MODIFICATION HISTORY:
;  08/29/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro script_0831
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'
;field_idx_list = [0,1,2,5,11]
field_name = 'ra23h30dec-55'

plot_raw_maps, field_name, dates

end


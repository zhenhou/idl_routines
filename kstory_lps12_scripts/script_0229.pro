;;;
; NAME: script_0227
; PURPOSE:
;   Make final runlists for ra5h30dec-55_2011
;
; NOTES:
;   calls make_runlist_for_lps12_single_field from make_runlist_for_lps12.pro
;
; MODIFICATION HISTORY:
;  02/27/2012: (KTS) Created
;;;

;...................................................................
; Make runlists
;
pro make_runlist_0229
compile_opt IDL2, HIDDEN

field_idx = 20
field_arr_ = lps12_fieldstruct()
print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx].name
make_runlist_for_lps12_single_field, field_idx

end

;...................................................................
; Make plot summaries
;
pro plot_0229
compile_opt IDL2, HIDDEN

idx = 20
field_arr_ = lps12_fieldstruct()

name = (*field_arr_[idx]).name
print, "************* make raw maps, " + name
plot_raw_maps, idx, scale=0.17
make_summary_pdf, [idx]

end

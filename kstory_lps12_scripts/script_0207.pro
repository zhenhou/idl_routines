;;;
; NAME: script_0207
; PURPOSE:
;   Make final runlists for all fields
;
; NOTES:
;   calls make_runlist_for_lps12_single_field from make_runlist_for_lps12.pro
;
; MODIFICATION HISTORY:
;  09/21/2011: (KTS) Created
;  02/07/2012: (KTS) Modified to include all runlists
;;;

;...................................................................
; Make runlists
;
pro script_0207_0to4
compile_opt IDL2, HIDDEN
field_idx_list = [0,1,2,3,4]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end

pro script_0207_5to9
compile_opt IDL2, HIDDEN
field_idx_list = [5,6,7,8,9]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end

pro script_0207_10to14
compile_opt IDL2, HIDDEN
field_idx_list = [10,11,12,13,14]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end

pro script_0207_15to19
compile_opt IDL2, HIDDEN
field_idx_list = [15,16,17,18,19]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end




;...................................................................
; Plot raw maps
;
pro plot_0207_11to14
compile_opt IDL2, HIDDEN

field_list = [11,12,13,14]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name
    print, "************* make raw maps, " + name
    plot_raw_maps, idx, scale=0.17
    make_summary_pdf, [idx]
endfor
end

pro plot_0207_15to19
compile_opt IDL2, HIDDEN

;field_list = [15,16,17,18,19]
field_list = [16,17,18,19]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name
    print, "************* make raw maps, " + name
    plot_raw_maps, idx, scale=0.17
    make_summary_pdf, [idx]
endfor
end


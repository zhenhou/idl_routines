;;;
; NAME: script_0921
; PURPOSE:
;   re-run runlists for dirs that have 2010 obs in 2011 dirs
;
; NOTES:
;   calls make_runlist_for_lps12_single_field from make_runlist_for_lps12.pro
;
; MODIFICATION HISTORY:
;  09/20/2011: (KTS) Created
;;;

;...................................................................
; Make runlists
;
pro script_0921_0to4
compile_opt IDL2, HIDDEN
field_idx_list = [0,1,2,3,4]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end

pro script_0921_5to9
compile_opt IDL2, HIDDEN
field_idx_list = [5,6,7,8,9]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end

pro script_0921_10to14
compile_opt IDL2, HIDDEN
field_idx_list = [10,11,12,13,14]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end

pro script_0921_15to16
compile_opt IDL2, HIDDEN
field_idx_list = [15,16]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end




;...................................................................
; Plot raw maps
;
pro plot_0921_0to4
compile_opt IDL2, HIDDEN

field_list = [0,1,2,3,4]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name
    print, "************* make raw maps, " + name
    plot_raw_maps, idx, scale=0.17
    make_summary_pdf, [idx]
endfor
end

pro plot_0921_5to9
compile_opt IDL2, HIDDEN


; already did 7
field_list = [5,6,8,9]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name
    print, "************* make raw maps, " + name
    plot_raw_maps, idx, scale=0.17
    make_summary_pdf, [idx]
endfor
end

pro plot_0921_10to14
compile_opt IDL2, HIDDEN

field_list = [10,11,12,13,14]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name
    print, "************* make raw maps, " + name
    plot_raw_maps, idx, scale=0.17
    make_summary_pdf, [idx]
endfor
end

pro plot_0921_15to16
compile_opt IDL2, HIDDEN

field_list = [15,16]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name
    print, "************* make raw maps, " + name
    plot_raw_maps, idx, scale=0.17
    make_summary_pdf, [idx]
endfor
end


;;;
; NAME: script_0920
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
; Main function
;
pro script_0920_12to16
compile_opt IDL2, HIDDEN
field_idx_list = [12,13,14,15,16]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end


pro script_0920_6to11
compile_opt IDL2, HIDDEN
field_idx_list = [6,7,8,9,10,11]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end


pro script_0920_11
compile_opt IDL2, HIDDEN
field_idx_list = [11]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** make_runlist_for_lps12_single_field, ',field_arr_[field_idx_list[ii]].name
    make_runlist_for_lps12_single_field, field_idx_list[ii]
endfor
end


;...................................................................
; Make maps
;
pro plot_0920_0to4
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

pro plot_0920_5to9
compile_opt IDL2, HIDDEN

field_list = [5,6,7,8,9]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name
    print, "************* make raw maps, " + name
    plot_raw_maps, idx, scale=0.17
    make_summary_pdf, [idx]
endfor
end

pro plot_0920_10to14
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

pro plot_0920_15to16
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


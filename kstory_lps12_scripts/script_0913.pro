;;;
; NAME: script_0913
; PURPOSE:
;   make runlists with new program
;
; NOTES:
;
; MODIFICATION HISTORY:
;  09/12/2011: (KTS) Created
;;;

;field_idx_list = [1, 5, 6, 12] - > already done

;...................................................................
; Main function
;
pro script_0913_0to4
compile_opt IDL2, HIDDEN
field_idx_list = [0,2,3,4]
for ii=0, n_elements(field_idx_list) -1 do begin
    ;print, 'make_runlist_single_field_0913, field_idx_list[',ii
    make_runlist_single_field_0913, field_idx_list[ii]
endfor
end

pro script_0913_5to9
compile_opt IDL2, HIDDEN
field_idx_list = [5,6,7,8,9]
for ii=0, n_elements(field_idx_list) -1 do begin
    make_runlist_single_field_0913, field_idx_list[ii]
endfor
end

pro script_0913_10to14
compile_opt IDL2, HIDDEN
field_idx_list = [10,11,12,13,14]
for ii=0, n_elements(field_idx_list) -1 do begin
    make_runlist_single_field_0913, field_idx_list[ii]
endfor
end

pro script_0913_15to19
compile_opt IDL2, HIDDEN
field_idx_list = [15,16,17,18,19]
for ii=0, n_elements(field_idx_list) -1 do begin
    make_runlist_single_field_0913, field_idx_list[ii]
endfor
end



;...................................................................
; Run on a single field
;
pro make_runlist_single_field_0913, idx
field_arr_ = lps12_fieldstruct()

; The final cutsetup
cutSet = lps12_runlist_cutset()

print, '*** Make runlist for field ' + field_arr_[idx].name
print, '    map dirs: ', field_arr_[idx].map_dirs

make_runlist_for_lps12, 150, dates150, field_arr_[idx].map_dirs, islt=field_arr_[idx].lead_trail,$
  fac_high_medwt=cutSet.high_medwt, fac_low_medwt=cutSet.low_medwt, $
  fac_high_rms=cutSet.high_rms, fac_low_rms=cutSet.low_rms, $
  fac_high_wtrms=cutSet.high_wtrms, fac_low_wtrms=cutSet.low_wtrms, $
  fac_high_tweight=cutSet.high_tweight, fac_low_tweight=cutSet.low_tweight, $
  outfile='/data/kstory/projects/lps12/runlists/runlist_lps12_'+field_arr_[idx].name+'.txt',$
  savfile='/data/kstory/projects/lps12/runlists/sav_files/runlist_lps12_'+field_arr_[idx].name+'.sav';, /stopit

end

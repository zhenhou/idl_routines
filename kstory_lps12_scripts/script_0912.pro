;;;
; NAME: script_0912
; PURPOSE:
;   make runlists with new program
;
; NOTES:
;
; MODIFICATION HISTORY:
;  09/12/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro script_0912
compile_opt IDL2, HIDDEN

field_arr_ = lps12_fieldstruct()

; The final cutsetup
cutSet = {high_medwt:1.25,  $
          low_medwt:2.0,    $
          high_rms:2.0,     $
          low_rms:3.0,      $
          high_wtrms:2.0,   $
          low_wtrms:2.0,    $
          high_tweight:1.2, $
          low_tweight:1.5}


exclude_list = [-1]
;field_idx_list = setdifference(lindgen(n_elements(all_fields)), exclude_list)
field_idx_list = [1, 5, 6, 12]

for ii=0, n_elements(field_idx_list)-1 do begin
    idx = field_idx_list[ii]

    print, '*** Make runlist for field ' + field_arr_[idx].name
    print, '    map dirs: ', field_arr_[idx].map_dirs

    make_runlist_for_lps12, 150, dates150, field_arr_[idx].map_dirs, islt=field_arr_[idx].lead_trail,$
      fac_high_medwt=cutSet.high_medwt, fac_low_medwt=cutSet.low_medwt, $
      fac_high_rms=cutSet.high_rms, fac_low_rms=cutSet.low_rms, $
      fac_high_wtrms=cutSet.high_wtrms, fac_low_wtrms=cutSet.low_wtrms, $
      fac_high_tweight=cutSet.high_tweight, fac_low_tweight=cutSet.low_tweight, $
      outfile='/data/kstory/projects/lps12/runlists/runlist_lps12_'+field_arr_[idx].name+'.txt',$
      savfile='/data/kstory/projects/lps12/runlists/sav_files/runlist_lps12_'+field_arr_[idx].name+'.sav', /stopit
endfor

print, '*** Done with runlist-making'
    
end



;...................................................................
; Run on a single field
;
pro make_runlist_single_field_0912, idx
field_arr_ = lps12_fieldstruct()

; The final cutsetup
cutSet = {high_medwt:1.25,  $
          low_medwt:2.0,    $
          high_rms:2.0,     $
          low_rms:3.0,      $
          high_wtrms:2.0,   $
          low_wtrms:2.0,    $
          high_tweight:1.2, $
          low_tweight:1.5}

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

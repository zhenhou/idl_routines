;;;
; NAME: make_runlist_xspec_lps12
; PURPOSE:
;   Make runlists of maps that we plan to cross in the power spectrum analysis
;
; INPUTS:
;   field_idx,       field index from lps12_fieldstruct()
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  03/01/2012: (KTS) Created
;;;

pro run_all
for ii=0, 19 do begin
    make_runlist_dateOnly_lps12, ii
endfor
end

;...................................................................
;
pro make_runlist_xspec_lps12, field_idx

field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
is_lt = fst.lead_trail
runlist = fst.obs_runlist

runlist_dir = '/data/kstory/projects/lps12/runlists/'
;runlist_dir = '/data/kstory/projects/lps12/runlists/obsolete/'
new_runlist = runlist_dir + 'runlist_xspec_lps12_'+field_name+'.txt'


; If this is a normal field, runlist_dateOnly will be the same as runlist_xspec
if ~is_lt then begin

    print, "Make xspec runlist for field: ", new_runlist
    spawn, 'cp ' + runlist + ' ' + new_runlist

; if this is lead-trail, then use the maps in the lt directory
endif else begin
    print, "Make xspec runlist for LT field: ", new_runlist

    xspec_dir = fst.xspec_map_dir
    spawn, 'ls ' + xspec_dir + '*.fits', maps
    nmaps = n_elements(maps)
    date_list = extract_date_from_filename(maps)

    ; write out new runlist
    get_lun, lun1
    openw, lun1, new_runlist
    for ii=0, n_elements(date_list)-1 do begin
        printf, lun1, date_list[ii]
    endfor
    close, lun1
    free_lun,lun1

endelse

end

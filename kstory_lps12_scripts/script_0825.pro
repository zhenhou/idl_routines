;;;
; NAME: script_0825
; PURPOSE:
;   Try getting total(data.map.map)
;
; MODIFICATION HISTORY:
;  08/25/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro script_0825, field_name
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'

field_idx = where(all_fields eq field_name, nwh)
if(nwh le 0) then begin
    print, "field_name ["+field_name+"] does not exist!  exit..."
    return
endif

files = file_search(all_map_dirs[field_idx] + '/map_*150_*.fits',count=nfiles)
if nfiles eq 0 then begin
    print,"No map files for 150 GHz in directory "+all_map_dirs[field_idx]+".  Quitting."
endif
print,"SCRIPT_0825.PRO: "+string(nfiles)+" map files found for 150 GHz in directory "+all_map_dirs[field_idx]+"."

total_map = fltarr(nfiles)
dates_new = fltarr(nfiles)

for i=0,nfiles-1 do begin
    print, 'Restore map ',strtrim(string(i),1),' out of ', strtrim(string(nfiles-1),1), ': '+files[i]
    mstemp = expand_fits_struct(read_spt_fits(files[i]))
    dates_new[i] = mstemp.mapinfo.start_time
    total_map[i] = total(mstemp.weight.map)

    ;tv_spt_map, mstemp.weight.map
endfor

restore, '/data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_'+field_name+'.sav'

save, total_map, dates_new, $
  DATES, FIELD_NAME, MAXWTS_STRIPE, MEDWTS, RMS_1AM_MIDMAP, $
  file = '/data/kstory/projects/lps12/runlists/cut_map_total/runlist_cut_map_total_'+field_name+'.sav'

end

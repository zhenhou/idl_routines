;;;
; NAME: make_final_runlists
; PURPOSE:
;   Make runlist for lps12
;
; OUTPUT:
;   Output file: /data/kstory/projects/lps12/runlists/final/runlist_final_[field_name].txt
;
; NOTES:
;   - This code expects sav files to be in
;        '/data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_'+field_name+'.sav'
;   - If we need to re-make the sav files, set remake_sav
;
; MODIFICATION HISTORY:
;  08/31/2011: (KTS) Created
;;;

;...................................................................
; Run the main program
;
pro run_make_final_runlists
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'
;exclude_list = [1,2,5,11]
;exclude_list = [0,1,2,5,11]
exclude_list = [-1]
field_idx_list = setdifference(lindgen(n_elements(all_fields)), exclude_list)
;field_idx_list = [15]

for ii=0, n_elements(field_idx_list)-1 do begin
    idx = field_idx_list[ii]
    field_name = all_fields[idx]
    print, "************* make_final_runlists, "+field_name
    make_final_runlists, field_name
endfor
stop
end

;...................................................................
; Main function
;
pro make_final_runlists, field_name, remake_sav=remake_sav, remake_tweight_sav=remake_tweight_sav
compile_opt IDL2, HIDDEN

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Step 0: Setup
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; The final cutsetup
cutSet = {high_medwt:1.25,  $;high_medwt:2.0,   $
          low_medwt:2.0,    $
          high_rms:2.0,     $
          low_rms:3.0,      $
          high_wtrms:2.0,   $
          low_wtrms:2.0,    $
          high_tweight:1.2, $
          low_tweight:1.5}

restore, 'lps12_fieldnames.sav'
;cutSetNames = cut_set_names()
idx = where(all_fields eq field_name, nwh)
if(nwh le 0) then begin
    print, "field_name ["+field_name+"] does not exist!  exit..."
    return
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Step 1: Make default sav files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; If we need to re-make the .sav file, do so
if keyword_set(remake_sav) then begin

    make_default_sav_files, field_name

    ; We will also need to remake the tweight sav
    remake_tweight_sav = 1 

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Step 2: Make tweight sav files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; If we need to re-make tweight_sav
if keyword_set(remake_tweight_sav) then begin
    add_tweight_to_sav, field_name
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Step 3: Make tweight sav files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

make_runlist_from_tweight_sav, field_name, cutSet;, outdir='/data/kstory/projects/lps12/runlists/medwts_1.25/'

end

;...................................................................
; Get the actual runlists
;
pro make_runlist_from_tweight_sav, field_name, cuts, outdir=outdir
compile_opt IDL2, HIDDEN

print, "--- Make final runlist for field: ",field_name
restore, '/data/kstory/projects/lps12/runlists/tweight_sav/runlist_tweight_'+field_name+'.sav'
rms = rms_1am_midmap
wtrms = medwts*rms^2.

;help, medwts, rms, wtrms, tweight ; DEBUGGING
; Apply cuts
wh = where( medwts lt median(medwts)*cuts.high_medwt and $
            medwts gt median(medwts)/cuts.low_medwt and $
            rms   lt median(rms)*cuts.high_rms and $
            rms   gt median(rms)/cuts.low_rms and $
            wtrms lt median(wtrms)*cuts.high_wtrms and $
            wtrms gt median(wtrms)/cuts.low_wtrms and $
            tweight lt median(tweight)*cuts.high_tweight and $
            tweight gt median(tweight)/cuts.low_tweight, nwh)

; Get dates
dates150 = dates[wh]

; Make runfile
if ~keyword_set(outdir) then outdir = '/data/kstory/projects/lps12/runlists/final'

outname = 'runlist_final_'+field_name+'.txt'
outname_sav = 'runlist_final_'+field_name+'.sav'

outfile = outdir + '/'+outname
get_lun, lun1
openw, lun1, outfile
for ii=0, n_elements(dates150)-1 do begin
    printf, lun1, dates150[ii]
endfor
close, lun1
free_lun, lun1

; Make .sav file

print, "   Total number of maps in this directory: ", n_elements(dates)
print, "   Number of maps that pass cuts:          ", n_elements(dates150)
print, "   --- Save data to " + outdir
savfile = outdir + '/' + outname_sav
save, dates150, filename=savfile

end

;...................................................................
; Make default_sav files
; OUTPUT: saved in /data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_[field_name].sav
;
pro make_default_sav_files, field_name
compile_opt IDL2, HIDDEN
restore, 'lps12_fieldnames.sav'
idx = where(all_fields eq field_name, nwh)

if (lead_trail[idx]) then begin
    print, 'Make default runlist for field ', field_name, ' with Lead-Trail'
    make_runlist_with_cuts, 150, dates150, mapdir=all_map_dirs[idx],$
      fac_high_medwt=2., fac_high_wtrms=2., /islt,$
      outfile='/data/kstory/projects/lps12/runlists/default_sav/runlist_wtrms_'+field_name+'.txt',$
      savfile='/data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_'+field_name+'.sav'
endif else begin
    print, 'Make default runlist for field ', field_name, ', Lead-Trail OFF'
    make_runlist_with_cuts, 150, dates150, mapdir=all_map_dirs[idx],$
      fac_high_medwt=2., fac_high_wtrms=2., $
      outfile='/data/kstory/projects/lps12/runlists/default_sav/runlist_wtrms_'+field_name+'.txt',$
      savfile='/data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_'+field_name+'.sav';,/stopit
endelse
end

;...................................................................
; Add tweight info to runlist.sav files
; OUTPUT: saved in /data/kstory/projects/lps12/runlists/tweight_sav/runlist_tweight_[field_name].sav
;
pro add_tweight_to_sav, field_name
compile_opt IDL2, HIDDEN
restore, 'lps12_fieldnames.sav'
field_idx = where(all_fields eq field_name, nwh)

files = file_search(all_map_dirs[field_idx] + '/map_*150_*.fits',count=nfiles)
if nfiles eq 0 then begin
    print,"No map files for 150 GHz in directory "+all_map_dirs[field_idx]+".  Quitting."
endif
print,"SCRIPT_0825.PRO: "+string(nfiles)+" map files found for 150 GHz in directory "+all_map_dirs[field_idx]+"."

dates_new = fltarr(nfiles)
tweight   = fltarr(nfiles)

for i=0,nfiles-1 do begin
    print, 'Restore map ',strtrim(string(i),1),' out of ', strtrim(string(nfiles-1),1), ': '+files[i]
    mstemp = expand_fits_struct(read_spt_fits(files[i]))
    dates_new[i] = mstemp.mapinfo.start_time
    tweight[i] = total(mstemp.weight.map)

    ;tv_spt_map, mstemp.weight.map
endfor

restore, '/data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_'+field_name+'.sav'

save, tweight, dates_new, $
  DATES, FIELD_NAME, MAXWTS_STRIPE, MEDWTS, RMS_1AM_MIDMAP, $
  file = '/data/kstory/projects/lps12/runlists/tweight_sav/runlist_tweight_'+field_name+'.sav'

end


;...................................................................
; Make list of all observations of this field
;
pro make_obs_list, field_name, remake_sav=remake_sav, remake_tweight_sav=remake_tweight_sav
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'

field_idx = where(all_fields eq field_name, nwh)

files = file_search(all_map_dirs[field_idx] + '/map_*150_*.fits',count=nfiles)

save, files, $
  file = '/data/kstory/projects/lps12/runlists/obs_list_'+field_name+'.sav'


end

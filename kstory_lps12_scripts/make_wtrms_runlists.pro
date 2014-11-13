;;;
; NAME: make_wtrms_runlists
; PURPOSE:
;   Make runlist with specified cuts
;
; NOTES:
;   - This code expects sav files to be in
;        '/data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_'+field_name+'.sav'
;   - If we need to re-make the sav files, set remake_sav
;
;
; MODIFICATION HISTORY:
;  08/21/2011: (KTS) Created
;;;

;...................................................................
; Main function to make runlists
;
pro make_wtrms_runlists, field_name, remake_sav=remake_sav, remake_tweight_sav=remake_tweight_sav
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'
cutSetNames = cut_set_names()
idx = where(all_fields eq field_name)

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

for cc=0, n_elements(cutSetNames)-1 do begin
    make_runlist_with_specific_cut, field_name, cutSetNames[cc]
endfor
end

;...................................................................
; Get the actual runlists
;
pro make_runlist_with_specific_cut, field_name, cut_name
compile_opt IDL2, HIDDEN

print, "Make runlist for field: ",field_name," --- cut: ",cut_name

;restore, '/data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_'+field_name+'.sav'
restore, '/data/kstory/projects/lps12/runlists/tweight_sav/runlist_tweight_'+field_name+'.sav'
rms = rms_1am_midmap
wtrms = medwts*rms^2.

; Get full cut set
cutSet = cut_set()
if (cut_name eq 'default2') then cuts = cutSet.default2
if (cut_name eq 'high_medwt_test1p5') then cuts = cutSet.high_medwt_test1p5
if (cut_name eq 'high_medwt_test1p5') then cuts = cutSet.high_medwt_test1p5
if (cut_name eq 'low_medwt_test3p0') then cuts = cutSet.low_medwt_test3p0
if (cut_name eq 'low_medwt_test1p5') then cuts = cutSet.low_medwt_test1p5
if (cut_name eq 'low_medwt_test1p3') then cuts = cutSet.low_medwt_test1p3

if (cut_name eq 'high_rms_test1p5') then cuts = cutSet.high_rms_test1p5
if (cut_name eq 'low_rms_test3p0') then cuts = cutSet.low_rms_test3p0
if (cut_name eq 'low_rms_test1p5') then cuts = cutSet.low_rms_test1p5

if (cut_name eq 'high_wtrms_test1p5') then cuts = cutSet.high_wtrms_test1p5
if (cut_name eq 'high_wtrms_test1p3') then cuts = cutSet.high_wtrms_test1p3
if (cut_name eq 'low_wtrms_test1p5') then cuts = cutSet.low_wtrms_test1p5
if (cut_name eq 'low_wtrms_test1p3') then cuts = cutSet.low_wtrms_test1p3

if (cut_name eq 'high_tweight_test1p5') then cuts = cutSet.high_tweight_test1p5
if (cut_name eq 'high_tweight_test1p3') then cuts = cutSet.high_tweight_test1p3
if (cut_name eq 'high_tweight_test1p2') then cuts = cutSet.high_tweight_test1p2
if (cut_name eq 'low_tweight_test1p5') then cuts = cutSet.low_tweight_test1p5
if (cut_name eq 'low_tweight_test1p3') then cuts = cutSet.low_tweight_test1p3

; Apply cuts
wh = where( medwts lt median(medwts)*cuts.high_medwt and $
            medwts gt median(medwts)/cuts.low_medwt or $
            rms   lt median(rms)*cuts.high_rms and $
            rms   gt median(rms)/cuts.low_rms and $
            wtrms lt median(wtrms)*cuts.high_wtrms and $
            wtrms gt median(wtrms)/cuts.low_wtrms and $
            tweight lt median(tweight)*cuts.high_tweight and $
            tweight gt median(tweight)/cuts.low_tweight, nwh)

; Get dates
dates150 = dates[wh]

; Make runfile
outname = 'runlist_wtrms_'+cut_name+'_'+field_name+'.txt'
outname_sav = 'runlist_wtrms_'+cut_name+'_'+field_name+'.sav'

outfile = '/data/kstory/projects/lps12/runlists/'+field_name+'/'+outname
get_lun, lun1
openw, lun1, outfile
for ii=0, n_elements(dates150)-1 do begin
    printf, lun1, dates150[ii]
endfor
close, lun1
free_lun, lun1

; Make .sav file
savfile = '/data/kstory/projects/lps12/runlists/'+field_name+'/'+outname_sav
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
      savfile='/data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_'+field_name+'.sav'
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

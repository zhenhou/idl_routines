;;;
; NAME: plots_0809
; PURPOSE:
;   Test medwt cut on 2011 fields
;
; MODIFICATION HISTORY:
;  08/07/2011: (KTS) Created
;;;

;...................................................................
; Call runlist and coadd functions on all maps of interest
;
pro runlist_wtrms_cut_test, save_plots=save_plots, stopit=stopit
compile_opt IDL2, HIDDEN

; Field names
restore, 'lps12_fieldnames.sav'

; One field from each year
field_list = [ 1, 2, 5, 11]

;-------------------
; high medwt test
;-------------------
high_medwt = [1.5, 2., 3.]

for ii=0, n_elements(field_list) - 1 do begin 

    field_idx = field_list[ii]
    field_name = all_fields[idx]
    field_dir = all_map_dirs[idx]
    

    for jj=0, 2 do begin
        cuts = {high_medwt:high_medwt,$
                low_medwt:2.,$
                high_rms:2.,$
                low_rms:2.,$
                high_wtrms:2.,$
                low_wtrms:2.}
        output_fname = 'runlist_cuts_medwt_high'+string(cuts.high_medwt, format='(F3.1)')+'_low'+string(cuts.low_medwt, format='(F3.1)')+'_' +field_name
        output_coadd_fname = 'runlist_cuts_coadd_medwt_high'+string(cuts.high_medwt, format='(F3.1)')+'_low'+string(cuts.low_medwt, format='(F3.1)')+'_' +field_name
        
        ; Run
        print_test, field_idx, cuts, output_fname, output_coadd_fname
        ;make_runlist_and_coadd, field_idx, cuts, output_fname, output_coadd_fname
    endfor
endfor

;-------------------
; low medwt test
;-------------------
low_medwt = [1.5, 2., 3.]

for ii=0, n_elements(field_list) - 1 do begin 

    field_idx = field_list[ii]

    for jj=0, 2 do begin
        cuts = {high_medwt:2.,$
                low_medwt:low_medwt,$
                high_rms:2.,$
                low_rms:2.,$
                high_wtrms:2.,$
                low_wtrms:2.}
        output_fname = 'runlist_cuts_medwt_high'+string(cuts.high_medwt, format='(F3.1)')+'_low'+string(cuts.low_medwt, format='(F3.1)')+'_' +field_name
        output_coadd_fname = 'runlist_cuts_coadd_medwt_high'+string(cuts.high_medwt, format='(F3.1)')+'_low'+string(cuts.low_medwt, format='(F3.1)')+'_' +field_name
        
    ; Run
        print_test, field_idx, cuts, output_fname, output_coadd_fname
    ;make_runlist_and_coadd_field_idx, cuts, output_fname, output_coadd_fname
    endfor
endfor

;--------------------------------------------------------------------------------

;-------------------
; high rms test
;-------------------
high_rms = [1.5, 2., 3.]

for ii=0, n_elements(field_list) - 1 do begin 

    field_idx = field_list[ii]

    for jj=0, 2 do begin
        cuts = {high_medwt:2.,$
                low_medwt:2.,$
                high_rms:high_rms,$
                low_rms:2.,$
                high_wtrms:2.,$
                low_wtrms:2.}
        output_fname = 'runlist_cuts_rms_high'+string(cuts.high_rms, format='(F3.1)')+'_low'+string(cuts.low_rms, format='(F3.1)')+'_' +field_name
        output_coadd_fname = 'runlist_cuts_coadd_rms_high'+string(cuts.high_rms, format='(F3.1)')+'_low'+string(cuts.low_rms, format='(F3.1)')+'_' +field_name
        
    ; Run
        print_test, field_idx, cuts, output_fname, output_coadd_fname
    ;make_runlist_and_coadd_field_idx, cuts, output_fname, output_coadd_fname
    endfor
endfor
;-------------------
; low rms test
;-------------------
low_rms = [1.5, 2., 3.]

for ii=0, n_elements(field_list) - 1 do begin 

    field_idx = field_list[ii]

    for jj=0, 2 do begin
        cuts = {high_medwt:2.,$
                low_medwt:2.,$
                high_rms:2.,$
                low_rms:low_rms,$
                high_wtrms:2.,$
                low_wtrms:2.}
        output_fname = 'runlist_cuts_rms_high'+string(cuts.high_rms, format='(F3.1)')+'_low'+string(cuts.low_rms, format='(F3.1)')+'_' +field_name
        output_coadd_fname = 'runlist_cuts_coadd_rms_high'+string(cuts.high_rms, format='(F3.1)')+'_low'+string(cuts.low_rms, format='(F3.1)')+'_' +field_name
        
    ; Run
        print_test, field_idx, cuts, output_fname, output_coadd_fname
    ;make_runlist_and_coadd_field_idx, cuts, output_fname, output_coadd_fname
    endfor
endfor

;--------------------------------------------------------------------------------

;-------------------
; high wtrms test
;-------------------
high_wtrms = [1.5, 2., 3.]

for ii=0, n_elements(field_list) - 1 do begin 

    field_idx = field_list[ii]

    for jj=0, 2 do begin
        cuts = {high_medwt:2.,$
                low_medwt:2.,$
                high_rms:2.,$
                low_rms:2.,$
                high_wtrms:high_wtrms,$
                low_wtrms:2.}
        output_fname = 'runlist_cuts_rms_high'+string(cuts.high_rms, format='(F3.1)')+'_low'+string(cuts.low_rms, format='(F3.1)')+'_' +field_name
        output_coadd_fname = 'runlist_cuts_coadd_rms_high'+string(cuts.high_rms, format='(F3.1)')+'_low'+string(cuts.low_rms, format='(F3.1)')+'_' +field_name
        
    ; Run
        print_test, field_idx, cuts, output_fname, output_coadd_fname
    ;make_runlist_and_coadd_field_idx, cuts, output_fname, output_coadd_fname
    endfor
endfor

;-------------------
; low rms test
;-------------------
low_wtrms = [1.5, 2., 3.]

for ii=0, n_elements(field_list) - 1 do begin 

    field_idx = field_list[ii]

    for jj=0, 2 do begin
        cuts = {high_medwt:2.,$
                low_medwt:2.,$
                high_rms:2.,$
                low_rms:2.,$
                high_wtrms:2.,$
                low_wtrms:low_wtrms}
        output_fname = 'runlist_cuts_rms_high'+string(cuts.high_rms, format='(F3.1)')+'_low'+string(cuts.low_rms, format='(F3.1)')+'_' +field_name
        output_coadd_fname = 'runlist_cuts_coadd_rms_high'+string(cuts.high_rms, format='(F3.1)')+'_low'+string(cuts.low_rms, format='(F3.1)')+'_' +field_name
        
    ; Run
        print_test, field_idx, cuts, output_fname, output_coadd_fname
    ;make_runlist_and_coadd_field_idx, cuts, output_fname, output_coadd_fname
    endfor
endfor

end


;...................................................................
; Make the runlist and coadd
;
pro make_runlist_and_coadd, field_idx, cuts, output_fname, output_coadd_fname
compile_opt IDL2, HIDDEN

idx = field_idx                 ; get index in all_fields
field_name = all_fields[idx]
field_dir = all_map_dirs[idx]

print, field_name
restore, '/data/kstory/projects/lps12/runlists/runlist_cuts_wtrms_'+field_name+'.sav'

                                ; Make runlist
if (lead_trail[idx]) then begin
    print, 'Make runlist for field ', field_name, ' with Lead-Trail'
    make_runlist_with_cuts, 150, dates150, mapdir=all_map_dirs[idx],$
      fac_high_medwt=cuts.high_medwt, fac_low_medwt=cuts.low_medwt, $
      fac_high_rms=cuts.high_rms, fac_low_rms=cuts.low_rms, $
      fac_high_wtrms=cuts.high_wtrms, fac_low_wtrms=cuts.low_wtrms, $
      /islt,$
      outfile='/data/kstory/projects/lps12/runlists/'+output_fname+'.txt',$
      savfile='/data/kstory/projects/lps12/runlists/'+output_fname+'.sav'
    
endif else begin
    print, 'Make runlist for field ', field_name, ', Lead-Trail OFF'
    make_runlist_with_cuts, 150, dates150, mapdir=all_map_dirs[idx],$
      fac_high_medwt=cuts.high_medwt, fac_low_medwt=cuts.low_medwt, $
      fac_high_rms=cuts.high_rms, fac_low_rms=cuts.low_rms, $
      fac_high_wtrms=cuts.high_wtrms, fac_low_wtrms=cuts.low_wtrms, $
      outfile='/data/kstory/projects/lps12/runlists/'+output_fname+'.txt',$
      savfile='/data/kstory/projects/lps12/runlists/'+output_fname+'.sav'
endelse

    ; Make coadd
coadd_fits_maps,dates=dates, dir=field_dir, $
  band=150,fileout='/data/kstory/projects/lps12/runlists/coadds/'+output_coadd_fname+'.fits'

end


;...................................................................
; Print out what is going to be called
;
pro print_test
    print_test, field_idx, cuts, output_fname, output_coadd_fname

    print, "Field"
    print, "   ", all_fields[field_idx]
    print, "cuts:"
    print, "   cuts.high_medwt = ",cuts.high_medwt
    print, "   cuts.low_medwt = ",cuts.low_medwt
    print, "   cuts.high_rms = ",cuts.high_rms
    print, "   cuts.low_rms = ",cuts.low_rms
    print, "   cuts.high_wtrms = ",cuts.high_wtrms
    print, "   cuts.low_wtrms = ",cuts.low_wtrms
    print, "output_fname:"
    print, "   ", output_fname
    print, "output_coadd_fname:"
    print, "   ", output_coadd_fname
end

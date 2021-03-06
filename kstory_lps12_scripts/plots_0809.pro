;;;
; NAME: plots_0809
; PURPOSE:
;   Test medwt cut on 2011 fields
;
; MODIFICATION HISTORY:
;  08/07/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro plots_0809, save_plots=save_plots, stopit=stopit
compile_opt IDL2, HIDDEN

; Field names
restore, 'lps12_fieldnames.sav'

; 2009 fields

for idx=0, 2 do begin 
    field_name = fields_2011[idx]
    field_dir = map_dirs_2011[idx]

    print, field_name
    restore, '/data/kstory/projects/lps12/runlists/runlist_cuts_wtrms_'+field_name+'.sav'


    ;;; medwt test
    high_medwt = [1.5, 2., 3.]
    
    for jj = 0, 2 do begin
        fac_high_medwt = high_medwt[jj]
        fac_low_medwt = 2.

        output_fname = 'runlist_cuts_medwt_high'+fac_high_medwt+'_low'+fac_low_medwt+'_' +field_name
        output_coadd_fname = 'runlist_cuts_coadd_medwt_high'+fac_high_medwt+'_low'+fac_low_medwt+'_' +field_name
        fac_high_rms=2. & fac_low_rms=2. & fac_high_wtrms=2. & fac_low_wtrms=2.

        ; Make runlist
        if (lead_trail[idx]) then begin
            print, 'Make runlist for field ', field_name, ' with Lead-Trail'
            make_runlist_with_cuts, 150, dates150, mapdir=all_map_dirs[idx],$
              fac_high_medwt=fac_high_medwt, fac_low_medwt=fac_low_medwt, $
              fac_high_rms=fac_high_rms, fac_low_rms=fac_low_rms, $
              fac_high_wtrms=fac_high_wtrms, fac_low_wtrms=fac_low_wtrms, $
              /islt,$
              outfile='/data/kstory/projects/lps12/runlists/'+output_fname+'.txt',$
              savfile='/data/kstory/projects/lps12/runlists/'+output_fname+'.sav'
            
        endif else begin
            print, 'Make runlist for field ', field_name, ', Lead-Trail OFF'
            make_runlist_with_cuts, 150, dates150, mapdir=all_map_dirs[idx],$
              fac_high_medwt=fac_high_medwt, fac_low_medwt=fac_low_medwt, $
              fac_high_rms=fac_high_rms, fac_low_rms=fac_low_rms, $
              fac_high_wtrms=fac_high_wtrms, fac_low_wtrms=fac_low_wtrms, $
              outfile='/data/kstory/projects/lps12/runlists/'+output_fname+'.txt',$
              savfile='/data/kstory/projects/lps12/runlists/'+output_fname+'.sav'
        endelse

        ; Make coadd
        coadd_fits_maps,dates=dates, dir=field_dir, $
          band=150,fileout='/data/kstory/projects/lps12/runlists/coadds/'+output_coadd_fname+'.fits'

        if keyword_set(stopit) then stop
    endfor
endfor

end

;;;
; NAME: coadd_maps_0305
; PURPOSE:
;   Coadd maps
;
; INPUTS
;   coadd_dir,             ouput directory for where to write the coadds
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/27/2012: (KTS) Created
;  03/12/2012: (KTS) Use xspec_maps runlist
;;;


;...................................................................
; coadd all maps from one field from the maps/20120305/ directories
;
PRO coadd_maps_0305, field_idx, $
                     coadd_dir = coadd_dir

compile_opt IDL2, HIDDEN

;------------------------
; Get the field information
;------------------------
field_arr = lps12_fieldstruct()
field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name

;------------------------
; 50 mJy ptsrc cut
;------------------------
coadd_dir = '/data/kstory/projects/lps12/maps/20120305/coadds/'
coadd_name = coadd_dir+'coadd_'+field_name+'_50mJy.fits'

; If coadd already exists, return
if file_test(coadd_name) then begin
    print, "Coadd already exists: ", coadd_name
    print, "Returning."
    return
endif


;------------------------
; Make a list of the fits files to coadd
;------------------------

maps = get_lps12_runlist(field_idx, /xspec_maps)
ndates = n_elements(maps)

; write the file list out to a temporary file
get_lun, lun1
fname ='tmp_runlist.txt' 
openw, lun1, fname
for ii=0, ndates-1 do begin
    printf, lun1, maps[ii]
endfor
close, lun1
free_lun, lun1


;------------------------
; Coadd the maps
;------------------------
print, 'COADD_MAPS_0305: output = ', coadd_name
coadd_fits_maps, fits_list=fname, fileout=coadd_name


; clean up
spawn, 'rm -f '+fname
end




;...................................................................
; make coadds for all fields
;
pro run_0305, field_idx
field_list = indgen(20)

; skip fields in the exclude list
if 0 then begin
exclude_list = [0,1,3]
for iex=0, n_elements(exclude_list) - 1 do begin
    field_list = field_list( where(field_list ne exclude_list[iex]) )
endfor
endif


;for ii=0, n_elements(field_list) -1 do begin
for ii=2, n_elements(field_list) -1 do begin
    if (ii eq 0) then begin
        coadd_maps_pre_29Apr2008_ra5h30dec_55
    endif else coadd_maps_0305, field_list[ii]
endfor
end


;...................................................................
; Make a coadd of ra5h30dec-55_2008
; ra5h30dec-55 had two different scans with different coverage.  use
; the more conservative coverage area.
;
pro coadd_maps_pre_29Apr2008_ra5h30dec_55

;------------------------
; Get the field information
;------------------------
field_idx = 0
field_arr = lps12_fieldstruct()
field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name

;------------------------
; 50 mJy ptsrc cut
;------------------------
coadd_dir = '/data/kstory/projects/lps12/maps/20120305/coadds/'
coadd_name = coadd_dir+'pre_29Apr2008_coadd_'+field_name+'_50mJy.fits'

; If coadd already exists, return
if file_test(coadd_name) then begin
    print, "Coadd already exists: ", coadd_name
    print, "Returning."
    return
endif


;------------------------
; Make a list of the fits files to coadd
;------------------------

dates  = get_lps12_runlist(field_idx, /xspec_dates)
maps   = get_lps12_runlist(field_idx, /xspec_maps)
ndates = n_elements(dates)
wh_good_date = intarr(ndates)

; cut dates down to only before 29Apr2008
mjd_cutoff = date_toolkit('20080429_000000', 'mjd')
for ii=0, ndates -1 do begin
    mjd_tmp = date_toolkit(dates[ii], 'mjd')

    if mjd_tmp le mjd_cutoff then wh_good_date[ii] = 1
endfor

wh_good = where(wh_good_date eq 1)
dates = dates[ wh_good ]
maps  = maps[ wh_good ]
ndates = n_elements(dates)

; write the file list out to a temporary file
get_lun, lun1
fname ='tmp_runlist.txt' 
openw, lun1, fname
for ii=0, ndates-1 do begin
;for ii=0, 9 do begin
    printf, lun1, maps[ii]
endfor
close, lun1
free_lun, lun1


;------------------------
; Coadd the maps
;------------------------
coadd_fits_maps, fits_list=fname, fileout=coadd_name


; clean up
spawn, 'rm -f '+fname
end


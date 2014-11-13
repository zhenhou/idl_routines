;;;
; NAME: coadd_maps_0224
; PURPOSE:
;   Coadd maps
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/27/2012: (KTS) Created
;;;


;...................................................................
; coadd all maps from one field from the maps/20120224/ directories
;
pro coadd_maps_0224, field_idx
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
map_dir = '/data/kstory/projects/lps12/maps/20120224/'+field_name+'/'
coadd_dir = '/data/kstory/projects/lps12/maps/20120224/coadds/'
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

dates  = get_lps12_runlist(field_idx, /obs_dates)
ndates = n_elements(dates)
files  = strarr(ndates)
wh_has_map = intarr(ndates) + 1

for ii=0, ndates-1 do begin
    fits_map = map_dir+'map_'+field_dir_name+'_150_'+dates[ii]+'.fits'

    ; skip maps that don't exist
    if ~file_test(fits_map) then begin
        print, "Failed to find map:", fits_map
        wh_has_map[ii] = 0
    endif

    files[ii] = fits_map
endfor

wh = where(wh_has_map eq 1, nwh)
if (nwh ne 0) then files = files[wh]

; write the file list out to a temporary file
get_lun, lun1
fname ='tmp_runlist.txt' 
openw, lun1, fname
for ii=0, ndates-1 do begin
;for ii=0, 9 do begin
    printf, lun1, files[ii]
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




;...................................................................
; make coadds for all fields
;
pro run_0224, field_idx
field_list = indgen(20)
exclude_list = [0,1,3]

; skip fields in the exclude list
for iex=0, n_elements(exclude_list) - 1 do begin
    field_list = field_list( where(field_list ne exclude_list[iex]) )
endfor

for ii=0, n_elements(field_list) -1 do begin
    coadd_maps_0224, field_list[ii]
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
map_dir = '/data/kstory/projects/lps12/maps/20120224/'+field_name+'/'
coadd_dir = '/data/kstory/projects/lps12/maps/20120224/coadds/'
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

dates  = get_lps12_runlist(field_idx)
ndates = n_elements(dates)
wh_good_date = intarr(ndates)

; cut dates down to only before 29Apr2008
mjd_cutoff = date_toolkit('20080429_000000', 'mjd')
for ii=0, ndates -1 do begin
    mjd_tmp = date_toolkit(dates[ii], 'mjd')

    if mjd_tmp le mjd_cutoff then wh_good_date[ii] = 1
endfor

dates = dates[ where(wh_good_date eq 1) ]
ndates = n_elements(dates)
stop



ndates = n_elements(dates)
files  = strarr(ndates)
wh_has_map = intarr(ndates) + 1

for ii=0, ndates-1 do begin
    fits_map = map_dir+'map_'+field_dir_name+'_150_'+dates[ii]+'.fits'

    ; skip maps that don't exist
    if ~file_test(fits_map) then begin
        print, "Failed to find map:", fits_map
        wh_has_map[ii] = 0
    endif

    files[ii] = fits_map
endfor

wh = where(wh_has_map eq 1, nwh)
if (nwh ne 0) then files = files[wh]

; write the file list out to a temporary file
get_lun, lun1
fname ='tmp_runlist.txt' 
openw, lun1, fname
for ii=0, ndates-1 do begin
;for ii=0, 9 do begin
    printf, lun1, files[ii]
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


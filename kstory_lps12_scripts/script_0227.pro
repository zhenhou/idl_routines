;;;
; NAME: script_0227
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
; main
pro coadd_maps_0224, field_idx

;------------------------
; Get the field information
;------------------------
field_arr = lps12_fieldstruct()
field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name
;idf_dir = field_arr[field_idx].idf_lpfds_dirs
idf_dir = '/home/rkeisler/lowellfits/'+field_name+'/'


;------------------------
; 50 mJy ptsrc cut
;------------------------
map_dir = '/data/kstory/projects/lps12/maps/20120224/'+field_name+'/'
coadd_dir = '/data/kstory/projects/lps12/maps/20120224/coadds/'

;;;
dates  = get_lps12_runlist(field_idx)
ndates = n_elements(dates)
files  = strarr(ndates)
wh_has_map = intarr(ndates) + 1

for ii=0, ndates-1 do begin
    fits_map = map_dir+'map_'+field_name+'_150_'+dates[ii]+'.fits'

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

stop
coadd_fits_maps, fits_list=fname, fileout=coadd_dir+'coadd_'+file_name+'_50mJy.fits'


; clean up
spawn, 'rm '+fname
end





;;;
; NAME: script_0228
; PURPOSE:
;   Investigate map overlap
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/28/2012: (KTS) Created
;;;

;...................................................................
; Make nmiss arrays by field
pro sss;, field_idx

;------------------------
; Setup
;------------------------
map_dir = '/data/kstory/projects/lps12/maps/20120224/'
outdir = '/data/kstory/projects/lps12/masks/apod_0228/'

field_idx = 6


;------------------------
; Get field information
;------------------------
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name

;------------------------
; get the idf fits files
;------------------------
date_list = get_lps12_runlist(field_idx)
nfiles_allidfs = n_elements(date_list)

files  = strarr(nfiles_allidfs) ; input IDF's
has_map = intarr(nfiles_allidfs) + 1

for ii=0, nfiles_allidfs-1 do begin
    if (ii mod 50 eq 0) then print, 'finding files: ' + strtrim( string(ii), 2) + '/' + strtrim( string(nfiles_allidfs), 2)
    files[ii] = map_dir+field_name+'/map_'+field_dir_name+'_150_'+date_list[ii]+'.fits'
    
    ; check for missing files
    if ~file_test(files[ii]) then begin
        print, "Missing file: ", files[ii]
        has_map[ii] = 0
    endif

endfor

; Skip missing maps
files = files[where(has_map eq 1)]
nfiles = n_elements(files)

;------------------------
; Get nmiss for this observation
;------------------------
info = expand_fits_struct(read_spt_fits(files[0])) ; get map shape
nmiss = 0*info.map.map

;for ii=0, nfiles-1 do begin
for ii=0, 9 do begin
    d = expand_fits_struct(read_spt_fits(files[ii])) ; get map shape
    weight = d.weight.map
    weight /=max(weight)

    ; find pixels that have zero obs
    wh_0 = where(weight eq 0, nwh, complement=wh_1)

    nmiss[wh_0] += 1

endfor

; Save nmiss
outname = outdir+'nmiss_'+field_name+'.sav'
save, nmiss, filename=outname
stop    
end


;...................................................................
; Investigate nmiss in observations
pro script_0228

; ; Get field information
; field_arr = lps12_fieldstruct()
; fst = field_arr[field_idx]
; field_name = fst.name
; field_dir_name = fst.dir_name

; r5h30dec-55 [0]
;myfile = '/data/kstory/projects/lps12/maps/20120224/coadds/coadd_ra5h30dec-55_2008_50mJy.fits'a
;myfile = '/data/kstory/projects/lps12/maps/20120224/ra5h30dec-55_2008/map_ra5h30dec-55_150_20080529_091936.fits'
;myfile = '/data/kstory/projects/lps12/maps/20120224/ra5h30dec-55_2008/map_ra5h30dec-55_150_20080528_042504.fits'

; ra23h30dec-55_2008 [1]
;myfile = '/data/kstory/projects/lps12/maps/20120224/coadds/coadd_ra23h30dec-55_2008_50mJy.fits'
myfile = '/data/kstory/projects/lps12/maps/20120224/ra23h30dec-55_2008/map_ra23h30dec-55_150_20080527_070222.fits'
myfile = '/data/kstory/projects/lps12/maps/20120224/ra23h30dec-55_2008/map_ra23h30dec-55_150_20080527_073643.fits'
;myfile = '/data/kstory/projects/lps12/maps/20120224/ra23h30dec-55_2008/map_ra23h30dec-55_150_20080527_153901.fits'
;myfile = '/data/kstory/projects/lps12/maps/20120224/ra23h30dec-55_2008/map_ra23h30dec-55_150_20080610_014537.fits'

d0full = expand_fits_struct(read_spt_fits(myfile))
d0full_ss = d0full.weight.map
d0full_ss /=max(d0full_ss)

; find pixels that have zero obs
wh_0 = where(d0full_ss eq 0, nwh, complement=wh_1)

; set to either 0 or 1
tmp = d0full_ss
tmp[wh_1] = 1

tmp = -1*tmp + 1 ; invert
tv_spt_map, tmp, /norms, scale=1, /forcesize


stop
end





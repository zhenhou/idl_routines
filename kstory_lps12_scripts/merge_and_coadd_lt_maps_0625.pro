;;;
; NAME: merge_and_coadd_lt_maps.pro
; PURPOSE:
;   Merge lead-trail pairs, and coadd multiple pairs in LT fields.
;   The result will be the maps that we intend to use for xspec
;
; INPUTS:
;   field_idx,       field index from lps12_fieldstruct()
;   ncoadd,          number of lt-pairs to coadd into xspec maps
;
; OUTPUTS:
;
; NOTES:
; 1) Current map directory is /data/kstory/projects/lps12/maps/20120305/
; 2) Resulting coadded maps are stored in /data/kstory/projects/lps12/maps/20120305/*field*_lt/
;
; MODIFICATION HISTORY:
;  03/01/2012: (KTS) Created
;  03/01/2012: (KTS) Change output dir to 20120305
;  03/10/2012: (KTS) Change output map name to '4Xpair_map_...', add skip_missing option
;  06/25/2012: (KTS) Re-write to use runlists and maps from 20120620
;;;


;...................................................................
; Merge and coadd maps
pro merge_and_coadd_lt_maps, field_idx, ncoadd=ncoadd, skip_missing=skip_missing
compile_opt IDL2, HIDDEN

; Parse command line arguments
if ~keyword_set(ncoadd) then ncoadd=1

; Don't use this procedure for ra23h30dec-55_2008
if field_idx eq 1 then begin
    print, 'To make coadded lt maps from ra23h30dec-55_2008, use merge_and_coadd_lt_maps_ra23h30dec-55_2008.pro'
    return
endif

;------------------------
; Get field information
;------------------------
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name

;------------------------
; get the lead-trail lists
;------------------------
lt_list_save = '/data/kstory/projects/lps12/runlists/lt/list_'+field_name+'_150.sav'
if keyword_set(skip_missing) then lt_list_save = '/data/kstory/projects/lps12/runlists/lt/list_skip_missing_'+field_name+'_150.sav'
restore, lt_list_save
nobs = n_elements(lead_list)

n1 = nobs - (nobs mod ncoadd)
nset = n1 / ncoadd
for ii=0, nset - 1 do begin

    ; Get a set of indicies to use in lead_set, trail_set
    wh = ii*ncoadd + indgen(ncoadd)
    print, "*** Coadding ", wh[0], " / ", n1

    ; Set the date of the filename to the first date
    fdate = extract_date_from_filename(lead_list[wh[0]])

    map_dir = '/data/kstory/projects/lps12/maps/20120305/'
    s_ncoadd = strtrim(string(2*ncoadd), 2)
    coadd_name = map_dir+field_name+'_lt/'+s_ncoadd+'Xpair_map_'+field_dir_name+'_150_'+fdate+'.fits'

    ; Coadd the maps
    filelist = [lead_list[wh], trail_list[wh]]
    coadd_fits_maps, filelist=filelist, fileout=coadd_name
endfor

; Print the maps that were dropped because there weren't sets of them
ndrop = nobs - n1
if ndrop gt 0 then begin
    print, "----- Dropped maps : "
    for jj=0, ndrop-1 do begin
        print, lead_list[n1+jj]
        print, trail_list[n1+jj]
    endfor
endif


end


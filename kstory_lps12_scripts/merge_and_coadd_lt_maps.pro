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
; if field_idx eq 1 then begin
;     print, 'To make coadded lt maps from ra23h30dec-55_2008, use merge_and_coadd_lt_maps_ra23h30dec-55_2008.pro'
;     return
; endif

;------------------------
; Get field information
;------------------------
f = lps12_fieldstruct()
fst = f[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name
runlist = '/home/kstory/lps12/runlists/runlist_preaz_'+field_name+'.txt'

;------------------------
; get the lead-trail lists
;------------------------
; get ncol, nrow
obslist = (read_ascii(runlist)).field1
case (size(obslist))[0] of
    1: begin
        ; should never get here if this is a LT field
        print, 'MERGE_AND_COADD_LT_MAPS: your runlist should have multiple columns, but has 1.  Returning.'
        RETURN
    end
    2: begin
        nrow=n_elements(obslist[0,*])
        ncol=n_elements(obslist[*,0])
    end
endcase

; Read the runlist of dates
case ncol of 
    2: begin
        readcol,/silent,runlist,a,b,format='a,a'
        date_list = strarr(ncol, nrow)
        date_list[0,*] = a
        date_list[1,*] = b
        ;date_list = date_list[sort(date_list)]
    end
    4: begin
        readcol,/silent,runlist,a,b,c,d,format='a,a,a,a'
        date_list = strarr(ncol, nrow)
        date_list[0,*] = a
        date_list[1,*] = b
        date_list[2,*] = c
        date_list[3,*] = d
        ;date_list = date_list[sort(date_list)]
    end
endcase

ncoadd = nrow

; Get runlists of observation map names
obs_maps  = get_lps12_runlist(field_idx, /obs_maps, /preaz)
obs_dates = get_lps12_runlist(field_idx, /obs_dates, /preaz)

for ii=0, ncoadd - 1 do begin

    ; Get a set of indicies to use in lead_set, trail_set
    
    wh = intarr(ncol)
    for icol=0, ncol-1 do begin
        wh[icol] = where( obs_dates eq date_list[icol,ii])
    endfor
    print, "*** Coadding ", ii, " / ", ncoadd

    ; Set the date of the filename to the first date
    fdate = date_list[0,ii]

    map_dir = '/data/kstory/projects/lps12/maps/20120620/'
    coadd_name = map_dir+field_name+'_lt/lt_map_'+field_dir_name+'_150_'+fdate+'.fits'

    ; Coadd the maps
    filelist = [obs_maps[wh]]
    coadd_fits_maps, filelist=filelist, fileout=coadd_name

endfor

end


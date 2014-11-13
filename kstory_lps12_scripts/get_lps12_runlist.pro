;;;
; NAME: get_lps12_runlist.pro
; PURPOSE:
;   Get a list of dates from lps12 runlists
;   You must choose one of the return formats.
;
; INPUTS:
;   field_idx,       field index from lps12_fieldstruct()
;   idfs,            return a list of low-pass filtered, downsampled idf's
;   obs_maps,        return a list of single observation map files
;   obs_dates,       return a list of single observation dates
;   xspec_maps,      return a list of xspec map files
;   xspec_dates,     return a list of xspec dates
;   cut_name,        use this to get runlists for things like
;                    'azrms_90'.  This option only works in conjunction with /xspec_maps
;
; OPTIONAL KEYWORDS:
;   preaz,           Add the /preaz keyword to any of the above formats to get the runlist before the azrms_95 cut was applied. 
;
; OUTPUTS:
;   date_list,       a strarr of the dates from the runlist for this field
;
; NOTES:
; 1) Current map directory is lps12_fieldstruct(), xspec_map_dir
; 2) Current idf directory is from lps12_fieldstruct()
;
; MODIFICATION HISTORY:
;  02/23/2012: (KTS) Created
;  03/01/2012: (KTS) Add /maps and /idfs options
;  03/02/2012: (KTS) Add /obs_maps and /xspec_maps options, force return format specification
;  03/16/2012: (KTS) Add /full_idfs
;  05/31/2012: (KTS) Add cut_name option to get azrms_90 runlists
;  06/12/2012: (KTS) Change completely to use new runlists (with azrms_95 cut), add FUNCTION get_date_list 
;  06/12/2012: (KTS) Add /preaz_dates option to get the runlist before the pre-az cut.
;  06/23/2012: (KTS) Change to /preaz option, which can be applied to any of the other return types
;;;


;-------------------------------------
; Utility function to get list of dates from runlist with columns
FUNCTION get_date_list, idx, runlist=runlist

f = lps12_fieldstruct()
if (n_elements(runlist) eq 0) then runlist = f[idx].runlist

; get ncol, nrow
obslist = (read_ascii(runlist)).field1
case (size(obslist))[0] of
    1: begin
        nrow=n_elements(obslist[0,*])
        ncol = 1
    end
    2: begin
        nrow=n_elements(obslist[0,*])
        ncol=n_elements(obslist[*,0])
    end
endcase

; Read the runlist of dates
case ncol of 
    1: begin
        readcol,/silent,runlist,date_list,format='a'
    end
    2: begin
        readcol,/silent,runlist,a,b,format='a,a'
        date_list = [a,b]
        date_list = date_list[sort(date_list)]
    end
    4: begin
        readcol,/silent,runlist,a,b,c,d,format='a,a,a,a'
        date_list = [a,b,c,d]
        date_list = date_list[sort(date_list)]
    end
endcase

ndates_all = n_elements(date_list)

return, date_list
END



;...................................................................
; Main Function
FUNCTION get_lps12_runlist, field_idx, idfs=idfs, full_idfs=full_idfs, obs_maps=obs_maps, obs_dates=obs_dates, $
                            xspec_maps=xspec_maps, xspec_dates=xspec_dates, preaz=preaz,$
                            cut_name=cut_name

; Error checking
if ( ~keyword_set(idfs) and $
     ~keyword_set(full_idfs) and $
     ~keyword_set(obs_maps) and $
     ~keyword_set(obs_dates) and $
     ~keyword_set(xspec_maps) and $
     ~keyword_set(xspec_dates) ) then begin
    print, "Error in GET_LPS12_RUNLIST: You must specify the return format.  Return -1."
    return, -1
endif

;------------------------
; Setup
;------------------------

f = lps12_fieldstruct()
fst = f[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name
runlist_dir = '/home/kstory/lps12/runlists/'

;------------------------
; If requested, get list of xspec maps
;  - Get dates from the runlist, then look for those files in the directory
;------------------------
if keyword_set(xspec_maps) then begin

    runlist = fst.runlist
    map_dir = fst.xspec_map_dir

    ; get preaz runlist if requested
    if keyword_set(preaz) then runlist = runlist_dir+'runlist_preaz_'+field_name+'.txt'

    ; add extra cuts if requested
    if n_elements(cut_name) ne 0 then begin
        runlist = '/home/kstory/lps12/runlists/runlist_'+cut_name+'_'+field_name+'.txt'
        print, 'GET_LPS12_RUNLIST: using cut_name = '+cut_name+', runlist = '+runlist
    endif

    ; Read only the first column from the runlist
    readcol,/silent,runlist,date_list,format='a'
    ndates_all = n_elements(date_list)

    files  = strarr(ndates_all) ; map names
    has_map = intarr(ndates_all) + 1

    ; grab all maps in the directory
    spawn, 'ls ' + map_dir+'*.fits', all_files
    dates = extract_date_from_filename(all_files)
    for ii=0, ndates_all-1 do begin

        ; check if a map is in our date list
        wh = where(dates eq date_list[ii], nwh)

        if nwh gt 0 then begin
            files[ii] = all_files[wh]
        ; check for missing files
        endif else begin
            print, "Missing file: ", date_list[ii]
            has_map[ii] = 0
        endelse

    endfor
    
    ; Skip missing maps
    files = files[where(has_map eq 1)]
    nfiles = n_elements(files)

    return, files
endif


;------------------------
; If requested, get list of observation maps
;  - Get dates from the runlist, then look for those files in the directory
;------------------------
if keyword_set(obs_maps) then begin

    runlist = fst.runlist
    map_dir = fst.obs_map_dir

    ; get preaz runlist if requested
    if keyword_set(preaz) then runlist = runlist_dir+'runlist_preaz_'+field_name+'.txt'

    ; Read the runlist of dates
    date_list = get_date_list(field_idx, runlist=runlist)
    ndates_all = n_elements(date_list)

    files  = strarr(ndates_all) ; map names
    has_map = intarr(ndates_all) + 1

    for ii=0, ndates_all-1 do begin
        files[ii] = map_dir+'/map_'+field_dir_name+'_150_'+date_list[ii]+'.fits'
        
        ; check for missing files
        if ~file_test(files[ii]) then begin
            print, "Missing file: ", files[ii]
            has_map[ii] = 0
        endif
    endfor
    
    ; Skip missing maps
    files = files[where(has_map eq 1)]
    nfiles = n_elements(files)
    return, files
endif


;------------------------
; If requested, get list of obs idfs
;  - Get dates from the runlist, then look for those files in the directory
;------------------------
if keyword_set(idfs) then begin

    idf_dir = fst.idf_lpfds_dirs
    runlist = fst.runlist

    ; get preaz runlist if requested
    if keyword_set(preaz) then runlist = runlist_dir+'runlist_preaz_'+field_name+'.txt'

    ; Read the runlist of dates
    date_list = get_date_list(field_idx, runlist=runlist)
    ndates_all = n_elements(date_list)

    files  = strarr(ndates_all) ; map names
    has_map = intarr(ndates_all) + 1
        
    for ii=0, ndates_all-1 do begin
        files[ii] = idf_dir+'field_scan_150_'+date_list[ii]+'.fits'
        
        ; check for missing files
        if ~file_test(files[ii]) then begin
            print, "Missing file: ", files[ii]
            print, ii
            has_map[ii] = 0
        endif
    endfor

    ; Skip missing maps
    files = files[where(has_map eq 1)]
    nfiles = n_elements(files)
    return, files
endif


;------------------------
; If requested, get list of obs full-idfs, 
; i.e. not down-sampled
;  - Get dates from the runlist, then look for those files in the directory
;------------------------
if keyword_set(full_idfs) then begin

    idf_dir = fst.idf_dirs
    runlist = fst.runlist

    ; get preaz runlist if requested
    if keyword_set(preaz) then runlist = runlist_dir+'runlist_preaz_'+field_name+'.txt'

    ; Read the runlist of dates
    date_list = get_date_list(field_idx, runlist=runlist)
    ndates_all = n_elements(date_list)

    files  = strarr(ndates_all) ; map names
    has_map = intarr(ndates_all) + 1
        
    for ii=0, ndates_all-1 do begin
        files[ii] = idf_dir+'field_scan_150_'+date_list[ii]+'.fits'
        
        ; check for missing files
        if ~file_test(files[ii]) then begin
            print, "Missing file: ", files[ii]
            has_map[ii] = 0
        endif
    endfor

    ; Skip missing maps
    files = files[where(has_map eq 1)]
    nfiles = n_elements(files)
    
    return, files
endif


;------------------------
; If requested, return list of observation dates
;------------------------
if keyword_set(obs_dates) then begin

    runlist = fst.runlist
    
    ; get preaz runlist if requested
    if keyword_set(preaz) then runlist = runlist_dir+'runlist_preaz_'+field_name+'.txt'

    ; Read the runlist of dates
    date_list = get_date_list(field_idx, runlist=runlist)
    ndates_all = n_elements(date_list)

    return, date_list
endif


;------------------------
; If requested, return list of xspec dates
;------------------------
if keyword_set(xspec_dates) then begin

    runlist = fst.runlist
    
    ; get preaz runlist if requested
    if keyword_set(preaz) then runlist = runlist_dir+'runlist_preaz_'+field_name+'.txt'

    ; Read the first column from the runlist
    readcol,/silent,runlist,date_list,format='a'

    return, date_list
endif

;------------------------
; If requested, return list of dates before the azrms_95 cut
;------------------------
; if keyword_set(preaz_dates) then begin

;     runlist = '/home/kstory/lps12/runlists/runlist_preaz_'+fst.name+'.txt'
;     date_list = get_date_list(field_idx, runlist=runlist)

;     return, date_list
; endif

print, "GET_LPS12_RUNLIST: How did we get here? Return -1."
return, -1
end



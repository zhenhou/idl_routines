;;;
; NAME: print_cut_summary
; PURPOSE:
;   print a summary of what my cuts are doing
;
; OUTPUT: Command line
;
; MODIFICATION HISTORY:
;  09/02/2011: (KTS) Created
;  09/13/2011: (KTS) Modify to work with lps12_fieldstruct
;  06/12/2012: (KTS) Modify to work with get_lps12_runlist
;;;

;...................................................................
; Main function
;
pro print_cut_summary

f = lps12_fieldstruct()


field_list = indgen(20)
;field_list = [0,1,2,5,6,7,10,11,12,13,14,15,16]
;field_list = indgen(5)+12

for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = f[idx].name
    map_dir = f[idx].autoprocessed_map_dirs

    ;--------------------
    ; Get the number of observations that passed runlist cuts
    ;--------------------

    ;restore, '/data/kstory/projects/lps12/runlists/sav_files/runlist_lps12_'+name+'.sav'
    dates_preaz = get_lps12_runlist(ii, /preaz_dates)
    npreaz = n_elements(dates_preaz)

    dates_final = get_lps12_runlist(ii, /obs_dates)
    nDatesFinal = n_elements(dates_final)

    ;--------------------
    ; Get the number of auto-processed maps
    ;--------------------
    files = file_search(map_dir + '/map_*150_*.fits',count=nfiles)
    nmaps_auto = n_elements(files)

    ;final_runlist = '/data/kstory/projects/lps12/runlists/runlist_lps12_'+name+'.txt'
    ;nDatesFinal = FILE_LINES(final_runlist)


    ;--------------------
    ; Get the number of observations made
    ;--------------------

    nobs = n_elements(read_runlist('/data/sptdat/autotools/runlists/scans_' + f[idx].dir_name + '.txt'))

    ; Special case for ra23h30dec-55, 2008 or 2010
    if ((idx eq 1) or (idx eq 2)) then begin
        nobs_2008 = 0
        nobs_2010 = 0
        r = read_runlist('/data/sptdat/autotools/runlists/scans_' + f[idx].dir_name + '.txt')
        for jj=0, n_elements(r) - 1 do begin
            if r[jj].START_MJD lt date_string_to_mjd('1-Jan-2010:00:00:00') then nobs_2008 = nobs_2008 + 1 $
            else nobs_2010 = nobs_2010 + 1
        endfor
        if (idx eq 1) then nobs = nobs_2008 else nobs = nobs_2010
    endif


;     print, "*** " + name + " ["+strtrim(string(idx),2)+"]:     # Maps that pass cuts = ", npreaz
;     print, "    (# obs that pass cuts) / (# obs that have maps made) = ", float(npreaz) / float(nfiles)
;     print, "    (# obs that pass cuts) / (# obs in auto-gen runlist) = ", float(npreaz) / float(nobs)
    print, name + " ["+strtrim(string(idx),2)+"]: ", npreaz, ', ', float(npreaz) / float(nfiles), ', ', float(npreaz) / float(nobs)


;     print, "    # obs that have maps made  =      ", long(nfiles)
;     print, "    # obs in auto-gen runlist  =      ", long(nobs)
;     print, "    # of maps in final runlists:      ", long(nDatesFinal)

;     print, "    Number of dates in sav files:     --- ", long(n_elements(dates))
;     print, "    Number of medwt in sav files:     --- ", long(n_elements(medwts))
;     print, "    Number of rms in sav files:       --- ", long(n_elements(rms_1am_midmap))
;     print, "    Number of tweights in sav files:  --- ", long(n_elements(tweight))
;     print, ""
    ;stop
endfor

stop
end


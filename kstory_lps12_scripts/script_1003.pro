;;;
; NAME: script_1003
; PURPOSE:
;   print a cut summary
;
; MODIFICATION HISTORY:
;  10/03/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro script_1003

;restore, '/data/kstory/projects/lps12/scripts/lps12_fieldnames.sav'
field_arr_ = lps12_fieldstruct()


;field_list = [0,1,indgen(13)+3]
;field_list = [0,1,2,5,6,7,10,11,12,13,14,15,16]
field_list = indgen(16)

for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name
    map_dirs = field_arr_[idx].map_dirs

    ;restore, '/data/kstory/projects/lps12/runlists/tweight_sav/runlist_tweight_'+field_name+'.sav'
    restore, '/data/kstory/projects/lps12/runlists/sav_files/runlist_lps12_'+name+'.sav'

    nMapFiles = 0.
    for jj=0, n_elements(map_dirs) - 1 do begin
        files = file_search(map_dirs[jj] + '/map_*150_*.fits',count=tmp)
        nMapFiles = nMapFiles + tmp
    endfor

    final_runlist = '/data/kstory/projects/lps12/runlists/runlist_lps12_'+name+'.txt'
    nDatesFinal = FILE_LINES(final_runlist)

    ; Number of observations made
    nobs = n_elements(read_runlist('/data/sptdat/autotools/runlists/scans_' + field_arr_[idx].dir_name + '.txt'))

    ; Special case for ra23h30dec-55, 2008 or 2010
    if ((idx eq 1) or (idx eq 2)) then begin
        nobs_2008 = 0
        nobs_2010 = 0
        r = read_runlist('/data/sptdat/autotools/runlists/scans_' + field_arr_[idx].dir_name + '.txt')
        for jj=0, n_elements(r) - 1 do begin
            if r[jj].START_MJD lt date_string_to_mjd('1-Jan-2010:00:00:00') then nobs_2008 = nobs_2008 + 1 $
            else nobs_2010 = nobs_2010 + 1
        endfor
        if (idx eq 1) then nobs = nobs_2008 else nobs = nobs_2010
    endif


    print, "*** " + name + ":    (# obs that have maps made) / (# obs in auto-gen runlist)  = ", float(nMapFiles)/float(nobs)
    print, "    # obs in auto-gen runlist = ", float(nobs)
    print, "    # obs that have maps made = ", float(nMapFiles)

;     print, "    # obs that have maps made  =      ", long(nMapFiles)
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


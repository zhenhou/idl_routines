;;;
; NAME: check_idf_progress.pro
; PURPOSE:
;   Check the progress of the idf Down-sampling job
;
; NOTES:
;
; MODIFICATION HISTORY:
;  04/22/2012: (KTS) Created
;;;

;...................................................................
; Check progress of current idf job
PRO check_idf_progress, stopit=stopit
compile_opt IDL2, HIDDEN

;--------------------
; setup
;--------------------
f = lps12_fieldstruct()

not_done=[-1]

for ii=0, 19 do begin

    fst = f[ii]
    idf_dir = '/home/rkeisler/lowellfits/'+fst.name+'/'

    rdates = get_lps12_runlist(ii, /obs_dates)
    nr = n_elements(rdates)

    files = file_search(idf_dir+'field_scan_150_*.fits')
    nf = n_elements(files)

    ; compare against runlist
    tmp = intarr(nr)
    for jj=0, nr-1 do begin
        tmp[jj] = file_test(idf_dir+'field_scan_150_'+rdates[jj]+'.fits')
    endfor

    whg = where(tmp ne 0, nwhg, complement=whb)

    gfiles = nwhg gt 0 ? files[ where(tmp ne 0)] : 0
    bfiles = whb[0] ne -1 ? files[ where(tmp eq 0)] : 0

    print, strtrim(string(ii),2), ' ' + fst.name + ': nf / nr = ', strtrim(string(nf),2) + '/'+ strtrim(string(nr),2)
    print, '          : ncomplete / nr = ', strtrim(string(nwhg),2) , '/', strtrim(string(nr),2)
    if (nwhg ge nr) then print, '**** Finished field idx: '+strtrim(string(ii),2)
    print, ' '

    ; make array of unfinished
    if (nwhg lt nr) then not_done = [not_done, ii]
endfor

not_done = not_done[1:*]
print, not_done

for jj=0, n_elements(not_done)-1 do begin
    print, 'runbatch_'+strtrim(string(not_done[jj]),2)+'.tcsh'
endfor

if keyword_set(stopit) then stop
END

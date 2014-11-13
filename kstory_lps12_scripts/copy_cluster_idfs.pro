;;;
; Procedure to make scripts to copy cluster idfs
; 05/01/2012 (KTS) Created
;;;

PRO copy_cluster_idfs

; directories
script_dir = '/home/kstory/lps12/scripts/cluster_cp_scripts/'
rk_dir = '/home/rkeisler/lowellfits/'
cluster_dir = '/sptcloud/data/rkeisler/lowellfits/'

;rk_dir = '/data/kstory/projects/lps12/scratch/tmp/' ; DEBUGGING

f = lps12_fieldstruct()

for ii=0, 19 do begin
;for ii=13, 13 do begin
    
    script_fname = script_dir + 'copy_cluster_idfs_'+strtrim(string(ii),2)+'.tcsh'
    field = f[ii].name

    ; get list of files on the cluster
    spawn, 'ls '+cluster_dir+field+'/', fcl
    ncl = n_elements(fcl)
    ;d_cl = extract_date_from_filename(fcl)

    ; get list of files in keisler
    spawn, 'ls '+rk_dir+f[ii].name+'/', frk
    nrk = n_elements(frk)
    ;d_rk = extract_date_from_filename(frk)

    ; Make a list of files not in rk directories
    need = intarr(n_elements(fcl))
    for jj=0, n_elements(fcl) -1 do begin
        need[jj] = where( frk eq fcl[jj], nwh) eq -1 ? 1 : 0
    endfor
    whneed = where(need eq 1, nwh)

    if nwh eq 0 then continue

    ; write the file
    print, ii, ": ", script_fname

    get_lun, lun1
    openw,lun1,script_fname
    printf,lun1, '#!/bin/tcsh'
    for jj=0, n_elements(whneed)-1 do begin
        printf,lun1, 'scp ' + cluster_dir+field +'/'+ fcl[whneed[jj]] + ' ' + rk_dir+field +'/' + fcl[whneed[jj]]
    endfor
    close, lun1
    free_lun,lun1
    spawn, 'chmod ug+x '+script_fname

endfor    

stop

END

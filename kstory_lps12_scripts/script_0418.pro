;;;
; NAME: script_0417
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) coadd jacklr_*.dat files
; 2) coadd run2 maps for ra1hdec-60
;
; MODIFICATION HISTORY:
;  04/17/2012: (KTS) Created
;;;

;...................................................................
; plot kmaps
PRO sss, $
         stopit=stopit
compile_opt IDL2, HIDDEN

;--------------------
; setup
;--------------------
f = lps12_fieldstruct()


for ii=0, 19 do begin
;ii=5
    if ii eq 9 then ii =10

    fst = f[ii]
    idf_dir = '/home/rkeisler/lowellfits/'+fst.name+'/'
    spawn, 'ls '+idf_dir+'field_scan_150_*.fits', newlist
    nnew = n_elements(newlist)

    runlist = get_lps12_runlist(ii, /obs_maps)
    nold = n_elements(runlist)

    print, ii, ' ' + fst.name + ': nnew / nold = ', nnew , '/', nold
endfor
stop

END

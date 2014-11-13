;;;
; NAME: coadd_all_lps12_sims.pro
; PURPOSE:
;   wrapper to coadd sims from all fields.
;
; INPUTS:
;   jacklr,           coadd jacklr's
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  06/04/2012: (KTS) Created from /home/cr/code/spt/lowell/coadd_all_lowell_sims.pro
;;;

pro coadd_all_lps12_sims,jacklr=jacklr, lmax4500=lmax4500, stub=stub

dojack = keyword_set(jacklr) ? 1 : 0

for i=0,20 do begin

    f = lps12_fieldstruct()
    field = f[i].name
    print, 'COADD_ALL_LPS12_SIMS: field = ' + field
    
    ; call on old sims
    if keyword_set(lmax4500) then begin
        if n_elements(stub) ne 0 then begin
            coadd_lps12_sims,field, /lmax4500, jacklr=dojack, stub=stub
        endif else begin
            coadd_lps12_sims,field, /lmax4500, jacklr=dojack
        endelse

    ; call on new sims
    endif else begin
        if n_elements(stub) ne 0 then begin
            coadd_lps12_sims,field, jacklr=dojack, stub=stub
        endif else begin
            coadd_lps12_sims,field, jacklr=dojack
        endelse
    endelse

endfor


end

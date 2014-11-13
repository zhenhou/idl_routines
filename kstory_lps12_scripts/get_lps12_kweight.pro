;;;
; NAME: get_lps12_kweight.pro
; PURPOSE:
;   Return a 2-D kweight
;
; INPUTS: field_idx,      index in lps12_fieldstruct()
;
; OUTPUTS:
;   kweight,         2-D weight mask
;
; NOTES:
;   1) give option to use box cut of all kx < ell=500
;
; MODIFICATION HISTORY:
;  03/07/2012: (KTS) Created
;  05/08/2012: (KTS) use actual twod_kweight
;  06/05/2012: (KTS) Re-arrange setup into box-mask section.
;;;


;...................................................................
; Return the kweight
function get_lps12_kweight, field_idx, box_kmask=box_kmask

f    = lps12_fieldstruct()

; Get the kweight
if ~keyword_set(box_kmask) then begin
    restore, '/home/kstory/lps12/twod_kweights/weight_2d_'+f[field_idx].name+'.sav'
    kweight = weight_2d

; if asked, just use a box mask
endif else begin
    reso_arcmin=1.0
    reso=reso_arcmin/60.*!dtor
    
    info = get_lps12_fieldinfo(field_idx)
    npix = info.npix
    nbig = info.nbig
    
    ell_cut = 300.

    dell    = 2*!pi/(nbig*reso)
    jj      = fix(ell_cut/dell)
    kweight = fltarr(nbig,nbig)+1
    
    kweight[0:jj,*]=0
    kweight[nbig-jj:nbig-1,*]=0
endelse

return, kweight
end


;;;
; NAME: make_field_weights
; PURPOSE:
;   Make a sav file with the relative weights of the different fields
;
; INPUTS: 
;   mask_dir,      directory with masks (apod + ptsrc)
;   ptsrc_stub,    specify which ptsrc mask was used (for name only)
;
; OUTPUTS:
;   saves a file in mask_dir/field_weights.sav
;
; NOTES:
; 1) Based purely on weighted area
; 2) Cunnently not taking into account that field 1 and 2 are both ra23h30dec-55
;
; MODIFICATION HISTORY:
;  03/21/2012: (KTS) Created
;  05/31/2012: (KTS) Change to use get_lps12_mask.pro
;;;


PRO make_field_weights, mask_dir=mask_dir
compile_opt IDL2, HIDDEN

field_arr_ = lps12_fieldstruct()
nfields = 20
field_weights = fltarr(nfields) ; return array

; get the weights
for ii=0, nfields-1 do begin
    field_name = field_arr_[ii].name

    if keyword_set(mask_dir) then begin
        restore, mask_dir+'mask_'+field_name+'.sav'
    endif else begin
        mask_dir = '/home/kstory/lps12/masks/masks_50mJy/'
        mask = get_lps12_mask(ii)
    endelse

    field_weights[ii] = total(mask)

endfor

; re-normalize
total_weight = total(field_weights)
field_weights /= total_weight

; save
savname = mask_dir+'field_weights.sav'
save, field_weights, total_weight, filename=savname
stop
END

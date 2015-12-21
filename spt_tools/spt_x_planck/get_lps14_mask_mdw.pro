;;;
; NAME: get_lps12_mask.pro
; PURPOSE:
;   Easily return the mask by field
;
; INPUTS:
;   field_idx,       field index from lps12_fieldstruct()
;   padded,          return the padded mask
;
; OUTPUTS:
;   The mask, padded to 4320
;
; NOTES:
;
; MODIFICATION HISTORY:
;  04/26/2012: (KTS) Created
;  04/30/2012: (KTS) add 'padded' option, change defaut to un-padded
;;;


;...................................................................
; main function
FUNCTION get_lps14_mask_mdw, field_idx, type=type, padded=padded

;if n_elements(type) eq 0 then type='50mJy'

home = getenv('HOME')

f = mk_field_arr()
mask_dir = home+'/data/spt_data/masks/lps14/' ;/home/kstory/lps12/masks/masks_'+type+'/'
restore, mask_dir+'mask_'+f[field_idx].name+'.sav'

ret_mask = keyword_set(padded) ? mask_padded : mask

RETURN, ret_mask
END

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
FUNCTION get_lps12_mask_zh, field_idx, padded=padded

host = getenv('HOSTNAME')
if host eq 'spt' then home='/home/hou/'
if host eq  'midway' then home='/home/zhenhou/'

f = lps12_fieldstruct()
mask_dir = home+'data/projects/spt_x_planck/lps12/masks_50mJy/'
restore, mask_dir+'mask_'+f[field_idx].name+'.sav'

ret_mask = keyword_set(padded) ? mask_padded : mask

RETURN, ret_mask
END

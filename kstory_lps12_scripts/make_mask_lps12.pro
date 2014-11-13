;;;
; NAME: make_masks_lps12.pro
; PURPOSE:
;   Make final masks from apodization masks
;
; INPUTS:
;   field_idx,            index from lps12_fieldstruct
;   apod_dir,             directory with apodization masks
;   ptsrc_dir,            directory with ptsrc masks
;
; NOTES:
;
; TODO:
;
; MODIFICATION HISTORY:
;  03/13/2012: (KTS) Created
;  03/21/2012: (KTS) finalized from make_final_mask.pro
;  04/30/2013: (KTS) add outdir as optional argument
;;;


;...................................................................
; make final apod*ptsrc_mask mask for one field
pro make_mask_lps12_single_field, field_idx, apod_mask_file, ptsrc_mask_file, outdir
compile_opt IDL2, HIDDEN

field_arr = lps12_fieldstruct()
fst = field_arr[field_idx] & field_name = fst.name
info = get_lps12_fieldinfo(field_idx)
nbig = info.nbig

print, "make final mask for field ", field_name

; get the apod and ptsrc mask
restore, apod_mask_file
restore, ptsrc_mask_file
mask = apod*ptsrc_mask

; make padded mask
mask_padded = dblarr(nbig, nbig)
nx = (size(mask))[1]
ny = (size(mask))[2]
mask_padded[0:nx-1, 0:ny-1] = mask

; save the output
print, '  Save mask: '+outdir+'mask_'+field_name+'.sav'
stop
save, mask, mask_padded, filename=outdir+'mask_'+field_name+'.sav'

end




;...................................................................
; Main Function: Wrapper for single_field.
pro make_mask_lps12, $
                     field_idx = field_idx, $
                     apod_dir  = apod_dir, $
                     ptsrc_dir = ptsrc_dir, $
                     outdir    = outdir
compile_opt IDL2, HIDDEN
                           
                           

;------------------------
; Setup
;------------------------
field_arr = lps12_fieldstruct()

if ~keyword_set(apod_dir) then apod_dir = '/home/kstory/lps12/masks/apod/'
if ~keyword_set(ptsrc_dir) then ptsrc_dir = '/home/kstory/lps12/masks/ptsrc/'

ptsrc_stub = '50mJy'
if ~keyword_set(outdir) then outdir = '/home/kstory/lps12/masks/masks_'+ptsrc_stub+'/'

;------------------------
; If we only want to make one apod mask
;------------------------
if keyword_set(field_idx) then begin

    fst = field_arr[field_idx] & field_name = fst.name
    ; cop out and hard-copy the extensions
    apod_mask_file = apod_dir + 'apod_'+field_name + '_60_0.0500_30.sav'
    ptsrc_mask_file = ptsrc_dir + 'ptsrc_mask_'+field_name + '_'+ptsrc_stub+'_proj5.sav'

    ; if the apod mask does not exist, make it
    if ~file_test(apod_mask_file) then begin
        print, 'MAKE_MASK_LPS12, re-make apod mask for field ', strtrim(string(field_idx),2)
        make_apod_mask_lps12, field_idx
    endif

    make_mask_lps12_single_field, field_idx, apod_mask_file, ptsrc_mask_file, outdir
    return
endif


;------------------------
; otherwise make all masks    
;------------------------
for field_idx=0, 19 do begin

    fst = field_arr[field_idx] & field_name = fst.name
    ; cop out and hard-copy the extensions
    apod_mask_file = apod_dir + 'apod_'+field_name + '_60_0.0500_30.sav'
    ptsrc_mask_file = ptsrc_dir + 'ptsrc_mask_'+field_name + '_'+ptsrc_stub+'_proj5.sav'

    ; if the apod mask does not exist, make it
    if ~file_test(apod_mask_file) then begin
        print, 'MAKE_MASK_LPS12, re-make apod mask for field ', strtrim(string(field_idx),2)
        make_apod_mask_lps12, field_idx
    endif

    make_mask_lps12_single_field, field_idx, apod_mask_file, ptsrc_mask_file, outdir
endfor

end



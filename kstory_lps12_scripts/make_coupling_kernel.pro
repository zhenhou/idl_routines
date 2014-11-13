;;;
; NAME: make_coupling_kernel
; PURPOSE:
;   Make coupling kernel for each field from the masks
;
; INPUTS:
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  03/21/2012: (KTS) Created
;  03/26/2012: (KTS) Change mask to be padded
;  04/03/2012: (KTS) Set maxell to 4650, what ryan used.
;  04/05/2012: (KTS) Fix winsize
;  05/01/2012: (KTS) use get_lps12_mask() function
;;;


PRO make_coupling_kernel, field_idx, $
                          mask_padded=mask_padded, $ ; apod mask (padded to 4320)
                          savfile=savfile, $   ; output sav filename
                          test_pad=test_pad, $ ; testing only
                          oversamp=oversamp

compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
info = get_lps12_fieldinfo(0)
npix_pad = info.nbig
if keyword_set(test_pad) then npix_pad = test_pad ; debugging only
if ~keyword_set(oversamp) then oversamp=8

field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name
if ~keyword_set(savfile) then $
  savfile = '/data/kstory/projects/lps12/masks/masks_50mJy/kern_'+field_name+'.sav'

print, 'Analyzing field ', field_idx, ', name = ', field_name

; inputs
reso_arcmin = 1.0
maxell = 4650
if ~keyword_set(mask_padded) then mask_padded = get_lps12_mask(field_idx, /padded)

; calculate the coupling kernel
kernel=coupling_kernel(mask_padded, reso_arcmin, maxell=maxell, /changevar, $
                       interp=1000, oversamp=oversamp, /cheby, ellkern=ellkern, $
                       curlyw=curlyw)

; re-normalize (not sure what this does exactly)
winsize = npix_pad

reso=double(reso_arcmin)/60*!dtor
kernsize=(size(kernel))[1]
u=(dindgen(kernsize)+0.5)#replicate(1./(reso*winsize)^4, kernsize)
kernel*=u    

; save file
print, 'save file to: ' + savfile
save, kernel, mask_padded, reso_arcmin, maxell, ellkern, curlyw, filename=savfile
END

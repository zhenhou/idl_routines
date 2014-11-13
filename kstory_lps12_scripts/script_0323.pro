;;;
; NAME: script_0323
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) make tweight jack runlists
; 2) run tweight jackknife
; 3) run moon jackknife
; 4) run sun jackknife
; 5) test script for making coupling kernels
; 6) make all coupling kernels
;
; MODIFICATION HISTORY:
;  03/23/2012: (KTS) Created
;;;

;...................................................................
; make tweight jackknife runlists
PRO run_make_tweight_jack_list
compile_opt IDL2, HIDDEN

;for ii=0, nfields-1 do begin
for ii=1, 19 do begin

    make_tweight_jack_list, ii

endfor
END

;...................................................................
; run tweight jackknife
PRO run_tweight_jack
compile_opt IDL2, HIDDEN

;for ii=0, nfields-1 do begin
for ii=1, 19 do begin

    lps12_jack, ii, 'tweight', savdir='/home/kstory/lps12/jacks/'

endfor
END

;...................................................................
; run moon jackknife
PRO run_moon_jack
compile_opt IDL2, HIDDEN

;for ii=0, nfields-1 do begin
for ii=1, 19 do begin

    lps12_jack, ii, 'moon', savdir='/home/kstory/lps12/jacks/'

endfor
END

;...................................................................
; run moon jackknife
PRO run_sun_jack
compile_opt IDL2, HIDDEN

restore, '/data/kstory/projects/lps12/masks/masks_50mJy/sun_weights.sav'
idx = where(p_used gt 0.25)
nfields = n_elements(idx)

for ii=0, nfields-1 do begin

    field_idx = idx[ii]
    lps12_jack, field_idx, 'sun', savdir='/home/kstory/lps12/jacks/'

endfor
END


;...................................................................
; Test making kernels
PRO mk_kern, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name
savfile = '/data/kstory/projects/lps12/masks/masks_50mJy/kern_'+field_name+'.sav'

print, 'Analyzing field ', field_idx, ', name = ', field_name

; inputs
reso_arcmin = 1.0
maxell = 3000
maskfile = '/home/kstory/lps12/masks/masks_50mJy/mask_'+field_name+'.sav'
restore, maskfile

kernel=coupling_kernel(mask, reso_arcmin, maxell=maxell, /changevar, $
                       interp=1000, oversamp=8, /cheby, ellkern=ellkern, $
                       curlyw=curlyw)

stop
; save file
print, 'save file to: ' + savfile
save, kernel, mask, reso_arcmin, maxell, ellkern, curlyw, filename=savfile
stop
END


;...................................................................
; make coupling kernels
PRO run_kern
compile_opt IDL2, HIDDEN
for ii=0, 19 do begin
    print, 'Make coupling kernel for field idx ', ii
    make_coupling_kernel, ii
endfor
END

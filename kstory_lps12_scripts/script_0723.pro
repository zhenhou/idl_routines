;;;
; NAME: script_0723
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) beam test
;
; MODIFICATION HISTORY:
;  07/23/2012: (KTS) Created
;;;

PRO end2end_run08b1_idx4
idx = 4
tbegin = systime(0, /seconds)
print, 'run end2end_08b1, field ', idx
t0 = systime(0, /seconds)

try_end2end_08b1, idx, /use_kweight, /resume

tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'

END



;; Make kernel files from run_08
PRO save_kerns_08
f = lps12_fieldstruct()
edir = '/home/kstory/lps12/end2end/'
kerndir = '/home/kstory/lps12/masks/masks_50mJy/'

for i=0, 19 do begin
    restore, edir+'end_'+f[i].name+'_08_kweight.sav'
    save, ellkern, kernel, mask, mask_padded, maxell, reso_arcmin, $
      filename=kerndir+'kern_'+f[i].name+'.sav'
endfor

END


;; Beam test run_08b1
PRO b_08b1
year = 2009
beamfiles_orig = get_lps12_beams(year, l_beam, bl_beam)
beamfiles = '/home/kstory/lps12/scripts/beamtest/beam_2009_run08b1.txt'

; Test changing the beam by 5% at ell=3000
fl = -2.128*10^(-5.) *(l_beam - 650) + 1
bl_new = bl_beam*fl
get_lun, lun1
openw,lun1,beamfiles
for i=0, n_elements(l_beam)-1 do begin
    printf, lun1, l_beam[i], bl_new[i]
endfor
close, lun1
free_lun, lun1
stop
END


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


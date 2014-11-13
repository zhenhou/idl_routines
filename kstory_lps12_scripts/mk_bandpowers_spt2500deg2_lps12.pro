;;;
; Script for making SZ template for the s12 data products website
;
; Notes:
;  1) Make txt file of bandpwoers
;;;

;;;;;;;;;;
;
; Make template for SZ sources
;
;;;;;;;;;;
PRO mk_bandpowers_spt2500deg2_lps12

script_fname = '/home/kstory/lps12/webpage/bandpowers_spt2500deg2_lps12.txt'

restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120828_170101_kweight.sav'

istart=9
istop=55
dl_all_lps12 = dl_all[istart:istop] * 1d12
diag_nobeam_lps12 = diag_nobeam[istart:istop] * 1d12
l = l[istart:istop]
nbin = n_elements(l)

print, 'write to: '+script_fname

get_lun, lun1
openw,lun1,script_fname
printf,lun1,'# SPT Bandpowers from Story et al 2012.'
printf,lun1,'# These bandpowers are provided for plotting purposes only.'
printf,lun1,'# The bandpower errors do not include uncertainty in the'
printf,lun1,'# SPT beam and calibration.'
printf,lun1,'# Version 1.0, Created 21Nov2012.'
printf,lun1,'#'
printf,lun1,'# ell_center, D_ell (uK2), sigma(D_ell) (uK2)'
printf,lun1,'#'

for i=0, nbin-1 do begin
    printf,lun1, round(l[i]),'   ', dl_all_lps12[i],'   ', diag_nobeam_lps12[i], format = '(i6, A3, f6.1, A3, f4.1)'
endfor
close, lun1
free_lun,lun1

stop
END


;;;
; NAME: script12_0828
; PURPOSE:
;   General script for today (2012)
;
; NOTES:
; 1) Compare spectra
;;;


;;;;;;;;;;;;;;;;;
;
; Compare spectra
;
;;;;;;;;;;;;;;;;;


PRO compare_spec, idx
idx=5
edir='/home/kstory/lps12/end2end/'
f=lps12_fieldstruct()

restore, edir+'run_08/combined_spectrum_20120717_174249_kweight.sav'
dl08 = dl_all*1d12
cov08 = cov_all
diag08 = diag*1d12
diag_nobeam08 = diag_nobeam*1d12

restore, edir+'run_09/combined_spectrum_20120828_170101_kweight.sav'
dl09 = dl_all*1d12
cov09 = cov_all
diag09 = diag*1d12
diag_nobeam09 = diag_nobeam*1d12

vec=indgen(47)+8

fdir='/home/kstory/public_html/notebook/spt_lps12/'
; PLOTS
window, 1
plot, l[vec], dl08[vec], /yl, ytitle='Dl [uK^2]'
errplot, l[vec], (dl08-diag_nobeam08)[vec],(dl08+diag_nobeam08)[vec]
oplot, l[vec], dl09[vec], color=!red
errplot, l[vec], (dl09-diag_nobeam09)[vec],(dl09+diag_nobeam09)[vec],color=!red
legend,['ps run_08','ps run_09'],colors=[!black,!red],linestyle=[0,0]
err=tvread(/png,/nodialog,filename=fdir+'pspec_comp_0828')


;window, 2
plot, l[vec], ((dl09-dl08)/dl09)[vec], yr=[-0.01,0], ytitle='(dl09-dl08)/dl09', title='Combined Spectrum'
err=tvread(/png,/nodialog,filename=fdir+'pspec_frac_0828')

stop
END


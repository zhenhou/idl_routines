;;;
; NAME: script13_0501.pro
;
; NOTES:
;  1) Check ptscor 2500d spectrum
;;;

;;;;;;;;;;;;;;;;;;;
;
; 
;
;;;;;;;;;;;;;;;;;;;
PRO comp_ptscor
restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_ApJ2012.sav'
d0 = dl_all

restore, '/home/kstory/lps12/end2end/run_ptscor/combined_spectrum_20130501_175332_kweight.sav'
d1 = dl_all

vec = indgen(47)+9

window, 1
!x.margin=[12,3]
dl = 'D!D!12l!X!N'
plot, l[vec],((d1-d0)/d0)[vec], xtitle='!12l!X!N', $
  ytitle='('+dl+'_ptscor - '+dl+'_S12) / '+dl+'_S12', $
  title='Correct for missing point sources'
;err = tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/ptscor_ratio_0501')

window, 2
plot, l[vec], (d1-d0)[vec]*1d12, xtitle='!12l!X!N', ytitle=dl+'_ptscor - '+dl+'_S12 [uK^2]',$
  title='Difference in Power'
;err = tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/ptscor_power_0501')

window, 3
plot, l[vec], ((d1-d0)/diag_nobeam)[vec],xtitle='!12l!X!N',$
  ytitle='Shift in units of !4r!X!N'
err = tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/ptscor_sigma_0501')

stop
END

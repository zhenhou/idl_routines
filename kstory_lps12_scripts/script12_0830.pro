;;;
; NAME: script12_0830
; PURPOSE:
;   General script for today (2012)
;
; NOTES:
; 1) Check 2010-2011 only bandpowers
;;;


;;;;;;;;;;;;;;;;;;
;
; check 2010, 2011-only bandpowers
;
;;;;;;;;;;;;;;;;;;
PRO check_bp
restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120828_170101_kweight.sav'
dl_tot = dl_all*1d12
diag_tot = diag_nobeam*1d12
wf2use = wf_all

restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120830_210134_kweight_1011.sav'
dl_1011 = dl_all*1d12
diag_1011 = diag_nobeam*1d12

restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120830_154612_kweight_0809.sav'
dl_0809 = dl_all*1d12
diag_0809 = diag_nobeam*1d12

;-------------------
; get w7s12 best-fit
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', l_w7s12, dl_ML_w7s12_total

mnmx = fltarr(2)
mnmx[0] = max( [min(l_w7s12), min(l_wf)] )
mnmx[1] = min( [max(l_w7s12), max(l_wf)] )

; cut down wf
istart = (where_closest(l_wf,mnmx[0]))[0]
istop = (where_closest(l_wf,mnmx[1]))[0]
l_wf    = l_wf[istart:istop]
wf2use  = wf2use[istart:istop,*]

; cut down theory spectrum
istart = (where_closest(l_w7s12,mnmx[0]))[0]
istop = (where_closest(l_w7s12,mnmx[1]))[0]
dl_ML_w7s12_total = dl_ML_w7s12_total[istart:istop]
l_w7s12           = l_w7s12[istart:istop]

; Apply window function
dl_ML = dl_ML_w7s12_total
nl = n_elements(dl_all)
dl_th = dblarr(nl)
for i=0,nl-1 do dl_ML[i] = total(dl_ML_w7s12_total*wf2use[*,i])
;-------------------

;; plots
vec=indgen(47)+9

window, 0
plot,    l[vec], dl_tot[vec], yr=[50,3000],/yst,/yl
errplot, l[vec], (dl_tot-diag_tot)[vec], (dl_tot+diag_tot)[vec]
oplot,   l[vec], dl_0809[vec], color=!red
errplot, l[vec], (dl_0809-diag_0809)[vec], (dl_0809+diag_0809)[vec], color=!red
oplot,   l[vec], dl_1011[vec], color=!blue
errplot, l[vec], (dl_1011-diag_1011)[vec], (dl_1011+diag_1011)[vec], color=!blue

window, 1
plot,    l[vec], (dl_ML-dl_tot)[vec]
errplot, l[vec], (dl_ML-dl_tot-diag_tot)[vec], (dl_ML-dl_tot+diag_tot)[vec]

window, 2
plot,    l[vec], ((dl_ML-dl_tot)/dl_ML)[vec]
errplot, l[vec], ((dl_ML-dl_tot-diag_tot)/dl_ML)[vec], ((dl_ML-dl_tot+diag_tot)/dl_ML)[vec]

window, 3,xsize=900,ysize=400
plot, l[vec], ((dl_ML-dl_0809)/dl_ML)[vec],yr=[-0.2,0.2],/yst
oplot, l[vec], ((dl_ML-dl_tot)/dl_ML)[vec]
errplot, l[vec], ((dl_ML-dl_tot-diag_tot)/dl_ML)[vec], ((dl_ML-dl_tot+diag_tot)/dl_ML)[vec]
oplot, l[vec], ((dl_ML-dl_0809)/dl_ML)[vec], color=!red
errplot, l[vec], ((dl_ML-dl_0809-diag_0809)/dl_ML)[vec], ((dl_ML-dl_0809+diag_0809)/dl_ML)[vec],color=!red
oplot, l[vec], ((dl_ML-dl_1011)/dl_ML)[vec], color=!blue
errplot, l[vec], ((dl_ML-dl_1011-diag_1011)/dl_ML)[vec], ((dl_ML-dl_1011+diag_1011)/dl_ML)[vec],color=!blue

stop

END

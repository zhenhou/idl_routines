;;;
; NAME: plots_s12_s12_k11_compare.pro
; PURPOSE:
;   Plot the differences between K11 and S12
; 
; Notes:
;  1) plots_s12_k11s12_diff_sigma: Plot \Delta C_l/\sigma C_l
;  2) plots_s12_k11s12_diff: Plot (dl_s12 - dl_k11)/dl_s12
;
; MODIFICATION HISTORY:
;  08/27/2012: (KTS) Created
;;;


;;;;;;;;;;;;;;;;;;
; Plot \Delta C_l/\sigma C_l
; oplot_ml if you want to overplot the best-fit spectra difference
;;;;;;;;;;;;;;;;;;

PRO plots_s12_k11s12_diff_sigma, oplot_ml=oplot_ml
;--------------------
; Compare best-fit
;--------------------
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', format='d,d', ell_s12, dl_s12

restore, '/home/rkeisler/ps09/bestfit_WMAP7_SPT_lowell.sav'
ell_k11 = l_cmb
dl_k11  = dl_total

lstart = max([ell_s12[0], l_cmb[0]]) ; starting ell
lstop  = min([max(ell_s12), max(l_cmb)]) ; ending ell
nl = lstop-lstart+1

istart = where(ell_s12 eq lstart)
istop  = where(ell_s12 eq lstop)
;dl_s12 = dl_s12[istart:istop]
xs = dl_s12[istart:istop]

istart = where(ell_k11 eq lstart)
istop  = where(ell_k11 eq lstop)
;dl_k11 = dl_k11[istart:istop]
xk = dl_k11[istart:istop]

ell_th = indgen(nl) + lstart

;--------------------
; Compare bandpowers
;--------------------
restore, '/data23/kstory/lps12/end2end/run_09/combined_spectrum_ApJ2012.sav'
vec = indgen(47)+9
dl_s12 = dl_all[vec] *1d12
ell_s12 = l[vec]
err_s12 = diag[vec]*1d12
lwf_s12 = l_wf
wf_s12  = wf_all


restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
;readcol, '/home/kstory/lps12/best_fit/dl_k11.txt',format='d,d,d',ell_k11,dl_k11,err_k11
dl_k11 = dl_all[vec] *1d12
ell_k11 = l[vec]
err_k11 = diag[vec]*1d12
lwf_k11 = l_wf
wf_k11  = wf_all

; Difference
diff = (dl_s12 - dl_k11)/err_k11
mean = mean(diff)
d_data = diff-mean ; Do I get to do this?

;--------------------
;bin the theory spectrum
;--------------------
i1 = min(ell_th)
i2 = max(ell_th)

istart = (where_closest(lwf_s12,i1))[0]
istop  = (where_closest(lwf_s12,i2))[0]
l_wf    = lwf_s12[istart:istop]
wf_s12 = wf_s12[istart:istop,*]
wf_k11 = wf_k11[istart:istop,*]

istart = (where_closest(ell_th,i1))[0]
istop  = (where_closest(ell_th,i2))[0]
ell_th = ell_th[istart:istop]
xs     = xs[istart:istop,*]
xk     = xk[istart:istop,*]

nl = n_elements(dl_all)
xsb = dblarr(nl)
xkb = dblarr(nl)
for i=0, nl-1 do xsb[i] = total(xs * wf_s12[*,i])
for i=0, nl-1 do xkb[i] = total(xk * wf_k11[*,i])

; difference
diff = (xsb[vec] - xkb[vec])/err_k11
mean = mean(diff)
d_ml = diff-mean


;--------------------
; Plots
;--------------------
plot, ell_s12, d_data, psym=5,title='Difference between K11 and S12',xtitle='!12l!X!N',$
  ytitle='(!12D!Dl!U!3S12!N - !12D!Dl!U!3K11 !X!N) / !4r!12D!Dl!U!3K11!X!N', chars=1.8,thick=2,charthick=1,$
  yr=[-2,2],/yst
oplot, ell_s12, d_data, psym=5,thick=2,color=!red
oplot, ell_s12, d_ml,thick=2
legend,['Bandpowers','ML spectra'],color=[!red,!black],linestyle=[0,0],psym=[5,0]

print, "Save figure"
pdir = '/home/kstory/lps12/scripts/plotting/'
err=tvread(/png,/nodialog,filename=pdir+'s12_k11_diff_sigma')
spawn, 'cp /home/kstory/lps12/scripts/plotting/s12_k11_diff_sigma.png ~/public_html/notebook/spt_lps12/.'


stop
END



;;;;;;;;;;;;;;;;;;
;
; Plot difference: (dl_s12 - dl_k11)/dl_s12
;
;;;;;;;;;;;;;;;;;;
PRO plots_s12_k11s12_diff
;--------------------
; Compare best-fit
;--------------------
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', format='d,d', ell_s12, dl_s12

restore, '/home/rkeisler/ps09/bestfit_WMAP7_SPT_lowell.sav'
ell_k11 = l_cmb
dl_k11  = dl_total

lstart = max([ell_s12[0], l_cmb[0]]) ; starting ell
lstop  = min([max(ell_s12), max(l_cmb)]) ; ending ell
nl = lstop-lstart

istart = where(ell_s12 eq lstart)
istop  = where(ell_s12 eq lstop)
;dl_s12 = dl_s12[istart:istop]
xs = dl_s12[istart:istop]

istart = where(ell_k11 eq lstart)
istop  = where(ell_k11 eq lstop)
;dl_k11 = dl_k11[istart:istop]
xk = dl_k11[istart:istop]

ell = indgen(nl) + lstart

; difference
dd = (xs - xk)/xs


;--------------------
; Compare bandpowers
;--------------------
restore, '/data23/kstory/lps12/end2end/run_09/combined_spectrum_ApJ2012.sav'
vec = indgen(47)+9
dl_s12 = dl_all[vec] *1d12
ell_s12 = l[vec]

readcol, '/home/kstory/lps12/best_fit/dl_k11.txt',format='d,d,d',ell_k11,dl_k11,err_k11
dd_data = (dl_s12 - dl_k11)/dl_s12


ys = '!12D!Dl!X!N!Us12!X!N'
yk = '!12D!Dl!X!N!Uk11!X!N'
plot, ell_s12, dd_data, psym=5,title='Compare K11 and S12',xtitle='!12l!X!N',$
  $;ytitle='!4d!12D!Dl!X!N (!7l!6K!3!U2!X!N)', chars=1.8,thick=2,charthick=2
  $;ytitle='(!12D!Dl!X!N-s12 - !12D!Dl!X!N-k11 )/!12D!Dl!X!N-s12  (!7l!6K!3!U2!X!N)', chars=1.8,thick=2,charthick=2
  ytitle='('+ys+' - '+yk+' )/'+ys+'  (!7l!6K!3!U2!X!N)', $
  chars=2,thick=2,charthick=2
oplot, ell_s12, dd_data, psym=5,color=!red,thick=2
oplot, ell, dd,thick=2
legend,['Best-fit','Bandpowers'],color=[!black,!red],linestyle=[0,0],psym=[0,5]

print, "Save figure"
pdir = '/home/kstory/lps12/scripts/plotting/'
err=tvread(/png,/nodialog,filename=pdir+'s12_k11_diff')
spawn, 'cp /home/kstory/lps12/scripts/plotting/s12_k11_diff.png ~/public_html/notebook/spt_lps12/.'

stop
END

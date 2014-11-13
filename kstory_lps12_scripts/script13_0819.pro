;;;
; NOTES:
;  1) S12 ref response; look at differences between K11 and S12 best-fit
;;;

PRO compare_best_fit
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

stop
END


;;;;;;;;;;;;;;;;;;
; Plot \Delta C_l/\sigma C_l
;;;;;;;;;;;;;;;;;;

PRO plots_s12_k11s12_diff_sigma
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
err_s12 = diag[vec]*1d12

readcol, '/home/kstory/lps12/best_fit/dl_k11.txt',format='d,d,d',ell_k11,dl_k11,err_k11
diff = (dl_s12 - dl_k11)/err_s12
mean = mean(diff)
d2 = diff-mean ; Do I get to do this?

plot, ell_s12, diff, psym=5,title='Difference between K11 and S12',xtitle='!12l!X!N',$
  ytitle='(!12D!Dl!U!3S12!N - !12D!Dl!U!3K11 !X!N) / !4r!12D!Dl!U!3S12!X!N (!7l!6K!3!U2!X!N)', chars=1.8,thick=2,charthick=1,$
  yr=[-3,3],/yst

print, "Save figure"
pdir = '/home/kstory/lps12/scripts/plotting/'
err=tvread(/png,/nodialog,filename=pdir+'s12_k11_diff_sigma_0819')
spawn, 'cp scripts/plotting/s12_k11_diff_sigma_0819.png ~/public_html/notebook/spt_lps12/.'

stop
END



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
err=tvread(/png,/nodialog,filename=pdir+'s12_k11_diff_0819')
spawn, 'cp /home/kstory/lps12/scripts/plotting/s12_k11_diff_0819.png ~/public_html/notebook/spt_lps12/.'

stop
END

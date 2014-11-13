;;;
; NAME: script12_0825
; PURPOSE:
;   General script for today (2012)
;
; NOTES:
; 1) Plotting residuals
;;;


;;;;;;;;;;;;;;;;;;
;
; Get a calibration factor, given data, theory, and icov
;
;;;;;;;;;;;;;;;;;;
FUNCTION get_cal, dl_all, dl_th, icov
dcal = 0.0001
mincal = 1./1.3
maxcal = 1.*1.3
ncal = ceil((maxcal-mincal)/dcal)
cal = findgen(ncal)*dcal + mincal
chisq = dblarr(ncal)
for i=0,ncal-1 do begin
    dl_tmp = dl_all*cal[i]
    delta = dl_tmp-dl_th
    chisq[i] = (delta ## (icov ## delta))[0]
endfor


; Get the ratio of data to theory
chisq_tmp = chisq
min_chisq = min(chisq_tmp)
wh_min=(where(chisq_tmp eq min_chisq))[0]
best_cal = cal[wh_min]
dcal = abs(best_cal - cal[where_closest(chisq_tmp,min_chisq+1.)])
print,'full cov, BEST_CAL: ',best_cal,' +/- ',dcal


RETURN, best_cal
END




;;;;;;;;;;;;;;;;;;
;
; Make sav-file with everything I need
;
;;;;;;;;;;;;;;;;;;
PRO mk_sav

; get WMAP best-fit
readcol,'/home/kstory/lps12/best_fit/lcdm_w7_ML_dl_fromZH.txt',l_w7,dl_w7,te,ee,bb

; get w7s12 best-fit
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', l_total, dl_total
l_vec = l_total
;---------------------
; Get the foregrounds
;---------------------

; Get ML values, from params_c4_lcdm_camb_w7s12_ML.ini
poisson3000 = 20.0976
clust_3000 = 5.46899
sz3000 = 4.18801

;---------------------
; get the residual Poisson power
help,poisson3000
dl_poisson = poisson3000*((l_vec/3000.)^2.)

;---------------------
; get the power from SZ
b=read_ascii('/home/cr/paramfits/cosmomc.r11/ptsrc/dl_shaw_tsz_s10_153ghz.txt')
yy = b.field1
; match the ell range with that of the cmb
istart = (where_closest(reform(yy[0,*]),2))[0]
istop = istart + n_elements(l_vec) -1
yy = yy[*,istart:istop]
l_sz = reform(yy[0,*])
dl_sz = reform(yy[1,*])
wh_3000 = (where_closest(l_sz,3000))[0]
help, sz3000
dl_sz = dl_sz/dl_sz[wh_3000]*sz3000

wh_lstart = where( l_sz eq l_vec[0])
l_sz = l_sz[wh_lstart:*]
dl_sz = dl_sz[wh_lstart:*]

;---------------------
; get the power from spatially correlated galaxies and tSZ and kSZ
help, clust_3000
dl_flat = clust_3000*((l_vec/3000.)^(0.8))

;---------------------
; Trim w7 vector to the same size
istart = (where_closest(l_w7,2))[0]
istop = istart + n_elements(l_vec) -1
l_w7  = l_w7[istart:istop]
dl_w7 = dl_w7[istart:istop]

;---------------------
; Back-out dl_w7s12
dl_foregrounds = dl_sz + dl_poisson + dl_flat
dl_w7s12 = dl_total - dl_foregrounds

;---------------------
; Save to file
;---------------------
filename = '/home/kstory/lps12/scripts/plotting/lcdm_ML_spectra_0825.sav'
print, 'saving to file: ', filename
save,  l_vec, dl_w7, dl_w7s12, dl_sz, dl_poisson, dl_flat, dl_foregrounds, $
  filename=filename

END





;;;;;;;;;;;;;;;;;;
;
; Plot residuals
;
;;;;;;;;;;;;;;;;;;
PRO resid_1
;---------------------
; Setup
;---------------------
pdir = '/home/kstory/lps12/scripts/plotting/'
;fdir = '/home/kstory/public_html/notebook/spt_lps12/'
fdir = '/home/kstory/lps12/scripts/figs/'

; Get best-fit lines
restore, '/home/kstory/lps12/scripts/plotting/lcdm_ML_spectra_0825.sav'
dl_w7s12_total = dl_w7s12 + dl_foregrounds
dl_w7_total = dl_w7 + dl_foregrounds

; get data
data = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
restore, data

; copies for plotting down to very low ell
dl_w7s12_totalo = dl_w7s12_total
dl_w7s12o = dl_w7s12
l_veco = l_vec

; Specify which spectrum to use for calculating the best-fit calibration
dl_4cal = dl_w7s12_total
;dl_4cal = dl_w7_total


;---------------------
; Load and process lps12 spectrum 
;---------------------

istart = (where_closest(l_vec,min(l_wf)))[0]
istop = (where_closest(l_vec,max(l_wf)))[0]
dl_w7s12 = dl_w7s12[istart:istop]
dl_poisson = dl_poisson[istart:istop]
dl_flat = dl_flat[istart:istop]
dl_w7s12_total = dl_w7s12_total[istart:istop]
dl_4cal = dl_4cal[istart:istop]
l_vec = l_vec[istart:istop]

; only use a subset of the bandpowers
istart=9
istop=55
cov_all = cov_all[istart:istop,istart:istop]
cov_all_nobeam = cov_all_nobeam[istart:istop,istart:istop]
cov_all_nocal = cov_all_nobeam
diag = diag[istart:istop]
diag_nobeam = diag_nobeam[istart:istop]
dl_all = dl_all[istart:istop]
l = l[istart:istop]
wf_all = wf_all[*,istart:istop]
;stop
nl = n_elements(dl_all)
dl_th = dblarr(nl)
for i=0,nl-1 do dl_th[i] = total(dl_4cal*wf_all[*,i])

; convert the data from K to uK2
dl_all *= (1d12)
diag *= (1d12)
diag_nobeam *= (1d12)
cov_all *= (1d24)
cov_all_nobeam *= (1d24)
cov_all_nocal *= (1d24)

icov = invert(cov_all,/double)
icov_nocal = invert(cov_all_nocal,/double)
icov_nobeam = invert(cov_all_nobeam,/double)


;---------------------
; find the calibration factor which is provides the best fit between
; the data and the theory spectrum.
;---------------------

dl_th_w7s12 = dblarr(nl)
for i=0,nl-1 do dl_th_w7s12[i] = total(dl_w7s12_total*wf_all[*,i])
best_cal_w7s12 = get_cal(dl_all, dl_th_w7s12, icov_nobeam)

dl_th_w7 = dblarr(nl)
for i=0,nl-1 do dl_th_w7[i] = total(dl_w7_total*wf_all[*,i])
best_cal_w7 = get_cal(dl_all, dl_th_w7, icov_nobeam)


best_cal = best_cal_w7s12

; Get ratios
ratio_w7s12 = dl_all*best_cal_w7s12/dl_th_w7s12
err_ratio_w7s12 = diag_nobeam/dl_th_w7s12

ratio_w7 = dl_all*best_cal_w7/dl_th_w7
err_ratio_w7 = diag_nobeam/dl_th_w7

; get residuals, all with w7s12 cal
residual_w7s12 = dl_all*best_cal_w7s12 - dl_th_w7s12
residual_w7    = dl_all*best_cal_w7s12 - dl_th_w7

;---------------------
; Plot
;---------------------

window,0
color_cmb = !eggplant
color_poisson = !darkgreen
color_flat = !orange
color_total = !black
color_data = !blue
thick1 = 1
thick2 = 1
plot,[0],[0],/nodata,xtitle='!12l!X!N', $
  ytitle='D!D!12l!X!N '+textoidl('(\muK^2)'), chars=1.8, $
  color=!black, xthick=1, ythick=1, $
  xr=[0,3e3],yr=[40,8000],/yst,/xst,/yl
;oplot,l,dl_all*best_cal,color=color_data
;oplot,l,dl_all*best_cal,color=color_data,ps=1
errplot,l,dl_all*best_cal-diag_nobeam,dl_all*best_cal+diag_nobeam,color=color_data,thick=2
;pause
;stop
oplot,l_vec,dl_w7s12,color=color_cmb,thick=thick1
oplot,l_vec,dl_poisson,color=color_poisson,thick=thick1
oplot,l_vec,dl_flat,color=color_flat,thick=thick1
oplot,l_vec,dl_w7s12_total,color=color_total,thick=thick2
errplot,l,dl_all*best_cal-diag_nobeam,dl_all*best_cal+diag_nobeam,color=color_data,thick=2




ratio = ratio_w7s12
err_ratio = err_ratio_w7s12

window,1
plot,l,ratio,/yn,/nodata,xtitle='!12l!X!N',ytitle='Data/Theory',yr=[-1,1]*0.2+1,/yst
oplot,[-1,1]*9999999999.,[1,1],thick=2
oplot,l,ratio
errplot,l,(ratio-err_ratio),(ratio+err_ratio)

save, l, ratio_w7, err_ratio_w7, ratio_w7s12, err_ratio_w7s12, residual_w7, residual_w7s12, $
  dl_th_w7s12, dl_th_w7, $
  filename=pdir+'lcdm_ratios_residuals_0825.sav'
stop
END






;;;;;;;;;;;;;;;;;;
;
; compare ML lines
;
;;;;;;;;;;;;;;;;;;
PRO comp_ml
pdir = '/home/kstory/lps12/scripts/plotting/'
;fdir = '/home/kstory/public_html/notebook/spt_lps12/'
fdir = '/home/kstory/lps12/scripts/figs/'

; Get best-fit lines
restore, '/home/kstory/lps12/scripts/plotting/lcdm_ML_spectra_0825.sav'
dl_w7s12_total = dl_w7s12 + dl_foregrounds
dl_w7_total = dl_w7 + dl_foregrounds

; plot spectra
window, 1
plot, l_vec, dl_w7s12_total, /ylog
oplot, l_vec, dl_w7_total, color=!red



;----------------
; Ratios plot
;----------------

window, 2, ysize=500,xsize=800

;-------------------
; plot difference
yr = [-40,40]
plot, l_vec, (dl_w7s12_total - dl_w7_total), thick=2, $
    ytitle="Ignore units", xtitle='ell'
legend,['ML_w7s12 - ML_w7','Ratio_w7', 'ratio_w7s12'],color=[!black,!red,!blue],linestyle=[0,0,0],thick=[2,1,1],$
  pos=[2000,35]
;stop
;-------------------
; plot ML lines
dd = alog(dl_w7s12) ; baseline
dd1 = dd - min(dd)
mxmn = maxmin(dd1)
; dvec = (dd1 / mxmn[0]) * (yr[1]-yr[0]) + yr[0]

dd_w7 = alog(dl_w7) & dd_w7 = dd_w7 - min(dd)
dvec_w7 = (dd_w7 / mxmn[0]) * (yr[1]-yr[0]) + yr[0]

dd_w7s12 = alog(dl_w7s12) & dd_w7s12 = dd_w7s12 - min(dd)
dvec_w7s12 = (dd_w7s12 / mxmn[0]) * (yr[1]-yr[0]) + yr[0]

oplot, l_vec, dvec_w7, color=!red
oplot, l_vec, dvec_w7s12, linestyle=2


;-------------------
; plot ML Ratios
restore, pdir+'lcdm_ratios_residuals_0825.sav'

;; need to scale this to -40,40
yr1=[0.5, 1.5]
rr_w7 = (ratio_w7-1) * (yr[1]-yr[0])/(yr1[1]-yr1[0])
err_rr_w7 = (err_ratio_w7) * (yr[1]-yr[0])/(yr1[1]-yr1[0])

rr_w7s12 = (ratio_w7s12-1) * (yr[1]-yr[0])/(yr1[1]-yr1[0])
err_rr_w7s12 = (err_ratio_w7s12) * (yr[1]-yr[0])/(yr1[1]-yr1[0])

oplot,l,rr_w7, color=!red, thick=2
errplot,l,(rr_w7-err_rr_w7),(rr_w7+err_rr_w7), color=!red

oplot,l,rr_w7s12, color=!blue
errplot,l,(rr_w7s12-err_rr_w7s12),(rr_w7s12+err_rr_w7s12), color=!blue

;stop

;-------------------
; Save figure
;-------------------
ff = 'lcdm_ratios_0825'
print, 'save figure to: ', fdir+ff+'.png'
err=tvread(/png,/nodialog,filename=fdir+ff)
spawn, 'cp '+fdir+ff+'.png /home/kstory/public_html/notebook/spt_lps12/.'





;----------------
; Residuals plot
;----------------

window, 3, ysize=500,xsize=800

;-------------------
; plot difference
yr = [-600,600]
plot, l, dl_th_w7s12 - dl_th_w7, thick=2, $
    ytitle="uK^2", xtitle='ell'
legend,['ML_w7s12 - ML_w7','Resid_w7', 'Resid_w7s12'],color=[!black,!red,!blue],linestyle=[0,0,0],thick=[2,1,1],$
  pos=[2000,400]
;stop
;-------------------
; plot ML lines
dd = alog(dl_w7s12) ; baseline
dd1 = dd - min(dd)
mxmn = maxmin(dd1)
; dvec = (dd1 / mxmn[0]) * (yr[1]-yr[0]) + yr[0]

dd_w7 = alog(dl_w7) & dd_w7 = dd_w7 - min(dd)
dvec_w7 = (dd_w7 / mxmn[0]) * (yr[1]-yr[0]) + yr[0]

dd_w7s12 = alog(dl_w7s12) & dd_w7s12 = dd_w7s12 - min(dd)
dvec_w7s12 = (dd_w7s12 / mxmn[0]) * (yr[1]-yr[0]) + yr[0]

oplot, l_vec, dvec_w7, color=!red
oplot, l_vec, dvec_w7s12, linestyle=2


;-------------------
; plot ML Residuals
restore, pdir+'lcdm_ratios_residuals_0825.sav'

;; need to scale this to -40,40
oplot,l,residual_w7, color=!red, thick=2
oplot,l,residual_w7s12, color=!blue

stop

;-------------------
; Save figure
;-------------------
ff = 'lcdm_residuals_0825'
print, 'save figure to: ', fdir+ff+'.png'
err=tvread(/png,/nodialog,filename=fdir+ff)
spawn, 'cp '+fdir+ff+'.png /home/kstory/public_html/notebook/spt_lps12/.'

stop
END




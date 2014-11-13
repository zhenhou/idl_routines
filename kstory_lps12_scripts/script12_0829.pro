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
readcol,'/home/kstory/lps12/best_fit/lcdm_w7_ML_dl_fromZH.txt',l_w7,dl_ML_w7,te,ee,bb

; get w7s12 best-fit
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', l_w7s12, dl_ML_w7s12_total
l_vec = l_w7s12

; cut down to l_vec
istart = (where_closest(l_w7, l_vec[0]))[0]
istop = istart + n_elements(l_vec) -1
l_w7 = l_w7[istart:istop]
dl_ML_w7 = dl_ML_w7[istart:istop]


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
istart = (where_closest(reform(yy[0,*]),l_vec[0]))[0]
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
; Back-out dl_w7s12
dl_foregrounds = dl_sz + dl_poisson + dl_flat
dl_ML_w7s12 = dl_ML_w7s12_total - dl_foregrounds

;---------------------
; Make best-fit spectra
;---------------------
dl_ML_w7_total = dl_ML_w7 + dl_foregrounds


; get data
data = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
restore, data

; copies for plotting down to very low ell
dl_ML_w7s12_totalo = dl_ML_w7s12_total
dl_ML_w7s12o = dl_ML_w7s12
l_veco = l_vec


;---------------------
; Load and process lps12 spectrum 
;---------------------

; Set the ell-range for comparison; this means cutting down both l_vec and l_wf
mnmx = fltarr(2)
mnmx[0] = max( [min(l_vec), min(l_wf)] )
mnmx[1] = min( [max(l_vec), max(l_wf)] )

istart = (where_closest(l_wf,mnmx[0]))[0]
istop = (where_closest(l_wf,mnmx[1]))[0]
l_wf    = l_wf[istart:istop]
wf_all  = wf_all[istart:istop,*]

; cut everything down to the ell-range of the window-function
istart = (where_closest(l_vec,mnmx[0]))[0]
istop = (where_closest(l_vec,mnmx[1]))[0]

dl_ML_w7s12       = dl_ML_w7s12[istart:istop]
dl_ML_w7s12_total = dl_ML_w7s12_total[istart:istop]
dl_ML_w7          = dl_ML_w7[istart:istop]
dl_ML_w7_total    = dl_ML_w7_total[istart:istop]
dl_poisson        = dl_poisson[istart:istop]
dl_flat           = dl_flat[istart:istop]
dl_sz             = dl_sz[istart:istop]
dl_foregrounds    = dl_foregrounds[istart:istop]
l_vec             = l_vec[istart:istop]

; Specify which spectrum to use for calculating the best-fit calibration
dl_4cal = dl_ML_w7s12_total
;dl_4cal = dl_w7_total


; only use a subset of the bandpowers
istart=9
istop=55
cov_all        = cov_all[istart:istop,istart:istop]
cov_all_nobeam = cov_all_nobeam[istart:istop,istart:istop]
cov_all_nocal  = cov_all_nobeam
diag           = diag[istart:istop]
diag_nobeam    = diag_nobeam[istart:istop]
dl_all         = dl_all[istart:istop]
l              = l[istart:istop]
wf_all         = wf_all[*,istart:istop]

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
; Make the binned theory spectra
;---------------------

dl_ML_w7s12_binned = dblarr(nl)
for i=0,nl-1 do dl_ML_w7s12_binned[i] = total(dl_ML_w7s12_total*wf_all[*,i])

dl_ML_w7_binned = dblarr(nl)
for i=0,nl-1 do dl_ML_w7_binned[i] = total(dl_ML_w7_total*wf_all[*,i])

;---------------------
; find the calibration factor which is provides the best fit between
; the data and the theory spectrum.
;---------------------

best_cal_w7s12 = get_cal(dl_all, dl_ML_w7s12_binned, icov_nobeam)
best_cal_w7 = get_cal(dl_all, dl_ML_w7_binned, icov_nobeam)

best_cal = best_cal_w7s12

; Get ratios
ratio_w7s12 = dl_all/dl_ML_w7s12_binned
err_ratio_w7s12 = diag_nobeam/dl_ML_w7s12_binned

ratio_w7 = dl_all/dl_ML_w7_binned
err_ratio_w7 = diag_nobeam/dl_ML_w7_binned

; get residuals, all with w7s12 cal
residual_w7s12 = dl_all - dl_ML_w7s12_binned
residual_w7    = dl_all - dl_ML_w7_binned




;---------------------
; Save to file
;---------------------
filename = '/home/kstory/lps12/scripts/plotting/lcdm_ML_spectra_0829.sav'
print, 'saving to file: ', filename
save,  l_vec, dl_ML_w7,dl_ML_w7_total, dl_ML_w7s12,dl_ML_w7s12_total, dl_sz, dl_poisson, dl_flat, dl_foregrounds, $
  filename=filename


save, l, ratio_w7, err_ratio_w7, ratio_w7s12, err_ratio_w7s12, residual_w7, residual_w7s12, diag_nobeam,$
  dl_ML_w7s12_binned, dl_ML_w7_binned, $
  dl_all, best_cal_w7s12, best_cal_w7, $
  filename='/home/kstory/lps12/scripts/plotting/lcdm_ratios_residuals_0829.sav'

stop
END





;;;;;;;;;;;;;;;;;;
;
; compare ML lines
;
;;;;;;;;;;;;;;;;;;
PRO plot_resid_ml
pdir = '/home/kstory/lps12/scripts/plotting/'
;fdir = '/home/kstory/public_html/notebook/spt_lps12/'
fdir = '/home/kstory/lps12/scripts/figs/'

; Get best-fit lines
restore, '/home/kstory/lps12/scripts/plotting/lcdm_ML_spectra_0829.sav'

; get residuals and ratios
restore, '/home/kstory/lps12/scripts/plotting/lcdm_ratios_residuals_0829.sav'

; plot spectra
window, 1
plot, l_vec, dl_ML_w7s12_total, /ylog
oplot, l_vec, dl_ML_w7_total, color=!red


;----------------
; Residuals plot, no error bars
;----------------

window, 2, ysize=500,xsize=800
;-------------------
; plot difference
yr = [-60,100]
;yr = [-100,150]
xr = [0,3200]
plot, l_vec, dl_ML_w7s12_total - dl_ML_w7_total, thick=2, $
    ytitle="Dl [uK^2]", xtitle='ell', xr=xr,yr=yr,/xst,/yst
oplot, l, dl_ML_w7s12_binned - dl_ML_w7_binned, thick=2, color=!darkgreen, psym=4
oplot, [xr],[0,0],linestyle=2
legend,['ML_w7s12 - ML_w7','Resid_w7', 'Resid_w7s12'],color=[!black,!red,!blue],psym=[-4,-3,-3],thick=[2,2,2],$
  pos=[2000,35]

;-------------------
; plot ML Residuals
oplot,l,residual_w7, color=!red, thick=2,psym=-4
oplot,l,residual_w7s12, color=!blue, thick=2,psym=-4

; Save figure
ff = 'lcdm_residuals_0829'
print, 'save figure to: ', fdir+ff+'.png'
;err=tvread(/png,/nodialog,filename=fdir+ff)
;spawn, 'cp '+fdir+ff+'.png /home/kstory/public_html/notebook/spt_lps12/.'

; Save figure
ff = 'lcdm_residuals_0829'
print, 'save figure to: ', fdir+ff+'.png'
;err=tvread(/png,/nodialog,filename=fdir+ff)
;spawn, 'cp '+fdir+ff+'.png /home/kstory/public_html/notebook/spt_lps12/.'


;----------------
; Residuals plot, ERROR BARS
;----------------

window, 3, ysize=500,xsize=800
;-------------------
; plot difference
;yr = [-60,100]
yr = [-100,180]
xr = [0,3200]
plot, l_vec, dl_ML_w7s12_total - dl_ML_w7_total, thick=2, $
    ytitle="Dl [uK^2]", xtitle='ell', xr=xr,yr=yr,/xst,/yst
oplot, l, dl_ML_w7s12_binned - dl_ML_w7_binned, thick=2, color=!darkgreen, psym=4
oplot, [xr],[0,0],linestyle=2
legend,['ML_w7s12 - ML_w7','Resid_w7', 'Resid_w7s12'],color=[!black,!red,!blue],psym=[-4,-3,-3],thick=[2,2,2],$
  pos=[2000,35]

;-------------------
; plot ML Residuals
oplot,l,residual_w7, color=!red, thick=2,psym=-4
oplot,l,residual_w7s12, color=!blue, thick=2,psym=-4

; Save figure
ff = 'lcdm_residuals_0829'
print, 'save figure to: ', fdir+ff+'.png'
;err=tvread(/png,/nodialog,filename=fdir+ff)
;spawn, 'cp '+fdir+ff+'.png /home/kstory/public_html/notebook/spt_lps12/.'

errplot, l, (residual_w7-diag_nobeam), (residual_w7+diag_nobeam),color=!red, thick=2
errplot, l, (residual_w7s12-diag_nobeam), (residual_w7s12+diag_nobeam),color=!blue, thick=2

;-------------------
; Save figure
;-------------------
ff = 'lcdm_residuals_err_0829'
print, 'save figure to: ', fdir+ff+'.png'
;err=tvread(/png,/nodialog,filename=fdir+ff)
;spawn, 'cp '+fdir+ff+'.png /home/kstory/public_html/notebook/spt_lps12/.'


stop
END



;;;;;;;;;;;;;;;;;;
;
; compare ML lines
;
;;;;;;;;;;;;;;;;;;
PRO plot_frac_resid_ml
pdir = '/home/kstory/lps12/scripts/plotting/'
;fdir = '/home/kstory/public_html/notebook/spt_lps12/'
fdir = '/home/kstory/lps12/scripts/figs/'

; Get best-fit lines
restore, '/home/kstory/lps12/scripts/plotting/lcdm_ML_spectra_0829.sav'
; get residuals and ratios
restore, '/home/kstory/lps12/scripts/plotting/lcdm_ratios_residuals_0829.sav'


;-------------------
; plot ML Residuals
;-------------------

window, 2, ysize=500,xsize=800
xr = [0,3200]
plot,l,residual_w7/dl_ML_w7_binned, $
  ytitle="residual/dl_ML", xtitle='ell', xr=xr,/xst
oplot,l,residual_w7/dl_ML_w7_binned, color=!red, thick=2,psym=-4
oplot,l,residual_w7s12/dl_ML_w7s12_binned, color=!blue, thick=2,psym=-4
oplot, [xr],[0,0],linestyle=2
legend,['residual_w7/dl_ML_w7_total','residual_w7s12/dl_ML_w7s12_total'],color=[!red,!blue],psym=[-3,-3],thick=[2,2],$
  pos=[1500,0.06]

; Save figure
ff = 'lcdm_frac_resid_0829'
print, 'save figure to: ', fdir+ff+'.png'
err=tvread(/png,/nodialog,filename=fdir+ff)
spawn, 'cp '+fdir+ff+'.png /home/kstory/public_html/notebook/spt_lps12/.'


stop
END







; ;----------------
; ; Ratios plot
; ;----------------

; window, 3, ysize=500,xsize=800

; ;-------------------
; ; plot difference
; yr = [-40,40]
; plot, l_vec, dl_ML_w7s12_total - dl_ML_w7_total, thick=2, $
;     ytitle="Dl [uK^2]", xtitle='ell'
; oplot, l, dl_ML_w7s12_binned - dl_ML_w7_binned, thick=2, color=!red, psym=4
; legend,['ML_w7s12 - ML_w7','Ratio_w7', 'ratio_w7s12'],color=[!black,!red,!blue],linestyle=[0,0,0],thick=[2,1,1],$
;   pos=[2000,35]
; ;stop

; ;-------------------
; ; plot ML Ratios

; ; ;; need to scale this to -40,40
; ; yr1=[0.5, 1.5]
; ; rr_w7 = (ratio_w7-1) * (yr[1]-yr[0])/(yr1[1]-yr1[0])
; ; err_rr_w7 = (err_ratio_w7) * (yr[1]-yr[0])/(yr1[1]-yr1[0])

; ; rr_w7s12 = (ratio_w7s12-1) * (yr[1]-yr[0])/(yr1[1]-yr1[0])
; ; err_rr_w7s12 = (err_ratio_w7s12) * (yr[1]-yr[0])/(yr1[1]-yr1[0])

; oplot,l,ratio_w7, color=!red, thick=2
; errplot,l,(ratio_w7-err_ratio_w7),(ratio_w7+err_ratio_w7), color=!red

; oplot,l,ratio_w7s12, color=!blue
; errplot,l,(ratio_w7s12-err_ratio_w7s12),(ratio_w7s12+err_ratio_w7s12), color=!blue

; ;stop

; ;-------------------
; ; Save figure
; ;-------------------
; ff = 'lcdm_ratios_0829'
; print, 'save figure to: ', fdir+ff+'.png'
; ;err=tvread(/png,/nodialog,filename=fdir+ff)
; ;spawn, 'cp '+fdir+ff+'.png /home/kstory/public_html/notebook/spt_lps12/.'

; stop
; END




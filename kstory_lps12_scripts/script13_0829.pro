;;;
; NOTES:
;  1) exploring ns constraints
;  2) mk_chisq: calculate chisq of k11 data to k11_ML and S12_ML spectra
;  3) make best-fit k11 spectrum and s12 spectrum
;;;

;--------------------
; Compare errors on ns
;--------------------
PRO pp
ns_w = 0.966850
err_w = 0.0135156

ns_15 = 0.961549
err_15 =  0.0112503

ns_30 = 0.968078
err_30 = 0.0104815


print, "Change in err w->15: ",  1-(err_15/err_w)
print, "Change in err w->30: ",  1-(err_30/err_w)

stop
END


PRO check_w7_ns
cdir = '/data23/kstory/lps12/chains/'

; WMAP-only
dir_w7 = cdir+'c1_lcdm_pico_w7_newKp/chains/'
files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7_newKp*.txt')
pname_w7 = dir_w7+'c1_lcdm_pico_w7_newKp.paramnames'

plot_like1dname,files_w7,pname_w7,'ns',subsamp=1000,nskip=1000
stop

END


;--------------------
; Calculate the chisq of the K11 data to the S12 best-fit
;--------------------
PRO mk_chisq
;--------------------

;--------------------
; data from K11
vec = indgen(47)+9
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
;readcol, '/home/kstory/lps12/best_fit/dl_k11.txt',format='d,d,d',ell_k11,dl_k11,err_k11
dl_k11 = dl_all*1d12
l_bin = l
err_k11 = diag*1d12

istart=9
istop=55
cov_k11 = cov_all[istart:istop,istart:istop] * 1d24
icov = invert(cov_k11,/double)


;--------------------
;Get ML, bin theory spectrum: S12
;--------------------
; Best-fit spectrum from s12
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', format='d,d', ellml_s12, dlml_s12

; get s12 window functions
restore, '/data23/kstory/lps12/end2end/run_09/combined_spectrum_ApJ2012.sav'
vec = indgen(47)+9
lwf_s12 = l_wf
wf_s12  = wf_all

i1 = max( [min(ellml_s12), min(lwf_s12)] )
i2 = min( [max(ellml_s12), max(lwf_s12)] )

; cut down theory spectrum
istart = (where_closest(ellml_s12,i1))[0]
istop  = (where_closest(ellml_s12,i2))[0]
ellml_s12 = ellml_s12[istart:istop]
dlml_s12 = dlml_s12[istart:istop]

; cut down wf
istart  = (where_closest(lwf_s12,i1))[0]
istop   = (where_closest(lwf_s12,i2))[0]
lwf_s12 = lwf_s12[istart:istop]
wf_s12  = wf_s12[istart:istop,*]

; bin theory
nl = n_elements(dl_all)
dlmlb_s12 = dblarr(nl)
for i=0, nl-1 do dlmlb_s12[i] = total(dlml_s12 * wf_s12[*,i])


;--------------------
;Get ML, bin theory spectrum: K11
;--------------------
; Best-fit spectrum from k11
;readcol, '/home/kstory/lps12/best_fit/dl_k11.txt',format='d,d,d',ellml_k11,dlml_k11,errml_k11
;readcol, '/home/rkeisler/ps09/wmap_k11_bestfit_r0p00_summnu0p0eV.txt',ellml_k11,dlml_k11,format='d,d'
readcol, '/home/rkeisler/ps09/k11_bestfit_total_dl_lcdm.txt',ellml_k11,dlml_k11,format='d,d'

; get k11 window functions
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
vec = indgen(47)+9
lwf_k11 = l_wf
wf_k11  = wf_all

i1 = max( [min(ellml_k11), min(lwf_k11)] )
i2 = min( [max(ellml_k11), max(lwf_k11)] )

; cut down theory spectrum
istart = (where_closest(ellml_k11,i1))[0]
istop  = (where_closest(ellml_k11,i2))[0]
ellml_k11 = ellml_k11[istart:istop]
dlml_k11 = dlml_k11[istart:istop]

; cut down wf
istart  = (where_closest(lwf_k11,i1))[0]
istop   = (where_closest(lwf_k11,i2))[0]
lwf_k11 = lwf_k11[istart:istop]
wf_k11  = wf_k11[istart:istop,*]

; bin theory
nl = n_elements(dl_all)
dlmlb_k11 = dblarr(nl)
for i=0, nl-1 do dlmlb_k11[i] = total(dlml_k11 * wf_k11[*,i])



; for plotting
dl_k11 = dl_k11[vec]
err_k11 = err_k11[vec]
l_bin = l_bin[vec]
dlmlb_k11 = dlmlb_k11[vec]
dlmlb_s12 = dlmlb_s12[vec]


;--------------------
; remove a mean from the data
;--------------------


; ;---------------------
; ; find the calibration factor which is provides the best fit between
; ; the data and the theory spectrum.
; dcal = 0.0001
; mincal = 1./1.3
; maxcal = 1.*1.3
; ncal = ceil((maxcal-mincal)/dcal)
; cal = findgen(ncal)*dcal + mincal
; chisq_s12 = dblarr(ncal)
; chisq_k11 = dblarr(ncal)
; for i=0,ncal-1 do begin
;     dl_tmp = dl_k11*cal[i]
;     delta = dl_tmp - dlmlb_s12
;     chisq_s12[i] = (delta ## (icov ## delta))[0]
;     delta = dl_tmp - dlmlb_k11
;     chisq_k11[i] = (delta ## (icov ## delta))[0]
; endfor

; ; s12 ML spectrum
; chisq_tmp = chisq_s12
; min_chisq = min(chisq_tmp)
; wh_min=(where(chisq_tmp eq min_chisq))[0]
; best_cal = cal[wh_min]
; dcal = abs(best_cal - cal[where_closest(chisq_tmp,min_chisq+1.)])
; print,'full cov, BEST_CAL: ',best_cal,' +/- ',dcal
; print, 'K11 data to S12_ML: chisq = ', min_chisq

; ; k11 ML spectrum
; chisq_tmp = chisq_k11
; min_chisq = min(chisq_tmp)
; wh_min=(where(chisq_tmp eq min_chisq))[0]
; best_cal = cal[wh_min]
; dcal = abs(best_cal - cal[where_closest(chisq_tmp,min_chisq+1.)])
; print,'full cov, BEST_CAL: ',best_cal,' +/- ',dcal
; print, 'K11 data to K11_ML: chisq = ', min_chisq


; stop

;;;;;;
; mean_dd = total( (dl_k11 - dlmlb_s12)/err_k11 ) / total(1/err_k11)
; dl_k11_unscaled = dl_k11
; dl_k11 = dl_k11 - mean_dd


; chisq
delta = dl_k11 - dlmlb_s12
chisq_s12 = ( delta ## (icov ## delta))[0]
print, "K11 data to S12_ML: chisq = ", chisq_s12

delta = dl_k11 - dlmlb_k11
chisq_k11 = ( delta ## (icov ## delta))[0]
print, "K11 data to K11_ML: chisq = ", chisq_k11

print, chisq_k11 - chisq_s12
stop
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Make best-fit k11 spectrum
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO mk_bestfit_k11



stop
END

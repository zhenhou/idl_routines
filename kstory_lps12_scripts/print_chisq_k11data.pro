;;;
; NAME: print_chisq_k11data
; PURPOSE:calculate chisq of k11 data to k11_ML and S12_ML spectra
;
; NOTES:
;
; MODIFICATION HISTORY:
;  08/29/2013: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Calculate the chisq of the K11 data to the S12 best-fit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO print_chisq_k11data
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


; Use only the bandpowers we care about
dl_k11 = dl_k11[vec]
err_k11 = err_k11[vec]
l_bin = l_bin[vec]
dlmlb_k11 = dlmlb_k11[vec]
dlmlb_s12 = dlmlb_s12[vec]



;--------------------
; chisq
;--------------------

delta = dl_k11 - dlmlb_s12
chisq_s12 = ( delta ## (icov ## delta))[0]
print, "K11 data to S12_ML: chisq = ", chisq_s12

delta = dl_k11 - dlmlb_k11
chisq_k11 = ( delta ## (icov ## delta))[0]
print, "K11 data to K11_ML: chisq = ", chisq_k11

print, chisq_k11 - chisq_s12
stop
END


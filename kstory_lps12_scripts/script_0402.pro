;;;
; NAME: script_0401
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) test_distr; prob of getting 0.0012 from 120 draws from uniform distribution
;
; MODIFICATION HISTORY:
;  04/01/2012: (KTS) Created
;;;

;...................................................................
; test distribution
PRO test_distr
compile_opt IDL2, HIDDEN

ntries = 120
nsims = 100000L
undefine, seed

nexp = 100
exp_success = intarr(nexp)
for iexp=0, nexp-1 do begin
    nsuccess = 0L

    for iv=0, nsims-1 do begin
        vals=RANDOMU(seed,[ntries])
        ; check for value
        wh = where(vals lt 0.0012, nwh)
        if nwh ne 0 then nsuccess++

        ;plot, histogram(vals, nbins=10) & pause
    endfor

    print, iexp, nsuccess, nsuccess / float(nsims)

    exp_success[iexp] = nsuccess
endfor

yvec = exp_success/float(nsims) & mnmx = minmax(yvec) & xvec=indgen(10)*(mnmx[1]-mnmx[0])/10.+mnmx[0] & plot, xvec, histogram(yvec, nbins=10)
stop
        
END


;...................................................................
; lr jack failure
PRO sss, plotit=plotit, stopit=stopit
compile_opt IDL2, HIDDEN

;field_idx = 6
f = lps12_fieldstruct()

list = indgen(14)+6
nfields = n_elements(list)

; make arrays
nbf_mean    = fltarr(nfields)
nbf_median  = fltarr(nfields)
avgW_mean   = fltarr(nfields)
avgW_median = fltarr(nfields)


; hard-code chisq for lr jack from 0401
chisq = [10.16,  1.87,  2.06,  9.34, 3.62, 5.59,16.86, 2.74, 8.75, 2.24, 8.16, 9.19, 3.13, 0.84, 3.60, 5.07, 8.80, 4.96, 3.02,14.91]
pte   = [ 0.0708, 0.8670, 0.8405, 0.0963, 0.6057, 0.3481, 0.0048, 0.7406, 0.1193, 0.8145, 0.1477, 0.1017, 0.6801, 0.9742, 0.6087, 0.4074, 0.1175, 0.4208, 0.6962, 0.0107]

for ii=0, nfields-1 do begin

    field_idx = list[ii]
    fst = f[field_idx]
    field_name = fst.name

; print nice statement
    print, '---' + field_name + '-------------------------------------------'

; get jackknife file
    jackfile = '/home/kstory/lps12/scripts/sav_files/lr_test_0402_'+field_name+'.sav'
    restore, jackfile

    ss = size(nbolos_flagged_L)
    nobs    = ss[1]
    nscans  = ss[2]   ; number of left-going scans + right-going scans
    nsweeps = nscans/2 ; number of left-going scans, or right-going scans

; Because of how the arrays were created, find the indicies of R and L
; scans for NBOLOS_FLAGGED_L and WEIGHT_PER_SCAN_L
    idx_r = indgen(nsweeps)*2
    idx_l = indgen(nsweeps)*2 + 1

;-------------------------
; Number of flagged Bolometers
;-------------------------

; Find the number of flagged bolometers per observation, averaged over
; all scans for L and R
    nbf_L = fltarr(nobs) & nbf_R = fltarr(nobs)
    vbf_L = fltarr(nobs) & vbf_R = fltarr(nobs) ; variance
    for iobs=0, nobs-1 do begin
        nbf_R[iobs] = mean(nbolos_flagged_R[iobs,idx_r])
        nbf_L[iobs] = mean(nbolos_flagged_L[iobs,idx_l]) 

        vbf_R[iobs] = variance(nbolos_flagged_R[iobs,idx_r])
        vbf_L[iobs] = variance(nbolos_flagged_L[iobs,idx_l]) 
    endfor

    nbf_mean[ii]    = mean(nbf_r - nbf_l)
    nbf_median[ii]  = median(nbf_r - nbf_l)

    if keyword_set(plotit) then begin
        wset, 0
        plot, nbf_r - nbf_l, title='number of flagged bolos by observation, averaged over scans', xtitle='obs', ytitle='right - left'
    endif

    print, 'Number of flagged Bolos: <nbf_r - nbf_l>   = ', mean(nbf_r - nbf_l), ', median = ', median(nbf_r - nbf_l)
    print, 'Variance of flagged Bolos: <vbf_r - vbf_l> = ', mean(vbf_r - vbf_l), ', median = ', median(vbf_r - vbf_l)
    print, ' '

;-------------------------
; Weights of scans
;-------------------------

; Find the weights per observation, averaged over all scans for L and R
    avgW_L = fltarr(nobs) & avgW_R = fltarr(nobs)
    vavgW_L = fltarr(nobs) & vavgW_R = fltarr(nobs) ; variance
    for iobs=0, nobs-1 do begin
        avgW_R[iobs] = mean(weight_per_scan_R[iobs,idx_r])
        avgW_L[iobs] = mean(weight_per_scan_L[iobs,idx_l]) 

        vavgW_R[iobs] = variance(weight_per_scan_R[iobs,idx_r])
        vavgW_L[iobs] = variance(weight_per_scan_L[iobs,idx_l]) 
    endfor

    avgW_mean[ii]   = mean(avgW_r - avgW_l)
    avgW_median[ii] = median(avgW_r - avgW_l)

    if keyword_set(plotit) then begin
        wset, 1
        plot, avgW_r - avgW_l, title='summed bolowts by observation, averaged over scans', xtitle='obs', ytitle='right - left'
    endif

    print, 'Summed bolowts by scan: <avgW_r - avgW_l>       = ', mean(avgW_r - avgW_l), ', median = ', median(avgW_r - avgW_l)
    print, 'Variance of summed bolowts: <vavgW_r - vavgW_l> = ', mean(vavgW_r - vavgW_l), ', median = ', median(vavgW_r - vavgW_l)

endfor


; look for correlations


if keyword_set(stopit) then stop

END

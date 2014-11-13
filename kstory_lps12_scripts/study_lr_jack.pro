;;;
; NAME: study_lr_jack.pro
; PURPOSE:
;   study lr jack failure
;
; NOTES:
; 1) study_lr_jack, Study L-R leakage into 2010, 2011 maps
; 1.5) run_sss
; 2) plot_hist_single_field, plot histograms and print L-R number of scans asymmetry
; 2.5) plot_hist
; 3) ksrk, study scan flags between kts, rk, cr idfs
;
; MODIFICATION HISTORY:
;  04/02/2012: (KTS) Created from script_0317
;;;


;...................................................................
; investigate L-R leakage
pro study_lr_jack, field_idx, plot_scans=plot_scans, stopit=stopit
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

print, 'Analyzing field ', field_idx, ', name = ', field_name

; get list of idfs
dates = get_lps12_runlist(field_idx, /obs_dates)
idfs  = get_lps12_runlist(field_idx, /idfs)
maps  = get_lps12_runlist(field_idx, /obs_maps)
nidfs = n_elements(idfs)

; ; get full idf file names - OBSOLETE (0402)
; full_idfs = strarr(nidfs)
; spawn, 'ls '+fst.idf_dirs+'field_scan_150*', full_idfs_tmp
; full_dates = extract_date_from_filename(full_idfs_tmp)
; for ii=0, nidfs-1 do begin
;     wh = where(full_dates eq dates[ii], nwh)
;     if (nwh ne 0) then full_idfs[ii] = full_idfs_tmp[wh]
; endfor

; get initial info for arrays
;d = krf()
d = expand_fits_struct(read_spt_fits(idfs[0], include=['SCAN', 'OBSERVATION']))
nscans = d.observation.nscans
nsweep = nscans/2               ; number of left-going scans
nbolos = d.observation.nbolos


;------------------------
; make arrays
;------------------------
nbolos_flagged_R  = intarr(nidfs, nscans)
nbolos_flagged_L  = intarr(nidfs, nscans)
weight_per_scan_R = fltarr(nidfs, nscans)
weight_per_scan_L = fltarr(nidfs, nscans)
dmap_sign         = fltarr(nidfs, nscans)
nscans_R          = intarr(nidfs)
nscans_L          = intarr(nidfs)

;------------------------
; loop over all idfs
;------------------------
for iidf=0, nidfs-1 do begin
;for iidf=147, nidfs-1 do begin

    ; make some noise
    if (iidf mod 50 eq 2) then print, 'read file number '+ strtrim(string(iidf), 2) + '/' + strtrim(string(nidfs), 2)

    ; read in data
    d  = expand_fits_struct(read_spt_fits(idfs[iidf], include=['SCAN', 'ANTENNA0']))
    m  = krf(maps[iidf])

    ; by-scan data
    scan_bolo_flags = d.scan.bolo_flag ; arr[514, 222]

    scan_flag = d.scan.flag ; arr[222]

    ; get bolometer weights
    bolowts = m.processing.bolo_weights ; arr[514]

    ;------------------------
    ; iterate over all scans
    nscans = n_elements(d.scan.start_index)
    ;dmap_sign = intarr(nscans)
    ridx=0 & lidx=0
    print, 'nscans = ', nscans
    for iscan=0, nscans - 1 do begin

        ; Skip scans that have been flagged
        if scan_flag[iscan] eq 0 then begin

            ; determine if the scans are left or right going scans
            ra = d.antenna0.ra
            ra_bs = ra[d.scan[iscan].start_index:d.scan[iscan].stop_index]
            dmap_sign[iidf, iscan] =  (ra_bs[10]-ra_bs[0] > 0)? 1 : -1

            ; check for az-wrap over 0h border
            dmap_sign[iidf, iscan] =  abs(ra_bs[10]-ra_bs[0]) lt 180 ? dmap_sign[iidf, iscan] : -1.0*dmap_sign[iidf, iscan] ;

            ; Right going
            if (dmap_sign[iidf, iscan] gt 0) then begin
                whfail = where(scan_bolo_flags[*,iscan] ne 0, nwhfail, complement=whpass)
                nbolos_flagged_R[iidf, iscan] = nwhfail gt 0 ? nwhfail : 0
                weight_per_scan_R[iidf, iscan] = total( bolowts[whpass] )
                ridx += 1

            ; Left going
            endif else begin
                ; number of flagged bolometers
                whfail = where(scan_bolo_flags[*,iscan] ne 0, nwhfail, complement=whpass)
                nbolos_flagged_L[iidf, iscan] = nwhfail gt 0 ? nwhfail : 0
                weight_per_scan_L[iidf, iscan] = total( bolowts[whpass] )
                lidx += 1
            endelse

        endif ; scan_flag

    endfor ; iscan

    ; fill nscans_LR arrays
    whR = where(weight_per_scan_R[iidf,*] ne 0, nwhR)
    whL = where(weight_per_scan_L[iidf,*] ne 0, nwhL)
    nscans_R[iidf] = nwhR
    nscans_L[iidf] = nwhL
;     print, 'start_time = ', d.scan[0].start_time
    print, 'nwhR, nwhL = ', nwhR, nwhL
    
    ;------------------------
    ; Plot stuff by scan
    if keyword_set(plot_scans) then begin

        ; number of flagged bolometers
        if 0 then begin
            ptitle='number of flagged bolos'
            arrR = nbolos_flagged_R
            arrL = nbolos_flagged_L
        endif

        ; scan weights
        if 1 then begin
            ptitle='weights per scan'
            arrR = weight_per_scan_R
            arrL = weight_per_scan_L
        endif
        
        nbins = 50
        whR = where(arrR[iidf,*] ne 0, nwhR)
        whL = where(arrL[iidf,*] ne 0, nwhL)
        
        mnmx = minmax([ reform(arrL[iidf,whL]), reform(arrR[iidf,whR])])
        xvec = indgen(nbins)*(mnmx[1] - mnmx[0]) / nbins + mnmx[0]
        hhR = histogram(arrR[iidf,whR], nbins=nbins, min=mnmx[0], max=mnmx[1])
        hhL = histogram(arrL[iidf,whL], nbins=nbins, min=mnmx[0], max=mnmx[1])
        
        wset, 0
        plot,  xvec, hhR, psym=10, title='Histogram', xtitle=ptitle, ytitle='nscans'
        oplot, xvec, hhR, psym=10, color = !red
        oplot, xvec, hhL, psym=10, color = !blue
        
        wset, 1
        plot,  arrR[iidf, whR], title='By Scan', xtitle='scan number', ytitle=ptitle
        oplot, arrR[iidf, whR], color=!red
        oplot, arrL[iidf, whL], color=!blue
        
        stop
        
    endif

    
; DEBUGGING
    if iidf eq 147 then begin
        print,  ' iidf = 147.  Stopping' & stop
    endif

; DEBUGGING 2
    if iidf eq 148 then begin
        print, ' iidf = 148.  Stopping' & stop
    endif


endfor    
;stop
; plot a histogram of the number of flagged bolos

arrR = weight_per_scan_R
arrL = weight_per_scan_L

nbins = 50
whgR = where(arrR ne 0, nwhgR)
whgL = where(arrL ne 0, nwhgL)

mnmx = minmax([arrR[whgR], arrL[whgL]])
xvec = indgen(nbins)*(mnmx[1] - mnmx[0]) / nbins + mnmx[0]
hhR = histogram(arrR[whgR], nbins=nbins, min=mnmx[0], max=mnmx[1])
hhL = histogram(arrL[whgL], nbins=nbins, min=mnmx[0], max=mnmx[1])

wset, 0
plot,  xvec, hhR, color = !red, title='Histogram of weights', xtitle='weights', ytitle='number of observations'
oplot, xvec, hhL, color = !blue

; save output
savname='/home/kstory/lps12/scripts/sav_files/lr_test_0402_'+field_name+'.sav'
print, 'sav file: ', savname
save, nbolos_flagged_R, nbolos_flagged_L, weight_per_scan_R, weight_per_scan_L, nscans_R, nscans_L, dmap_sign, filename=savname

if keyword_set(stopit) then stop

END


pro run_sss
for ii=6, 19 do begin
    study_lr_jack, ii
endfor
end

;;;;;;;;;
pro plot_hist_single_field, field_idx, stopit=stopit

field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

savname='/home/kstory/lps12/scripts/sav_files/lr_test_0402_'+field_name+'.sav'
restore, savname

arrR = weight_per_scan_R
arrL = weight_per_scan_L

nbins = 80
whgR = where(arrR ne 0, nwhgR)
whgL = where(arrL ne 0, nwhgL)

mnmx = minmax([arrR[whgR], arrL[whgL]])
;mnmx = [0.0, 2e8]
xvec = indgen(nbins)*(mnmx[1] - mnmx[0]) / nbins + mnmx[0]
hhR = histogram(arrR[whgR], nbins=nbins, min=mnmx[0], max=mnmx[1])
hhL = histogram(arrL[whgL], nbins=nbins, min=mnmx[0], max=mnmx[1])

wset, 0
plot,  xvec, hhR, title = field_name+', Histogram of weights', xtitle='Bolowts of single scan', ytitle='nscans'
oplot, xvec, hhR, color = !red
oplot, xvec, hhL, color = !blue
legend,/top,/right,['Right','Left'],textc=[!red, !blue],chars=lchars,box=0

fname = '/home/kstory/lps12/figs/'+field_name+'_lr_weights_0402'
err = tvread(/png, filename=fname, /nodialog)

wset, 1
plot, nscans_R - nscans_L, title = field_name+', nscan asymmetry', ytitle='nscans_R - nscans_L', xtitle='observation number'

; make 'pass' array for lr jack 
pass = [1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0]
pass_str = strarr(n_elements(pass))
for ii=0, n_elements(pass) -1 do begin
    pass_str[ii] = pass[ii] ? '       PASS' : '       --- FAIL'
endfor

print, strtrim(string(field_idx),2) + ': ' + field_name, pass_str[field_idx]
print, 'total(Weight_left-scans) / total(Weight_right-scans) = ', total(arrL) / total(arrR)
print, 'mean(nscans_R - nscans_L) = ', mean(nscans_R - nscans_L)
print, ' '

if keyword_set(stopit) then stop
end


pro plot_hist
for ii=6, 19 do begin
    plot_hist_single_field, ii
end
end


;...................................................................
; Print Averaged results
; - mean and median of number of flagged bolos L - R
; - mean and median of number of avg bolowts L - R
PRO print_lr_avg, plotit=plotit, stopit=stopit
compile_opt IDL2, HIDDEN

;field_idx = 6
f = lps12_fieldstruct()

list = indgen(14)+6
nfields = n_elements(list)

; make arrays
nbf_mean    = fltarr(nfields) ; number of flagged bolos
nbf_median  = fltarr(nfields)
avgW_mean   = fltarr(nfields) ; average bolowts
avgW_median = fltarr(nfields)


; hard-code chisq for lr jack from 0401
chisq = [10.16,  1.87,  2.06,  9.34, 3.62, 5.59,16.86, 2.74, 8.75, 2.24, 8.16, 9.19, 3.13, 0.84, 3.60, 5.07, 8.80, 4.96, 3.02,14.91]
pte   = [ 0.0708, 0.8670, 0.8405, 0.0963, 0.6057, 0.3481, 0.0048, 0.7406, 0.1193, 0.8145, 0.1477, 0.1017, 0.6801, 0.9742, 0.6087, 0.4074, 0.1175, 0.4208, 0.6962, 0.0107]

for ii=0, nfields-1 do begin

    field_idx = list[ii]
    fst = f[field_idx]
    field_name = fst.name

; print nice statement
    print, ' '
    print, '---' + field_name + '-------------------------------------------'

; get jackknife file
    jackfile = '/home/kstory/lps12/scripts/sav_files/lr_test_0402_'+field_name+'.sav'
    restore, jackfile
    ;stop ; DEBUGGING
    ss = size(nbolos_flagged_L)
    nobs    = ss[1]
    nscans  = ss[2]   ; number of left-going scans + right-going scans
    nsweeps = nscans/2 ; number of left-going scans, or right-going scans

;-------------------------
; Number of flagged Bolometers
;-------------------------

; Find the number of flagged bolometers per observation, averaged over
; all scans for L and R
    nbf_L = fltarr(nobs) & nbf_R = fltarr(nobs)
    vbf_L = fltarr(nobs) & vbf_R = fltarr(nobs) ; variance
    for iobs=0, nobs-1 do begin
        idx_r = where(weight_per_scan_R[iobs,*] ne 0, nwhR)
        idx_l = where(weight_per_scan_L[iobs,*] ne 0, nwhL)

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
        idx_r = where(weight_per_scan_R[iobs,*] ne 0, nwhR)
        idx_l = where(weight_per_scan_L[iobs,*] ne 0, nwhL)

        avgW_R[iobs] = mean(weight_per_scan_R[iobs,idx_r])
        avgW_L[iobs] = mean(weight_per_scan_L[iobs,idx_l]) 

        vavgW_R[iobs] = variance(weight_per_scan_R[iobs,idx_r])
        vavgW_L[iobs] = variance(weight_per_scan_L[iobs,idx_l]) 
    endfor

;     avgW_mean[ii]   = mean(avgW_r - avgW_l)
;     avgW_median[ii] = median(avgW_r - avgW_l)
    avgW_mean[ii]   = mean( (avgW_r - avgW_l) / (avgW_r + avgW_l) )
    avgW_median[ii] = median( (avgW_r - avgW_l) / (avgW_r + avgW_l) )

    if keyword_set(plotit) then begin
        wset, 1
        plot, avgW_r - avgW_l, title='summed bolowts by observation, averaged over scans', xtitle='obs', ytitle='right - left'
    endif

    print, 'Summed bolowts by scan: <avgW_r - avgW_l>       = ', mean(avgW_r - avgW_l), ', median = ', median(avgW_r - avgW_l)
    print, 'Variance of summed bolowts: <vavgW_r - vavgW_l> = ', mean(vavgW_r - vavgW_l), ', median = ', median(vavgW_r - vavgW_l)

    ; DEBUGGING
    wh = where( [avgW_r + avgW_l] eq 0, nwh)
    print, 'nwh = ', nwh

endfor ; loop over fields


;-------------------------
; look for correlations
;-------------------------
hpl = [where(list eq 6), where(list eq 19)]
save_plots = 0

; cut down chisq, pte to size of list
chisq = chisq[list]
pte = pte[list]

; Number of Flagged Bolos
window, 2, xsize=600, ysize=600
xstr = '< <nbf_r>_scan - <nbf_l>_scan >_obs'
!p.multi = [0,1,2]
plot, nbf_mean, chisq, psym=4, title="Difference in Number of flagged Bolos, lr jack", xtitle=xstr, ytitle='chisq'
oplot, nbf_mean[hpl], chisq[hpl], psym=4, color=!red
plot, nbf_mean, pte, psym=4, title="Difference in Number of flagged Bolos, lr jack", xtitle=xstr, ytitle='pte'
oplot, nbf_mean[hpl], pte[hpl], psym=4, color=!red

if save_plots then err = tvread(/png, filename='/home/kstory/lps12/figs/lr_nbf_0403', /nodialog)

; Bolowts
window, 3, xsize=600, ysize=600
xstr = '< (<Wr>-<Wl>) / (<Wr>+<Wl>) >_obs, R - L'
!p.multi = [0,1,2]
plot, avgW_mean, chisq, psym=4, title="Difference in bolowts, lr jack", xtitle=xstr, ytitle='chisq'
oplot, avgW_mean[hpl], chisq[hpl], psym=4, color=!red
plot, avgW_mean, pte, psym=4, title="Difference in bolowts, lr jack", xtitle=xstr, ytitle='pte'
oplot, avgW_mean[hpl], pte[hpl], psym=4, color=!red
if save_plots then err = tvread(/png, filename='/home/kstory/lps12/figs/lr_avgW_0403', /nodialog)

!p.multi = 0
if keyword_set(stopit) then stop

END



;;;;;;;;;;;;;;;;;;
; study uneven scans, lps12 v.s. k11
pro ksrk

ik = read_spt_fits('/home/kstory/lps12/lowellfits/ra5h30dec-55_2008/field_scan_stub_150_20080602_014057.fits')
ir = read_spt_fits('/data/rkeisler/low_ell_fits_03Mar2011/ra5h30dec-55/fits/field_scan_150_20080602_014057.fits')
;io = read_spt_fits('/data/sptdat/run2b/ra5h30dec-55/fits/field_scan_stub_150_20080602_014057.fits')
;ikfull = read_spt_fits('/data/sptdat/run2b/ra5h30dec-55/fits/field_scan_150_20080602_014057.fits')

nscans = n_elements(ir.scan.start_index)
dmap_sign_r = intarr(nscans)
dmap_sign_k = intarr(nscans)

ridx_r = 0 & lidx_r = 0
ridx_k = 0 & lidx_k = 0

for iscan=0, nscans - 1 do begin

    ;;;;;;;;;;;;
    ; r
    ;;;;;;;;;;;;

    scan_flag_r = ir.scan.flag
    ; Skip scans that have been flagged
    if scan_flag_r[iscan] eq 0 then begin
            
;             dmap_sign1[iscan] = median(azrate[d.scan[iscan].start_index:d.scan[iscan].stop_index]) ge azrate_cut ? 1 : -1
        ra_bs = ir.antenna0.ra[ir.scan[iscan].start_index:ir.scan[iscan].stop_index]
        dmap_sign_r[iscan] =  (ra_bs[10]-ra_bs[0] > 0)? 1 : -1
        
        ; Right going
        if (dmap_sign_r[iscan] gt 0) then begin
            ridx_r += 1
        ; Left going
        endif else begin
            lidx_r += 1
        endelse
        
    endif                       ; scan_flag
    
    ;;;;;;;;;;;;
    ; k
    ;;;;;;;;;;;;

    scan_flag_k = ik.scan.flag
    ; Skip scans that have been flagged
    if scan_flag_k[iscan] eq 0 then begin
            
;             dmap_sign1[iscan] = median(azrate[d.scan[iscan].start_index:d.scan[iscan].stop_index]) ge azrate_cut ? 1 : -1
        ra_bs = ik.antenna0.ra[ik.scan[iscan].start_index:ik.scan[iscan].stop_index]
        dmap_sign_k[iscan] =  (ra_bs[10]-ra_bs[0] > 0)? 1 : -1
        
        ; Right going
        if (dmap_sign_k[iscan] gt 0) then begin
            ridx_k += 1
            weight_per_scan_R[iidf, ridx] = total( bolowts[whpass] )
        ; Left going
        endif else begin
            lidx_k += 1
        endelse
        
    endif                       ; scan_flag
    
endfor                          ; iscan

print, 'start_time = ', ik.scan[0].start_time
print, 'r: nwhR, nwhL = ', ridx_r, lidx_r
print, 'k: nwhR, nwhL = ', ridx_k, lidx_k

;;;;;;;;;;;;;;;;;;;


stop

end

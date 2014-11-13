;;;
; NAME: script_0316
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Study L-R leakage into 2010, 2011 maps
;
; MODIFICATION HISTORY:
;  03/16/2012: (KTS) Created
;;;


;...................................................................
; investigate L-R leakage
pro sss, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------

;field_idx = 7
;field_idx = 10
;field_idx = 16

field_arr_ = lps12_fieldstruct()

fst = field_arr_[field_idx]
field_name = fst.name

print, 'Analyzing field ', field_idx, ', name = ', field_name

; get list of idfs
dates = get_lps12_runlist(field_idx, /obs_dates)
idfs  = get_lps12_runlist(field_idx, /idfs)
maps  = get_lps12_runlist(field_idx, /obs_maps)
nidfs = n_elements(idfs)

; get full idf file names
full_idfs = strarr(nidfs)
spawn, 'ls '+fst.idf_dirs+'field_scan_150*', full_idfs_tmp
full_dates = extract_date_from_filename(full_idfs_tmp)
for ii=0, nidfs-1 do begin
    wh = where(full_dates eq dates[ii], nwh)
    if (nwh ne 0) then full_idfs[ii] = full_idfs_tmp[wh]
endfor

; get initial info for arrays
;d = krf()
d = expand_fits_struct(read_spt_fits(idfs[0], include=['SCAN', 'OBSERVATION']))
nscans = d.observation.nscans
nsweep = nscans/2               ; number of left-going scans
nbolos = d.observation.nbolos
Rscans = indgen(nsweep) * 2
Lscans = indgen(nsweep) * 2 + 1


;------------------------
; make arrays
;------------------------
nbolos_flagged_R  = intarr(nidfs, nscans)
nbolos_flagged_L  = intarr(nidfs, nscans)
weight_per_scan_R = fltarr(nidfs, nscans)
weight_per_scan_L = fltarr(nidfs, nscans)
weight_per_scan_L_tmp = fltarr(nidfs, nscans) ; DEBUGGING
nscans_R          = intarr(nidfs)
nscans_L          = intarr(nidfs)

;------------------------
; loop over all idfs
;------------------------
for iidf=0, nidfs-1 do begin
;for iidf=244, 245 do begin

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
    ; determine if these scans are left-going or right-going
    ;az_actual = df.antenna0.track_actual[*,0]
;     az_actual = d.antenna0.ra  ;; substitute RA for az since we are using down-sampled idfs
;     npts = n_elements(az_actual)
;     unwrap_az,az_actual
;     sample_rate = 100. / 6.
;     ttemp = findgen(npts)/sample_rate
;     azrate = deriv(ttemp,az_actual)
;     med_abs_azrate = median(abs(azrate))
;     azrate_cut = 0.2*med_abs_azrate


    ;------------------------
    ; iterate over all scans
    nscans = n_elements(d.scan.start_index)
    dmap_sign = intarr(nscans)
    dmap_sign1 = intarr(nscans)
    ridx=0 & lidx=0
    for iscan=0, nscans - 1 do begin

        ; Skip scans that have been flagged
        if scan_flag[iscan] eq 0 then begin

;             dmap_sign1[iscan] = median(azrate[d.scan[iscan].start_index:d.scan[iscan].stop_index]) ge azrate_cut ? 1 : -1
            ra = d.antenna0.ra
            unwrap_az,ra
            ra_bs = ra[d.scan[iscan].start_index:d.scan[iscan].stop_index]
            dmap_sign[iscan] =  (ra_bs[10]-ra_bs[0] > 0)? 1 : -1

            ; Right going
            if (dmap_sign[iscan] gt 0) then begin
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
    
    stop

    ;------------------------
    ; Plot stuff by scan
    if keyword_set(plot_scans) then begin

        ; number of flagged bolometers
        if 0 then begin
            arrR = nbolos_flagged_R
            arrL = nbolos_flagged_L
        endif

        ; scan weights
        if 1 then begin
            arrR = weight_per_scan_R
            arrL = weight_per_scan_L
        endif
        
        nbins = 50
        whR = where(arrR[iidf,*] ne 0, nwhR)
        whL = where(arrL[iidf,*] ne 0, nwhL)
        
        mnmx = minmax([arrL[iidf,whL], arrR[iidf,whR]])
        xvec = indgen(nbins)*(mnmx[1] - mnmx[0]) / nbins + mnmx[0]
        hhR = histogram(arrR[iidf,whR], nbins=nbins, min=mnmx[0], max=mnmx[1])
        hhL = histogram(arrL[iidf,whL], nbins=nbins, min=mnmx[0], max=mnmx[1])
        
        wset, 0
        plot,  xvec, hhR, color = !red
        oplot, xvec, hhL, color = !blue
        
        wset, 1
        plot,  arrR[iidf, whR], color=!red
        oplot, arrL[iidf, whL], color=!blue
        
        stop
        
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
plot,  xvec, hhR, color = !red, title='Histogram of weights'
oplot, xvec, hhL, color = !blue

savname='/home/kstory/lps12/scripts/sav_files/lr_test_0316_'+field_name+'.sav'
print, 'sav file: ', savname
save, nbolos_flagged_R, nbolos_flagged_L, weight_per_scan_R, weight_per_scan_L, nscans_R, nscans_L, filename=savname
;stop
end


pro run_sss
for ii=6, 19 do begin
    sss, ii
endfor
end

;;;;;;;;;
pro plot_hist_single_field, field_idx, stopit=stopit

field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

savname='/home/kstory/lps12/scripts/sav_files/lr_test_0316_'+field_name+'.sav'
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

fname = '/home/kstory/lps12/figs/'+field_name+'_lr_weights_0316'
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
print, 'total(arrL) / total(arrR) = ', total(arrL) / total(arrR)
print, 'mean(nscans_R - nscans_L) = ', mean(nscans_R - nscans_L)
print, ' '

if keyword_set(stopit) then stop
end


pro plot_hist
for ii=6, 19 do begin
    plot_hist_single_field, ii
end
end

;;;;;;;;;;;;;;;;;;
; study uneven scans
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

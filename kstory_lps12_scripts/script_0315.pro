;;;
; NAME: script_0315
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Study L-R leakage into 2010, 2011 maps
;
; MODIFICATION HISTORY:
;  03/15/2012: (KTS) Created
;;;


;...................................................................
; investigate L-R leakage
pro sss
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------

field_idx = 7
;field_idx = 10
;field_idx = 16

field_arr_ = lps12_fieldstruct()

fst = (*field_arr_[field_idx])
field_name = fst.name

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
nbolos_flagged_R = intarr(nidfs, nscans)
nbolos_flagged_L = intarr(nidfs, nscans)
weight_per_scan_R = intarr(nidfs, nscans)
weight_per_scan_L = intarr(nidfs, nscans)


;------------------------
; loop over all idfs
;------------------------

;for iidf=0, nidfs-1 do begin
for iidf=7, nidfs-1 do begin
    print, 'read file number ', iidf
    ; read in data
    ;d  = krf(idfs[iidf])
    d  = expand_fits_struct(read_spt_fits(idfs[iidf], include=['SCAN']))
    df = expand_fits_struct(read_spt_fits(full_idfs[iidf], include=['ANTENNA0']))
    m  = krf(maps[iidf])

    ; by-scan data
    scan_bolo_flags = d.scan.bolo_flag ; arr[514, 222]

    ; get bolometer weights
    bolowts = m.processing.bolo_weights ; arr[514]

    ;------------------------
    ; determine if these scans are left-going or right-going
    az_actual = df.antenna0.track_actual[*,0]
    npts = n_elements(az_actual)
    unwrap_az,az_actual
    sample_rate = 100.
    ttemp = findgen(npts)/sample_rate
    azrate = deriv(ttemp,az_actual)
    med_abs_azrate = median(abs(azrate))
    azrate_cut = 0.2*med_abs_azrate

    dmap_sign = intarr(nscans)
    ridx=0 & lidx=0

    ;------------------------
    ; iterate over all scans
    for iscan=0, nscans - 1 do begin
        dmap_sign[iscan] = median(azrate[d.scan[iscan].start_index:d.scan[iscan].stop_index]) ge azrate_cut ? 1 : -1

        ; Right going
        if (dmap_sign[iscan] gt 0) then begin
            wh = where(scan_bolo_flags[*,iscan] ne 0, nwh)
            nbolos_flagged_R[iidf, ridx] = nwh
            weight_per_scan_R[iidf, ridx] = total( bolowts[wh] )
            ridx += 1
            print, 'ridx = ', ridx
        ; Left going
        endif else begin
            ; number of flagged bolometers
            wh = where(scan_bolo_flags[*,iscan] ne 0, nwh)
            nbolos_flagged_L[iidf, lidx] = nwh
            weight_per_scan_L[iidf, ridx] = total( bolowts[wh] )
            lidx += 1
            print, 'lidx = ', lidx
        endelse
    endfor

    ;------------------------
    ; Plot stuff by scan

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

    mnmx = minmax([arrL[iidf,*], arrR[iidf,*]])
    nbins = 50
    xvec = indgen(nbins)*(mnmx[1] - mnmx[0]) / nbins + mnmx[0]
    whR = where(arrR[iidf,*] ne 0, nwhR)
    hhR = histogram(arrR[iidf,whR], nbins=nbins, min=mnmx[0], max=mnmx[1])
    whL = where(arrL[iidf,*] ne 0, nwhL)
    hhL = histogram(arrL[iidf,whL], nbins=nbins, min=mnmx[0], max=mnmx[1])

    print, 'nwhR, nwhL = ', nwhR, nwhL
    wset, 0
    plot,  xvec, hhR, color = !red
    oplot, xvec, hhL, color = !blue

    wset, 1
    plot,  arrR[iidf, whR], color=!red
    oplot, arrL[iidf, whL], color=!blue
    stop

    
endfor    
; plot a histogram of the number of flagged bolos
mnmx = minmax([nbolos_flagged_L, nbolos_flagged_R])
nbins = 50
xvec = indgen(nbins)*(mnmx[1] - mnmx[0]) / nbins + mnmx[0]
wh  = where(nbolos_flagged_R ne 0)
hhR = histogram(nbolos_flagged_R[wh], nbins=nbins, min=mnmx[0], max=mnmx[1])
wh  = where(nbolos_flagged_L ne 0)
hhL = histogram(nbolos_flagged_L[wh], nbins=nbins, min=mnmx[0], max=mnmx[1])

wset, 1
plot, xvec, hhR, title='hhR'

wset, 2
plot, xvec, hhL, title='hhL'


;save, nbolos_flagged_R, nbolos_flagged_L
stop
end

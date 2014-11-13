;;;
; NAME: make_runlist_for_lps12
; PURPOSE:
;   Make runlists for lps12
;
; CALLING SEQUENCE: make_lps_runlists(...)
;
; INPUTS: 
;
; OUTPUTS:
;   txt file with runlist dates
;   idl sav file with an array of structs
;
; NOTES: (from TC)
; 1) go through all single-obs maps in a directory, compute weights & rms
;    statistics on them, and make a runlist suitable for handing to
;    coadd_fits_maps.pro with outliers cut based on user-specified cut
;    parameter values.  MEDWT cuts on the median weights in the
;    uniform-coverage region (looks for anomalously high or low
;    weights). RMS cuts on the noise RMS in the center of the map (looks
;    for anomalously high or low values). WTRMS cuts on the product of
;    weights*rms^2, which should be a constant if our weighting scheme is
;    optimal (looks for anomalously high or low values). all cut
;    parameters are in units of the median value of the quantity over all
;    runs. 
;
; 2) BE SURE TO SET THE ISLT KEYWORD FOR LEAD-TRAIL OBS

; Check the example...
; EXAMPLE: make the data release 150 GHz coadd with this call:
;  make_runlist_with_cuts,150,dates150, $
;    mapdir='/data/sptdat/data_release/ra5h30dec-55/maps/ptsrc_unmasked/proj5/', $
;    fac_high_medwt=2., fac_high_wtrms=2., $
;    merge='/data/sptdat/data_release/ra5h30dec-55/runlists/coadd_runlist_ra5h30dec-55_coverage_150.txt'
;  coadd_fits_maps,dates=dates150, $
;    dir='/data/sptdat/data_release/ra5h30dec-55/maps/ptsrc_unmasked/proj5/', $
;    band=150,fileout='data_release_coadd.fits',/nonoise
;
; MODIFICATION HISTORY:
;   17/06/2011: (TC) make_runlist_with_cuts.pro created from data_release_runlist_check.pro, 17Jun11, TC.
;   09/12/2011: (KTS) Creaded from make_runlist_with_cuts.pro
;   09/20/2011: (KTS) Fix nfiles bug for fields with maps in multiple dirs
;   09/20/2011: (KTS) Add path_to_maps in runlist
;   03/10/2012: (KTS) Automatically make the dateOnly runlist also.
;;;



pro make_runlist_for_lps12, band, keepdates, mapdir_arr, outfile=outfile, $
                            fac_high_medwt=fac_high_medwt, $
                            fac_low_medwt=fac_low_medwt, $
                            fac_high_rms=fac_high_rms, $
                            fac_low_rms=fac_low_rms, $
                            fac_high_wtrms=fac_high_wtrms, $
                            fac_low_wtrms=fac_low_wtrms, $
                            fac_high_tweight=fac_high_tweight, $
                            fac_low_tweight=fac_low_tweight, $
                            minpix_stripe=minpix_stripe, $
                            savfile=savfile, $
                            use_saved_data=use_saved_data, $
                            overwrite_savfile=overwrite_savfile, $
                            merge_with_runfile=merge_with_runfile, $
                            islt=islt, $
                            stopit=stopit

;;;;;;;;;;;;;;;;;;;;;
; Deal with input parameters
;;;;;;;;;;;;;;;;;;;;;
if n_elements(fac_high_medwt) eq 0 then fac_high_medwt = 1e6
if n_elements(fac_low_medwt) eq 0 then fac_low_medwt = 1e6
if n_elements(fac_high_rms) eq 0 then fac_high_rms = 1e6
if n_elements(fac_low_rms) eq 0 then fac_low_rms = 1e6
if n_elements(fac_high_wtrms) eq 0 then fac_high_wtrms = 1e6
if n_elements(fac_low_wtrms) eq 0 then fac_low_wtrms = 1e6
if n_elements(fac_high_tweight) eq 0 then fac_high_tweight = 1e6
if n_elements(fac_low_tweight) eq 0 then fac_low_tweight = 1e6
if n_elements(minpix_stripe) eq 0 then minpix_stripe = 0

;;;;;;;;;;;;;;;;;;;;;
; Re-make sav file (only if necessary)
;;;;;;;;;;;;;;;;;;;;;
if keyword_set(use_saved_data) then begin
    restore,savfile
    nfiles = n_elements(dates)
endif else begin
    help, mapdir_arr
    if n_elements(mapdir_arr) eq 0 then begin
        print, "mapdir_arr is undefined! Exiting..." 
        return
    endif
    bandstr = strcompress(string(fix(band[0])),/rem)
    
    ; loop over all maps in map directories
    files = ['']
    for d=0, n_elements(mapdir_arr)-1 do begin
        tmp_files = file_search(mapdir_arr[d] + '/map_*' + bandstr + '_*.fits',count=nfiles)
        files = [files, tmp_files]
    endfor
    files = files[1:*]
    nfiles = n_elements(files)

    if nfiles eq 0 then begin
        print,"MAKE_RUNLIST_FOR_LPS12.PRO:  No map files for "+bandstr+" GHz in directory ", mapdir_arr, ".  Quitting."
        return
    endif
    print,"MAKE_RUNLIST_FOR_LPS12.PRO: "+string(nfiles)+" map files found for "+bandstr+" GHz in directory ", mapdir_arr,"."
    
    ; Initialize arrays
    medwts         = fltarr(nfiles)
    rms_1am_midmap = fltarr(nfiles)
    tweight        = fltarr(nfiles)

    ;;;;;;;;;;;;;;;;;;;;;
    ; step through files, getting some statistics from each
    ;;;;;;;;;;;;;;;;;;;;;

    dates = strarr(nfiles)
    path_to_maps = files

    for i=0,nfiles-1 do begin
    ;for i=10,13 do begin
        mstemp = expand_fits_struct(read_spt_fits(files[i]))
        
        ; Define the center portion of the map
        nx = n_elements(mstemp.map.map[*,0])
        ny = n_elements(mstemp.map.map[0,*])
        xymin = min([nx,ny])
        dxymin = xymin/4
        xmin = nx/2 - dxymin
        xmax = nx/2 + dxymin - 1
        if keyword_set(islt) then begin
            xmin = nx/3 - dxymin/2
            xmax = nx/3 + dxymin/2 - 1
            xmin_t = xmin + nx/3
            xmax_t = xmax + nx/3
        endif
        ymin = ny/2 - dxymin
        ymax = ny/2 + dxymin - 1
        maxwts_stripe = fltarr(ny,nfiles)


        dates[i] = mstemp.mapinfo.start_time
        restemp = mstemp.mapinfo[0].reso_arcmin
        if restemp eq 1. then begin
            map_1am = mstemp.map.map
            binfac = 1
        endif else begin 
            binfac = fix(1./restemp)
            map_big = mstemp.map.map[0L:(nx/binfac)*binfac-1L,0L:(ny/binfac)*binfac-1L]
            map_1am = rebin(map_big,nx/binfac,ny/binfac)
        endelse
        if keyword_set(islt) then begin
            is_trail = 0
            medwts[i] = median(mstemp.weight.map[xmin:xmax,ymin:ymax])
            if medwts[i] eq 0 then begin
                medwts[i] = median(mstemp.weight.map[xmin_t:xmax_t,ymin:ymax])
                is_trail = 1
            endif
            if is_trail then begin
                stripewts = mstemp.weight.map[2L*nx/3-200:2L*nx/3+200,*]
                midmap_1am = map_1am[2L*nx/binfac/3-100:2L*nx/binfac/3+100,ny/binfac/2-100:ny/binfac/2+100]
            endif else begin 
                stripewts = mstemp.weight.map[nx/3-200:nx/3+200,*]
                midmap_1am = map_1am[nx/binfac/3-100:nx/binfac/3+100,ny/binfac/2-100:ny/binfac/2+100]
            endelse
        endif else begin
            medwts[i] = median(mstemp.weight.map[xmin:xmax,ymin:ymax])
            stripewts = mstemp.weight.map[nx/2-200:nx/2+200,*]
            midmap_1am = map_1am[nx/binfac/2-100:nx/binfac/2+100,ny/binfac/2-100:ny/binfac/2+100]
        endelse
        maxwts_stripe[*,i] = max(stripewts,dim=1)
        rmstemp = stddev(midmap_1am)
        atemp = gaussfit_hist(midmap_1am,100,-5.*rmstemp,5.*rmstemp)
        rms_1am_midmap[i] = atemp[2]
        tweight[i] = total(mstemp.weight.map)
        if (i mod 50) eq 1 then print,i+1," maps read."
    endfor
endelse


;;;;;;;;;;;;;;;;;;;;;
; Apply Cuts
;;;;;;;;;;;;;;;;;;;;;

; 1) medwt
; note that the most pernicious thing we want to test
; for is stupidly *low* rms or, equivalently, stupidly *high*
; weights. plain old noisy, low-weight maps won't
mmw_raw = median(medwts)
atemp_mw = gaussfit_hist(medwts,100,mmw_raw/4.,4.*mmw_raw)
mmw = atemp_mw[1]
dmw = atemp_mw[2]
isgood_mw = medwts ge mmw_raw/fac_low_medwt[0] and medwts le mmw_raw*fac_high_medwt[0]

; 2) decl cut
; cut obs with not enough decl. coverage (presumably incomplete)
ngood_stripe = fltarr(nfiles)
for i=0,nfiles-1 do begin
    whtemp = where(maxwts_stripe[*,i] ge mmw/4.,ntemp)
    ngood_stripe[i] = ntemp
endfor
isgood_stripe = ngood_stripe ge minpix_stripe[0]

; 3) rms
; cut obs where the weights are inconsistent with the map rms (weights
; should scale as 1/rms^2)
mrms_raw = median(rms_1am_midmap)
; atemp_rms = gaussfit_hist(rms_1am_midmap,100,0.,2.*mrms_raw)
; mrms = atemp_rms[1]
; drms = atemp_rms[2]
isgood_rms = rms_1am_midmap ge mrms_raw/fac_low_rms[0] and rms_1am_midmap le mrms_raw*fac_high_rms[0]

; 4) wtrms
rw_fac = rms_1am_midmap^2*medwts
mrw_raw = median(rw_fac)
; atemp_rw = gaussfit_hist(rw_fac,100,mrw_raw/4.,4.*mrw_raw)
; mrw = atemp_rw[1]
; drw = atemp_rw[2]
isgood_wtrms = rw_fac ge mrw_raw/fac_low_wtrms[0] and rw_fac le mrw_raw*fac_high_wtrms[0]

; 5) tweight
; cut on the total of the weight map.  This should remove maps that
; have a few bolometers with stupid-high weights
mtw_raw = median(tweight)
isgood_tweight = tweight ge mtw_raw/fac_low_tweight and tweight le mtw_raw*fac_high_tweight


;;;;;;;;;;;;;;;;;;;;;
; write runlist of good files (merge with existing file if requested)
;;;;;;;;;;;;;;;;;;;;;

isinmergefile = bytarr(nfiles) + 1
if n_elements(merge_with_runfile) gt 0 then begin
    readcol,merge_with_runfile[0],mergelist,form='A',/sil
    for i=0,nfiles-1 do begin
        whmerge = where(mergelist eq dates[i],nmerge)
        if nmerge eq 0 then isinmergefile[i] = 0
    endfor
endif
isgood_all = isgood_mw and isgood_stripe and isgood_rms and $
  isgood_wtrms and isgood_tweight and isinmergefile
whgood_all = where(isgood_all,ngood_all)
if ngood_all eq 0 then begin
    print,"MAKE_RUNLIST_WITH_CUTS.PRO: No maps survived cuts!"
    return
endif
keepdates = dates[whgood_all]

if n_elements(outfile) gt 0 then begin
    if outfile[0] ne '' then begin
        my_table = transpose( [ [dates[whgood_all]], [path_to_maps[whgood_all]] ] )
        writetab, my_table, outfile
;         get_lun,lun1
;         openw,lun1,outfile
;         for i=0,nfiles-1 do begin
;             if isgood_mw[i] and isgood_stripe[i] and isgood_rms[i] and $
;               isgood_wtrms[i] and isgood_tweight[i] and $
;               isinmergefile[i] then printf,lun1,dates[i]
;         endfor
;         close,lun1
;         free_lun,lun1
    endif
endif

if keyword_set(stopit) then stop

if n_elements(savfile) gt 0 then begin
    if keyword_set(use_saved_data) then begin
        if keyword_set(overwrite_savfile) then save,keepdates,dates,medwts,maxwts_stripe,rms_1am_midmap,tweight,file=savfile[0]
    endif else begin
        save,keepdates,dates,medwts,maxwts_stripe,rms_1am_midmap,tweight,file=savfile[0]
    endelse
endif

end
    


;...................................................................
; Run on a single field
;
pro make_runlist_for_lps12_single_field, idx
field_arr_ = lps12_fieldstruct()

; The final cutsetup
cutSet = lps12_runlist_cutset()

print, '*** Make runlist for field ' + field_arr_[idx].name
print, '    map dirs: ', field_arr_[idx].autoprocessed_map_dirs

make_runlist_for_lps12, 150, dates150, field_arr_[idx].autoprocessed_map_dirs, islt=field_arr_[idx].lead_trail,$
  fac_high_medwt=cutSet.high_medwt, fac_low_medwt=cutSet.low_medwt, $
  fac_high_rms=cutSet.high_rms, fac_low_rms=cutSet.low_rms, $
  fac_high_wtrms=cutSet.high_wtrms, fac_low_wtrms=cutSet.low_wtrms, $
  fac_high_tweight=cutSet.high_tweight, fac_low_tweight=cutSet.low_tweight, $
  outfile='/data/kstory/projects/lps12/runlists/runlist_lps12_'+field_arr_[idx].name+'.txt',$
  savfile='/data/kstory/projects/lps12/runlists/sav_files/runlist_lps12_'+field_arr_[idx].name+'.sav';, /stopit

; Make dateOnly runlist also
make_runlist_dateOnly_lps12, idx

end

;;;
; NAME: make_sun_jack_list.pro
; PURPOSE:
;   Make runlists for sun up, sun down jack
;
; NOTES:
; 1) outputs: '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_sun_up.txt'
; 2) also makes a sav file with relative field weights, which are the
;    normal weights, scaled by the percent of scans used.
; 3) Inefficiencies: (a) takes fields in order, so does not consider
;    oldest maps in the set that has fewer maps. (b) if the number of
;    sun-contaminated obs is unbalanced, then the jack comparison
;    will not consider all maps.
;
; MODIFICATION HISTORY:
;  03/22/2012: (KTS) Created
;  03/26/2012: (KTS) fix bug in printing to .txt runlist, was
;              previously skipping half the data
;;;


;...................................................................
; Make list of sun jacks for all fields
PRO make_sun_jack_list;, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
outdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
;outdir  = '/data/kstory/projects/lps12/jacks/0322/'
nfields=20

; percent of files used in sun jack, needed for combining field jacks
p_used = fltarr(nfields)

for ii=0, nfields-1 do begin
    field_arr_ = lps12_fieldstruct()
    fst = field_arr_[ii]
    field_name = fst.name

    print, ''
    print, 'WRITE_SUN_JACK_LIST: processing field ' + field_name

    ; get runlist
    maps  = get_lps12_runlist(ii, /xspec_maps)
    dates = get_lps12_runlist(ii, /xspec_dates)
    nmaps = n_elements(maps)

    sunup = intarr(nmaps)      ; 1=> up, 0=> down

    ;------------------------
    ; Get sun position from all xspec maps
    ;------------------------
    for jj=0, nmaps-1 do begin
        
        sunpos, mjd2jd(date_string_to_mjd(dates[jj])), ra, dec
        if dec lt 0 then sunup[jj] = 1
        
    endfor
    
    ;------------------------
    ; Make the sun jack lists
    ;------------------------
    whup = where(sunup eq 1, complement=whdown)
    npair = min([n_elements(whup), n_elements(whdown)])
    whup = whup[0:npair-1] & whdown = whdown[0:npair-1] ; re-size
    uplist   = dates[whup]
    downlist = dates[whdown]

    ; print number of maps dropped
    p_used[ii] = npair*2. / nmaps
    print, 'using ' + strtrim(string(p_used[ii]),2) + '% of the data'

    ;------------------------
    ; Write out
    ;------------------------
    print, 'WRITE_SUN_JACK_LIST: output files: ', outdir+'goodfiles_'+field_name+'_sun_up.txt'
    openw,lun,/get_lun,outdir+'goodfiles_'+field_name+'_sun_up.txt'
    for i=0,npair-1 do $
      printf,lun,uplist[i]
    close,lun
    openw,lun,outdir+'goodfiles_'+field_name+'_sun_down.txt'
    for i=0,npair-1 do $
      printf,lun,downlist[i]
    free_lun,lun
endfor


; make sun_weights
mask_dir = '/data/kstory/projects/lps12/masks/masks_50mJy/'
restore, mask_dir+'field_weights.sav'
w = field_weights
tw = total_weight
area = w*tw

sun_weights = fltarr(nfields)

for ii=0, nfields-1 do begin
    sun_weights[ii] = p_used[ii]*area[ii]
endfor

total_sun_weight = total(sun_weights)
sun_weights /= total_sun_weight

; write out sav file
savname = mask_dir+'sun_weights.sav'
save, sun_weights, total_sun_weight, p_used, filename=savname
stop
END

;;;
; NAME: make_moon_jack_list.pro
; PURPOSE:
;   Make runlists for moon up, moon down jack
;
; NOTES:
; 1) outputs: '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_moon_up.txt'
; 2) also makes a sav file with relative field weights, which are the
;    normal weights, scaled by the percent of scans used.
; 3) Inefficiencies: (a) takes fields in order, so does not consider
;    oldest maps in the set that has fewer maps. (b) if the number of
;    moon-contaminated obs is unbalanced, then the jack comparison
;    will not consider all maps.
;
; MODIFICATION HISTORY:
;  03/22/2012: (KTS) Created
;  03/26/2012: (KTS) fix bug in printing to .txt runlist, was
;              previously skipping half the data
;;;


;...................................................................
; Make list of moon jacks for all fields
PRO make_moon_jack_list;, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
outdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
nfields=20

; percent of files used in moon jack, needed for combining field jacks
p_used = fltarr(nfields)

for ii=0, nfields-1 do begin
    field_arr_ = lps12_fieldstruct()
    fst = field_arr_[ii]
    field_name = fst.name

    print, ''
    print, 'WRITE_MOON_JACK_LIST: processing field ' + field_name

    outname = outdir + 'ground_rms_'+field_name+'.sav'

    ; get runlist
    maps  = get_lps12_runlist(ii, /xspec_maps)
    dates = get_lps12_runlist(ii, /xspec_dates)
    nmaps = n_elements(maps)

    moonup = intarr(nmaps)      ; 1=> up, 0=> down

    ;------------------------
    ; Get moon position from all xspec maps
    ;------------------------
    for jj=0, nmaps-1 do begin
        
        moonpos, mjd2jd(date_string_to_mjd(dates[jj])), ra, dec
        if dec lt 0 then moonup[jj] = 1
        
    endfor
    
    ;------------------------
    ; Make the moon jack lists
    ;------------------------
    whup = where(moonup eq 1, complement=whdown)
    npair = min([n_elements(whup), n_elements(whdown)])
    whup = whup[0:npair-1] & whdown = whdown[0:npair-1] ; re-size
    uplist   = dates[whup]
    downlist = dates[whdown]

    ; print number of maps dropped
    p_used[ii] = npair*2. / nmaps
    print, 'using ' + strtrim(string(p_used[ii]),2) + ' fraction of the data'

    ;------------------------
    ; Write out
    ;------------------------
    print, 'WRITE_MOON_JACK_LIST: output files: ', outdir+'goodfiles_'+field_name+'_moon_up.txt'
    openw,lun,/get_lun,outdir+'goodfiles_'+field_name+'_moon_up.txt'
    for i=0,npair-1 do $
      printf,lun,uplist[i]
    close,lun
    openw,lun,outdir+'goodfiles_'+field_name+'_moon_down.txt'
    for i=0,npair-1 do $
      printf,lun,downlist[i]
    free_lun,lun
endfor


; make moon_weights
mask_dir = '/data/kstory/projects/lps12/masks/masks_50mJy/'
restore, mask_dir+'field_weights.sav'
w = field_weights
tw = total_weight
area = w*tw

moon_weights = fltarr(nfields)

for ii=0, nfields-1 do begin
    moon_weights[ii] = p_used[ii]*area[ii]
endfor

total_moon_weight = total(moon_weights)
moon_weights /= total_moon_weight

; write out sav file
savname = mask_dir+'moon_weights.sav'
save, moon_weights, total_moon_weight, p_used, filename=savname
END

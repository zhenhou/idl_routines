;;;
; NAME: mk_az_quality
; PURPOSE:
;   Make sav files of the az ground quality
;
; NOTES:
; 1) outputs: '/data/kstory/projects/lps12/jacks/az_quality/goodfiles_'+field_name+'_az_lowrms.txt'
;
; MODIFICATION HISTORY:
;  03/19/2012: (KTS) Created, copied from /home/cr/code/spt/runlists/mk_az_quality.pro
;;;

;...................................................................
; Make sav files for az quality
PRO mk_az_quality_sav, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

savdir  = '/data/kstory/projects/lps12/jacks/az_quality/'
savname = savdir + 'ground_rms_'+field_name+'.sav'

;get runlist
dates = get_lps12_runlist(field_idx, /obs_dates)
ndates = n_elements(dates)
idfs = strarr(ndates)
for ii=0, ndates-1 do begin
    idfs[ii]=fst.idf_dirs+'field_scan_stub_150_'+dates[ii]+'.fits'
endfor

; get ground file
restore,'/data/rkeisler/ground2009/ground_rms.sav'
ground_var = ground_rms^2

;------------------------
; Find average ground quality
;------------------------
nn=n_elements(idfs)
ground_quality=fltarr(nn)
for ii=0,nn-1 do begin

    ; make some noise
    if (ii mod 50 eq 2) then print, 'reading file '+strtrim(string(ii),2)+'/'+strtrim(string(nn),2)

    a = read_spt_fits(idfs[ii],include = ['ANTENNA0'])
    
    ground_quality[ii] = mean(ground_var[ fix(a.antenna0.track_actual[0:*:100,0]) ])
endfor

; save the output
print, 'MK_AZ_QUALITY_SAV: file name = ', savname
save, dates, ground_quality, filename=savname
END


PRO mk_az_list, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

outdir  = '/data/kstory/projects/lps12/jacks/az_quality/'
savname = outdir + 'ground_rms_'+field_name+'.sav'

; get the az sav file
restore, savname

; If this field is LT, average the two scans
if (fst.lead_trail and field_idx ne 1) then begin
    ; get lead_list, trail_list
    lt_list_save = '/data/kstory/projects/lps12/runlists/lt/list_'+field_name+'_150.sav'
    restore, lt_list_save
    lead_dates  = extract_date_from_filename(lead_list)
    trail_dates = extract_date_from_filename(trail_list)

    nxspec = n_elements(lead_list)
    gqarr = fltarr(nxspec)
    dates_xspec = strarr(nxspec)

    for ii=0, nxspec-1 do begin
        d1 = lead_dates[ii]
        idx1 = where(dates eq d1)
        d2 = trail_dates[ii]
        idx2 = where(dates eq d2)

        gqarr[ii] = mean( [ground_quality[idx1], ground_quality[idx2] ] )
        dates_xspec[ii] = d1
    endfor

    ground_quality = gqarr
    dates = dates_xspec

    ; save a file with these new quantities
    savdir  = '/data/kstory/projects/lps12/jacks/az_quality/'
    savname = savdir + 'lt_ground_rms_'+field_name+'.sav'
    save, dates, ground_quality, filename=savname
endif

; special case for ra23h30dec-55_2008
if field_idx eq 1 then begin
    ; get FINAL_DATES
    restore, '/data/kstory/projects/lps12/runlists/lt/4Xpair_ra23h30dec-55_2008.sav'

    nxspec = n_elements(final_dates[*,0])
    gqarr = fltarr(nxspec)
    dates_xspec = strarr(nxspec)

    for ii=0, nxspec-1 do begin
        xspec_dates = final_dates[ii,*]
        ;idx=intarr(4)
        tmp = 0
        for jj=0, 3 do begin
            tmp += ground_quality[ where(dates eq xspec_dates[jj]) ]
        endfor

        gqarr[ii] = tmp / 4.
        dates_xspec[ii] = final_dates[ii,0]

    endfor

    ground_quality = gqarr
    dates = dates_xspec

    ; save a file with these new quantities
    savdir  = '/data/kstory/projects/lps12/jacks/az_quality/'
    savname = savdir + 'lt_ground_rms_'+field_name+'.sav'
    save, dates, ground_quality, filename=savname
endif

; Make the runlists
npair = n_elements(ground_quality)
ind = sort(ground_quality)
lowlist = dates[ind[0:npair/2-1]]
highlist = dates[reverse(ind[npair/2:*])]

; Write out
print, 'MK_AZ_LIST: output files: ', outdir+'goodfiles_'+field_name+'_az_lowrms.txt'
openw,lun,/get_lun,outdir+'goodfiles_'+field_name+'_az_lowrms.txt'
for i=0,npair/2-1 do $
  printf,lun,lowlist[i]
close,lun
openw,lun,outdir+'goodfiles_'+field_name+'_az_highrms.txt'
for i=0,npair-1-npair/2 do $
  printf,lun,highlist[i]
free_lun,lun

END


;...................................................................
; Make az_quality lists from scratch
; /resume if the sav file already exisits
PRO mk_az_quality, field_idx, resume=resume
compile_opt IDL2, HIDDEN

field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

; Make some noise
print, 'MK_AZ_QUALITY:  field ', field_name

; make the sav file
if ~keyword_set(resume) then begin
    print, 'MK_AZ_QUALITY:  call mk_az_quality_sav'
    mk_az_quality_sav, field_idx
endif

; make the txt file runlist
print, 'MK_AZ_QUALITY:  call mk_az_list'
mk_az_list, field_idx

END



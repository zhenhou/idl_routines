;;;
; NAME: make_azrms_jack_list.pro
; PURPOSE:
;   Make sav files of the az ground quality
;
; NOTES:
; 1) outputs: '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_lowrms.txt'
;
; MODIFICATION HISTORY:
;  03/19/2012: (KTS) Created, copied from /home/cr/code/spt/runlists/mk_az_quality.pro
;  03/22/2012: (KTS) change name
;  06/14/2012: (KTS) Update to work with new runlist scheme.
;  06/15/2012: (KTS) Add preaz option to use runlist before the azrms_95 cut
;;;

;...................................................................
; Make sav files for az quality
PRO make_azrms_jack_sav, field_idx, preaz=preaz
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

savdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
savname = savdir + 'ground_rms_'+field_name+'.sav'
;if n_elements(preaz) ne 0 then savname = savdir + 'ground_rms_preaz_'+field_name+'.sav'
if keyword_set(preaz) then savname = savdir + 'ground_rms_preaz_'+field_name+'.sav'

;get runlist
dates = get_lps12_runlist(field_idx, /obs_dates)
if keyword_set(preaz) then dates = get_lps12_runlist(field_idx, /obs_dates,/preaz)
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
print, 'MAKE_AZRMS_JACK_SAV: file name = ', savname
save, dates, ground_quality, filename=savname
END


PRO write_azrms_jack_list, field_idx, preaz=preaz
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

outdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
savname = outdir + 'ground_rms_'+field_name+'.sav'
if keyword_set(preaz) then savname = outdir + 'ground_rms_preaz_'+field_name+'.sav'

; get the az sav file
restore, savname

if (fst.lead_trail) then begin

    dates = get_lps12_runlist(field_idx, /obs_dates) ; all dates
    xspec_dates = get_lps12_runlist(field_idx, /xspec_dates) ; coadded dates
    nxspec = n_elements(xspec_dates)
    runlist = fst.runlist

    if keyword_set(preaz) then  begin
        dates = get_lps12_runlist(field_idx, /preaz)
        runlist = '/home/kstory/lps12/runlists/runlist_preaz_'+fst.name+'.txt'

        ; read the first column to get the xspec dates
        readcol,/silent,runlist,a,format='a'
        xspec_dates = a
        nxspec = n_elements(xspec_dates)
    endif

    gqarr = fltarr(nxspec)
    dates_out = strarr(nxspec) ; temp variable
    
    ; get ncol, nrow
    obslist = (read_ascii(runlist)).field1
    nrow=n_elements(obslist[0,*])
    ncol=n_elements(obslist[*,0])

    ; Read the runlist of dates
    case ncol of 
        2: begin
            date_list = strarr(2,nrow)
            readcol,/silent,runlist,a,b,format='a,a'
            date_list[0,*] = a
            date_list[1,*] = b
        end
        4: begin
            date_list = strarr(4,nrow)
            readcol,/silent,runlist,a,b,c,d,format='a,a,a,a'
            date_list[0,*] = a
            date_list[1,*] = b
            date_list[2,*] = c
            date_list[3,*] = d
        end
    endcase

    ; add up the ground quality of the LT pairs
    for ii=0, nxspec-1 do begin
        tmp = 0
        for jj=0, ncol-1 do begin
            tmp += ground_quality[ where(dates eq date_list[jj,ii]) ]
        endfor

        gqarr[ii] = tmp / ncol
        dates_out[ii] = date_list[0,ii]

    endfor

    ground_quality = gqarr
    dates = dates_out

    ; save a file with these new quantities
    savdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
    savname = savdir + 'lt_ground_rms_'+field_name+'.sav'
    if keyword_set(preaz) then savname = savdir + 'lt_ground_rms_preaz_'+field_name+'.sav'
    save, dates, ground_quality, filename=savname

endif

; Make the runlists
npair = n_elements(ground_quality) / 2
ind = sort(ground_quality)
lowlist = dates[ind[0:npair-1]]
highlist = dates[reverse(ind[npair:*])]

; Write out
stub = keyword_set(preaz) ? 'preaz_' : ''
out_file_low = outdir+'goodfiles_'+stub+field_name+'_az_lowrms.txt'
out_file_high = outdir+'goodfiles_'+stub+field_name+'_az_highrms.txt'
print, 'WRITE_AZRMS_JACK_LIST: output files: ', out_file_low

openw,lun,/get_lun,out_file_low
for i=0,npair-1 do $
  printf,lun,lowlist[i]
close,lun
openw,lun,out_file_high
for i=0,npair-1 do $
  printf,lun,highlist[i]
free_lun,lun

END


;...................................................................
; Make az_quality lists from scratch
; /resume if the sav file already exisits
PRO make_azrms_jack_list, field_idx, resume=resume, preaz=preaz
compile_opt IDL2, HIDDEN

field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name
do_preaz = keyword_set(preaz) ? 1 : 0

; Make some noise
print, 'MAKE_AZRMS_JACK_LIST:  field ', field_name

; make the sav file
if ~keyword_set(resume) then begin
    print, 'MAKE_AZRMS_JACK_LIST:  call make_az_jack_sav'
    make_azrms_jack_sav, field_idx, preaz=do_preaz
endif

; make the txt file runlist
print, 'MAKE_AZRMS_JACK_LIST:  call mk_az_list'
write_azrms_jack_list, field_idx, preaz=do_preaz

END




;;;;;;;;;;;;;;;;
; OBSOLETE
;;;;;;;;;;;;;;;;

; ; If this field is LT, average the two scans
; if (fst.lead_trail and field_idx ne 1) then begin
;     ; get lead_list, trail_list
;     lt_list_save = '/data/kstory/projects/lps12/runlists/lt/list_'+field_name+'_150.sav'
;     restore, lt_list_save
;     lead_dates  = extract_date_from_filename(lead_list)
;     trail_dates = extract_date_from_filename(trail_list)

;     nxspec = n_elements(lead_list)
;     gqarr = fltarr(nxspec)
;     dates_xspec = strarr(nxspec)

;     for ii=0, nxspec-1 do begin
;         d1 = lead_dates[ii]
;         idx1 = where(dates eq d1)
;         d2 = trail_dates[ii]
;         idx2 = where(dates eq d2)

;         gqarr[ii] = mean( [ground_quality[idx1], ground_quality[idx2] ] )
;         dates_xspec[ii] = d1
;     endfor

;     ground_quality = gqarr
;     dates = dates_xspec

;     ; save a file with these new quantities
;     savdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
;     savname = savdir + 'lt_ground_rms_'+field_name+'.sav'
;     save, dates, ground_quality, filename=savname
; endif

; ; special case for ra23h30dec-55_2008
; if field_idx eq 1 then begin
;     ; get FINAL_DATES
;     restore, '/data/kstory/projects/lps12/runlists/lt/4Xpair_ra23h30dec-55_2008.sav'

;     nxspec = n_elements(final_dates[*,0])
;     gqarr = fltarr(nxspec)
;     dates_xspec = strarr(nxspec)

;     for ii=0, nxspec-1 do begin
;         xspec_dates = final_dates[ii,*]
;         ;idx=intarr(4)
;         tmp = 0
;         for jj=0, 3 do begin
;             tmp += ground_quality[ where(dates eq xspec_dates[jj]) ]
;         endfor

;         gqarr[ii] = tmp / 4.
;         dates_xspec[ii] = final_dates[ii,0]

;     endfor

;     ground_quality = gqarr
;     dates = dates_xspec

;     ; save a file with these new quantities
;     savdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
;     savname = savdir + 'lt_ground_rms_'+field_name+'.sav'
;     save, dates, ground_quality, filename=savname
; endif


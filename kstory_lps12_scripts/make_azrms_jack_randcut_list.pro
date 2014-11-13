;;;
; NAME: make_azrms_jack_list.pro
; PURPOSE:
;   Make sav files of the az ground quality
;
; NOTES:
; 1) outputs: '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_lowrms.txt'
; 2) Assumes that ground_rms_preaz files have already been made.
;
; MODIFICATION HISTORY:
;  03/19/2012: (KTS) Created, copied from /home/cr/code/spt/runlists/mk_az_quality.pro
;  03/22/2012: (KTS) change name
;  06/14/2012: (KTS) Update to work with new runlist scheme.
;  06/15/2012: (KTS) Add preaz option to use runlist before the azrms_95 cut
;  06/16/2012: (KTS) Change to randcut
;;;

PRO write_azrms_jack_randcut_list, idx, seed, keep_pct
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[idx]
field_name = fst.name

skeep = strtrim(string(keep_pct), 2)
sseed = strtrim(string(seed), 2)

; Directories
outdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
runlist_dir  = '/home/kstory/lps12/runlists/'

;------------------------
; Get the ground_quality array for the preaz runlist
;------------------------
if (fst.lead_trail) then begin
    ; Get the dates and ground_quality for the coadded lt maps
    savname = outdir + 'lt_ground_rms_preaz_'+field_name+'.sav'

; For regular fields
endif else begin
    savname = outdir + 'ground_rms_preaz_'+field_name+'.sav'
endelse

restore, savname


;------------------------
; Apply the randcut runlist
;------------------------
runlist_randcut = runlist_dir+'runlist_randcut_'+skeep+'_seed'+sseed+'_'+field_name+'.txt' ; xspec_dates only
readcol,/silent,runlist_randcut,dates_randcut,format='a'
n_dates_randcut = n_elements(dates_randcut)

; Make new ground_quality array which has only the randcut dates
gq_randcut = fltarr(n_dates_randcut)
for i=0, n_dates_randcut-1 do begin
    gq_randcut[i] = ground_quality[where(dates eq dates_randcut[i])]
endfor

; Make the runlists
npair = n_elements(gq_randcut) / 2
ind = sort(gq_randcut)
lowlist = dates_randcut[ind[0:npair-1]]
highlist = dates_randcut[reverse(ind[npair:*])]

; Write out
stub = 'randcut_'+skeep+'_seed'+sseed
out_file_low = outdir+'goodfiles_'+field_name+'_'+stub+'_az_lowrms.txt'
out_file_high = outdir+'goodfiles_'+field_name+'_'+stub+'_az_highrms.txt'
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
PRO make_azrms_jack_randcut_list, field_idx, seed, keep_pct
compile_opt IDL2, HIDDEN

field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

; Make some noise
print, 'MAKE_AZRMS_JACK_LIST:  field ', field_name

; make the txt file runlist
print, 'MAKE_AZRMS_JACK_LIST:  call mk_az_list'
write_azrms_jack_randcut_list, field_idx, seed, keep_pct

END




;;;;;;;;;;;;;;;;
; OBSOLETE
;;;;;;;;;;;;;;;;
;...................................................................
; Make sav files for az quality
; PRO make_azrms_jack_sav, field_idx, preaz=preaz
; compile_opt IDL2, HIDDEN

; ;------------------------
; ; Setup
; ;------------------------
; field_arr_ = lps12_fieldstruct()
; fst = field_arr_[field_idx]
; field_name = fst.name

; savdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
; savname = savdir + 'ground_rms_'+field_name+'.sav'
; if n_elements(preaz) ne 0 then savname = savdir + 'ground_rms_preaz_'+field_name+'.sav'

; ;get runlist
; dates = get_lps12_runlist(field_idx, /obs_dates)
; if n_elements(preaz) ne 0 then dates = get_lps12_runlist(field_idx, /preaz)
; ndates = n_elements(dates)
; idfs = strarr(ndates)
; for ii=0, ndates-1 do begin
;     idfs[ii]=fst.idf_dirs+'field_scan_stub_150_'+dates[ii]+'.fits'
; endfor

; ; get ground file
; restore,'/data/rkeisler/ground2009/ground_rms.sav'
; ground_var = ground_rms^2

; ;------------------------
; ; Find average ground quality
; ;------------------------
; nn=n_elements(idfs)
; ground_quality=fltarr(nn)
; for ii=0,nn-1 do begin

;     ; make some noise
;     if (ii mod 50 eq 2) then print, 'reading file '+strtrim(string(ii),2)+'/'+strtrim(string(nn),2)

;     a = read_spt_fits(idfs[ii],include = ['ANTENNA0'])
    
;     ground_quality[ii] = mean(ground_var[ fix(a.antenna0.track_actual[0:*:100,0]) ])
; endfor

; ; save the output
; print, 'MAKE_AZRMS_JACK_SAV: file name = ', savname
; save, dates, ground_quality, filename=savname
; END



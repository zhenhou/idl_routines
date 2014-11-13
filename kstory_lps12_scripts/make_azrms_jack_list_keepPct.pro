;;;
; NAME: make_azrms_jack_list_keepPct.pro
; PURPOSE:
;   Make sav files of the az ground quality, cutting a percentage of
;   high azrms observations.
;
; CALLING SEQUENCE: make_azrms_jack_list_keepPct
;
; NOTES:
; 1) outputs: '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_lowrms.txt'
;
; MODIFICATION HISTORY:
;  03/19/2012: (KTS) Created, copied from /home/cr/code/spt/runlists/mk_az_quality.pro
;  03/22/2012: (KTS) change name
;  06/15/2012: (KTS) Update to work with new runlist scheme, but run
;  off of preaz runlists.
;;;

;...................................................................
; Make az_quality lists from scratch
; Assumes the sav files already exist from making the original azrms lists.
; keep_percent should be an integer between 1 and 100.
PRO make_azrms_jack_list_keepPct, field_idx, keep_percent, no_preaz=no_preaz
compile_opt IDL2, HIDDEN
resume = 1

field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

; Make some noise
print, 'MAKE_AZRMS_JACK_LIST:  field ', field_name

; write out the new list
savdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
outdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'

; get the az sav file, which contains dates, ground_quality
if (fst.lead_trail) then begin
    savname = savdir + 'lt_ground_rms_preaz_'+field_name+'.sav'
    if n_elements(no_preaz) ne 0 then savname = savdir + 'lt_ground_rms_'+field_name+'.sav'

endif else begin
    savname = savdir + 'ground_rms_preaz_'+field_name+'.sav'
    if n_elements(no_preaz) ne 0 then savname = savdir + 'ground_rms_'+field_name+'.sav'
endelse

restore, savname


; Make the runlists
npair = n_elements(ground_quality) / 2
ind = sort(ground_quality)
nkeep = floor(keep_percent/100.*npair) ; cut down to the percent specified
lowlist = dates[ind[0:nkeep-1]]
highlist = dates[reverse(ind[nkeep:nkeep*2-1])]

; Write out
kstr = strtrim(string(keep_percent), 2)
out_file_low = outdir+'goodfiles_'+field_name+'_az_lowrms_'+kstr+'.txt'
out_file_high = outdir+'goodfiles_'+field_name+'_az_highrms_'+kstr+'.txt'
print, 'WRITE_AZRMS_JACK_LIST: output files: ', out_file_low
stop

openw,lun,/get_lun,out_file_low
for i=0,nkeep-1 do $
  printf,lun,lowlist[i]
close,lun
openw,lun,out_file_high
for i=0,nkeep-1 do $
  printf,lun,highlist[i]
free_lun,lun

stop
END



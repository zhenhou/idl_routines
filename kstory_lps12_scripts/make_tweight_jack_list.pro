;;;
; NAME: make_tweight_jack_list.pro
; PURPOSE:
;   Make sav files of the az ground quality
;
; NOTES:
; 1) outputs: '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_lowrms.txt'
;
; MODIFICATION HISTORY:
;  03/19/2012: (KTS) Created, copied from /home/cr/code/spt/runlists/mk_az_quality.pro
;  03/22/2012: (KTS) change name
;;;

;...................................................................
; Make sav files for az quality
PRO make_tweight_jack_list, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

print, 'WRITE_TWEIGHT_JACK_LIST: processing field ' + field_name

outdir  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/'
outname = outdir + 'ground_rms_'+field_name+'.sav'

; get runlist
maps  = get_lps12_runlist(field_idx, /xspec_maps)
dates = get_lps12_runlist(field_idx, /xspec_dates)
nmaps = n_elements(maps)

tweight = fltarr(nmaps)

;------------------------
; Get tweights from all xspec maps
;------------------------
for ii=0, nmaps-1 do begin

    if (ii mod 50 eq 2) then print, 'reading map '+strtrim(string(ii),2)+'/'+strtrim(string(nmaps),2)
    d = krf(maps[ii])
    tweight[ii] = total(d.weight.map)

endfor

;------------------------
; Make the tweight jack lists
;------------------------
npair = n_elements(tweight)
ind = sort(tweight)
lowlist = dates[ind[0:npair/2-1]]
highlist = dates[reverse(ind[npair/2:*])]

;------------------------
; Write out
;------------------------
print, 'WRITE_TWEIGHT_JACK_LIST: output files: ', outdir+'goodfiles_'+field_name+'_tweight_low.txt'
openw,lun,/get_lun,outdir+'goodfiles_'+field_name+'_tweight_low.txt'
for i=0,npair/2-1 do $
  printf,lun,lowlist[i]
close,lun
openw,lun,outdir+'goodfiles_'+field_name+'_tweight_high.txt'
for i=0,npair-1-npair/2 do $
  printf,lun,highlist[i]
free_lun,lun

END



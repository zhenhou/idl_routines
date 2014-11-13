;;;
; NAME: make_noise_psd_lps12.pro
; PURPOSE:
;      Make a noise psd from jacknives for lps12
;
; CALLING SEQUENCE:
;      IDL> make_noise_psd_lps12
;
; INPUTS:
;   field_idx
; 
; OUTPUTS:
;    A .sav file with the psd saved in it.
;
; NOTES:
;
; MODIFICATION HISTORY:
;   09/22/2010: (KTS) created from make_noise_psd_0825.pro
;   10/21/2010: (KTS) use all available jacknifes
;   04/26/2012: (KTS) re-named from make_noise_psd_1021.pro, modified for lps12
;;;

pro make_noise_psd_lps12, field_idx
compile_opt IDL2, HIDDEN

; Setup
f = lps12_fieldstruct()
field_name = f[field_idx].name
info = get_lps12_fieldinfo(field_idx)
scratch_dir = '/home/kstory/lps12/scratch/'
out_dir     = '/home/kstory/lps12/noise_psds/'
reso_arcmin = 1.0

; Get list of maps
mapfiles = get_lps12_runlist(field_idx, /xspec_maps)

; Get mask (un-padded)
mask = get_lps12_mask(field_idx, /padded)

; Write out set_def
nn=n_elements(mapfiles)
fruns='tmp_150.txt'
get_lun,u
openw,u,fruns
printf,u,1,nn
printf,u,'Set 1:'
for i=0,nn-1 do printf,u,mapfiles[i]
free_lun,u

psd_150 = estimate_2d_psd_from_maps_multifield(fruns, njacks=100, apodmask=MASK, /overwrite, /expand_fits, $
                                               scrfile=scratch_dir+'noisePSD_scrfile.sav', $
                                               scrwtfile=scratch_dir+'noisePSD_scrwtfile.sav', $
                                               reso_arcmin=reso_arcmin, $
                                               filename=out_dir+'noise_psd_'+field_name+'.sav')

end

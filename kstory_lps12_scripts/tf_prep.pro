;;;
; NAME: tf_prep.pro
; PURPOSE:
;   Make sav file in preparation for making transfer function.
;
; INPUTS:
;   run,        identifier from end2end run, must match what was run in save_mode_coupling_lps12.pro
;
; NOTES:
; 1) expects file at tf_dir/save_mode_coupling.sav, as made by save_mode_coupling_lps12.pro
;
; MODIFICATION HISTORY:
;  04/25/2012: (KTS) Created from /home/rkeisler/ps09/create_twodim_tfs.pro
;  05/23/2012: (KTS) remove commented out mode-coupling code.
;  06/05/2012: (KTS) add PRO init_tf_prep
;;;

PRO tf_prep, idx, run, stopit=stopit, sim_lmax4500=sim_lmax4500
compile_opt IDL2, HIDDEN

field_list = indgen(20)
nfields = n_elements(field_list)

; retrieve default values for structures going into sav file
savfile = 'tf_prep'+run+'.sav'
if (file_test(savfile) eq 1) then restore,savfile


;------------------------
; Setup
;------------------------
f = lps12_fieldstruct()
field_name = f[idx].name
info = get_lps12_fieldinfo(idx) 
nx = info.npix[0] & ny = info.npix[1]
nbig = long(info.nbig)
n_sim_coadds = 100

; directories
sim_dir     = '/home/kstory/lps12/sims/'
mask_dir    = '/home/kstory/lps12/masks/masks_50mJy/'
tf_dir      = '/home/kstory/lps12/twod_tfs/'

bstub = 'coaddsim_lmax8000_'
if keyword_set(sim_lmax4500) then bstub = 'coaddsim_'

; Make some noise
print, 'TF_PREP for field : '+ field_name

; Theory Dl's
;readcol,'/data/rkeisler/low_ell_sims/input/dl_input_true_20101221_061503.txt',lth,dlth
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2

;------------------------
; Make struct 's#' for this field
;------------------------

; get the power
power = dblarr(nbig, nbig)

; get mask_padded
mask_padded = get_lps12_mask(idx, /padded)
mask_factor = mean(mask_padded^2.)
    
; Get coadded sim maps and average
mcfiles = sim_dir+bstub+field_name+'.dat'
maps = read_dat_map(mcfiles, nx, ny, n_sim_coadds)
print, 'TF_PREP: process file, ' + mcfiles
for imap=0, n_sim_coadds-1 do begin
    print, strtrim(string(imap),2)+'/'+strtrim(string(n_sim_coadds-1),2)
    map_padded = pad_array(maps[*,*,imap], nbig)
    power[*,*] += (abs(fft(map_padded*mask_padded))^2.)
endfor        
power *= (1./n_sim_coadds)

ex=execute('s'+strtrim(floor(idx),2)+'={mask_factor:mask_factor, field_name:field_name, nbig:nbig, power:power}')

; Add this information to the sav file
save,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,$
  filename=savfile

if keyword_set(stopit) then stop
END




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize the tf_prep file for a given run
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO init_tf_prep, run ; for example, run='05'
s0=-1 & s1=-1 & s2=-1 & s3=-1 & s4=-1 & s5=-1 & s6=-1 & s7=-1 & s8=-1 & s9=-1
s10=-1 & s11=-1 & s12=-1 & s13=-1 & s14=-1 & s15=-1 & s16=-1 & s17=-1 & s18=-1 & s19=-1

savfile = 'tf_prep'+run+'.sav'

save,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,$ ;l_mode_coupling,mode_couplings, $
  filename=savfile
END


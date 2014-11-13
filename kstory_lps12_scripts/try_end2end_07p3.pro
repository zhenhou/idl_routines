;;;
; NAME: try_end2end.pro
; PURPOSE:
;   Run end_to_end_powspec_lowell for lps12
;   Actually run CR's version
;
; INPUTS:
;
; NOTES:
; 1) This is probably going to be a temporary script that will be renamed
;
; MODIFICATION HISTORY:
;  04/20/2012: (KTS) Created from /home/rkeisler/ps09/try_end2end.pro
;  05/08/2012: (KTS) Add check for argument 'run'
;  05/23/2012: (KTS) Make this for run_04
;  05/24/2012: (KTS) Set correct theoretical Dl
;  05/31/2012: (KTS) Add 'cut_name' option for azrms_90, randcut
;  06/05/2012: (KTS) make this run_05
;  06/05/2012: (KTS) make this run_06: apply calibration factor of 0.825, new (cut azrms_95) runlists
;  06/15/2012: (KTS) make this run_07: Use new beams
;;;

PRO try_end2end_07p3, field_idx, $
                    stopit=stopit, $
                    resume=resume, $ ; generally set this to true
                    use_kweight=use_kweight, $
                    calib=calib, $
                    cut_name=cut_name, $     ; this should be something like cut_name='randcut_90_seed1' or 'azrms_90'
                    savename=savename_in, $
                    save_mc=save_mc, $
                    out_dir=out_dir, $
                    lr=lr       ; set to use DMAP instead of MAP

;-------------------------
; Setup
;-------------------------
f = lps12_fieldstruct()
fst = f[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name
info = get_lps12_fieldinfo(field_idx)
nbig = info.nbig

run = '07p3'

; directories
sim_dir     = '/home/kstory/lps12/sims/' ; soft-linked to /data23/
workdir     = '/home/kstory/lps12/scratch/'+field_name+'/' ; soft-linked to /data23/
if keyword_set(use_kweight) then workdir     = '/home/kstory/lps12/scratch/'+field_name+'_kweight_07p3/'
if (n_elements(cut_name) ne 0) then workdir     = '/home/kstory/lps12/scratch/'+field_name+'_'+cut_name+'/'
mask_dir    = '/home/kstory/lps12/masks/masks_50mJy/'
kern_dir    = mask_dir
if (n_elements(out_dir) eq 0 ) then out_dir     = '/home/kstory/lps12/end2end/'

print, 'workdir = '+workdir
print, 'out_dir = '+out_dir
spawn, 'mkdir '+workdir

flag_calib = 1 ; only for output sav name
if n_elements(calib) eq 0 then begin
    calib = 0.825
    flag_calib = 0 
endif

scalib = sigfig(calib,3)

reso_arcmin = 1.0
info = get_lps12_fieldinfo(field_idx)
npix = info.npix
f0   = info.lpf ; low-pass filter associated with downsampling


;-------------------------
; get MAPFILES
;-------------------------
mapfiles = get_lps12_runlist(field_idx, /xspec_maps)
if keyword_set(lr) then mapname = 'DMAP.MAP' else mapname = 'MAP.MAP'

; runlist cut option
if n_elements(cut_name) ne 0 then begin
    mapfiles = get_lps12_runlist(field_idx, /xspec_maps, cut_name=cut_name)
endif

;-------------------------
; get MCFILES
;-------------------------
mcfiles = [[[sim_dir+'coaddsim_lmax8000_'+field_name+'.dat']]]
mcmapname = 'coadd'

;-------------------------
; get MASK and KERNEL
;   restores 'kernel' and 'mask_padded'
;-------------------------
kernfile = kern_dir+'kern_'+field_name+'.sav'
restore, kernfile


;-------------------------
; get BANDDEF
;-------------------------
delta_l = 50.
min_l = 250.
max_l = 3150.
nl = floor((max_l-min_l)/delta_l)
lhi = dindgen(nl)*delta_l + min_l
banddef = lhi
; this can possibly be much lower:
maxell = max(banddef)*1.5


;-------------------------
; get BEAMFILES
;-------------------------
year = fst.year
beamfiles = get_lps12_beams(year, l_beam, bl_beam)
simbeamfiles = get_lps12_beams(year, l_beam, bl_beam, /sim_lmax8000)
print, 'beamfiles = ', beamfiles
print, 'simbeamfiles = ', simbeamfiles

;-------------------------
; get SPECTRUMFILES
;-------------------------
spectrumfiles = '/home/kstory/lps12/cls_theory/Dls_theory.txt'

;-------------------------
; if desired, get a KWEIGHT
;-------------------------
if keyword_set(use_kweight) then begin
    kmask = get_lps12_kweight(field_idx)
endif

;-------------------------
; Run end_to_end
;-------------------------


end_to_end_powspec, mapfiles, mcfiles, mask_padded, reso_arcmin, $;f0,$
  banddef=banddef, mapname=mapname, $
  kernel=kernel, ellkern=ellkern, $
  mcmapname=mcmapname, $
  maxell=maxell, $
  transferfunction=transfer, $
  beamfiles=beamfiles, $
  simbeamfiles=simbeamfiles,$
  theoryspectrumfiles=spectrumfiles, $
  pixfunc=pixfunc, $
  spectrum=spectrum, covariance=cov, $
  sample_covariance=sample_cov, $
  measurement_covariance=meas_cov, $
  measurement_chisq=meas_chisq, $
  npix=npix, workdir=workdir, $
  kmask=kmask,$
  convol_kern=convol_kern, $;/skip_window_func, $
  $;field=field_name, $
  $;iskips=iskips, $
  $;nbands=nbands, nspectra=nspectra, $
  resume=resume, $
  $
  $
  $
  curlyw=curlyw, $
  raw_spectrum_mc=spectrum_mc_raw, $
  fine_raw_spectrum_mc=spectrum_mc_raw_fine, $
  fine_banddef=banddef_fine, $
  theory_cl=simspec_interp, $
  cov_mc_raw=cov_mc_raw, $
  setdef_mc=setdef_mc, $
  raw_spectrum_data=spectrum_data_raw, $
  cov_data_raw=cov_data_raw, $
  setdef_data=setdef_data, $
  nmodes=nmodes, $
  fine_nmodes=nmodes_fine, $
  transfer_iter=transfer_iter, $
  binned_kern=superkern, $
  sim_binned_kern=sim_superkern, $
  invkern=invkern, $
  sim_invkern=sim_invkern, $
  beam_funcs=beam_interp, $
  $;simbeam_funcs=simbeam_interp, $
  windowfunc=windowfunc, $
  $;sim_windowfunc=sim_windowfunc, $
  all_data_spectra=all_data_spectra, $
  $;all_mc_spectra=all_mc_spectra, $
  est1_data_cov=est1_data_cov, $
  est2_data_cov=est2_data_cov;, $
  ;calib=calib

if n_elements(savename_in) eq 0 then begin
    savename = out_dir+'end_'+field_name+'_'+run
    if keyword_set(use_kweight) then savename = savename+'_kweight'
    if keyword_set(extra_masking) then savename = savename+'_extra_masking'
    if keyword_set(blob) then savename = savename+'_blob'
    if keyword_set(lr) then savename = savename+'_lr'
    if (flag_calib eq 1)  then savename = savename+'_calib'
    if keyword_set(cut_name) then savename = savename+'_'+cut_name
    savename = savename + '.sav'
endif else savename = savename_in
print, 'TRY_END2END_07: save out put in file: ' + savename
save,filename=savename

; remove scratch files
if keyword_set(save_mc) then begin
    dirname = 'mc_'+run
    if keyword_set(use_kweight) then dirname = dirname+'_kweight'
    if keyword_set(extra_masking) then dirname = dirname+'_extra_masking'
    if keyword_set(blob) then dirname = dirname+'_blob'
    spawn,'mkdir '+dirname,sout
    spawn,'cp -r '+workdir+'mc '+dirname+'/'+field_name,sout2
endif 
spawn,'rm -f '+workdir+'mc/*',soutput
spawn,'rm -f '+workdir+'data/*',soutput2


if keyword_set(stopit) then stop
end






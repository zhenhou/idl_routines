;;;
; NAME: try_end2end.pro
; PURPOSE:
;   Run end_to_end_powspec_lowell for lps12
;
; INPUTS:
;
; NOTES:
; 1) This is probably going to be a temporary script that will be renamed
;
; MODIFICATION HISTORY:
;  04/20/2012: (KTS) Created from /home/rkeisler/ps09/try_end2end.pro
;  05/08/2012: (KTS) Add check for argument 'run'
;;;

PRO try_end2end, field_idx, $
                 stopit=stopit, $
                 $;use_kmask=use_kmask, $
                 calib=calib, $
                 use_kweight=use_kweight, $
                 run=run, $
                 savename=savename_in, $
                 save_mc=save_mc, $
                 resume=resume, $ ; generally set this to true
                 lr=lr   ; set to use DMAP instead of MAP

;-------------------------
; Setup
;-------------------------
f = lps12_fieldstruct()
fst = f[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name
info = get_lps12_fieldinfo(field_idx)
nbig = info.nbig

; Check inputs
if n_elements(run) eq 0 then begin
    print, 'TRY_END2END: input argument "run" is mandatory.  usage:'
    print, "  IDL> try_end2end, 5, run='03', /resume, /use_kweight"
    print, 'Returning...'
    return
endif

; directories
sim_dir     = '/home/kstory/lps12/sims/'
workdir     = '/home/kstory/lps12/scratch/'+field_name+'/'
if keyword_set(use_kweight) then workdir     = '/home/kstory/lps12/scratch/'+field_name+'_kweight/'
mask_dir    = '/home/kstory/lps12/masks/masks_50mJy/'
kern_dir    = mask_dir
out_dir     = '/home/kstory/lps12/end2end/'
kweight_dir = '/home/kstory/lps12/twod_kweights/'

spawn, 'mkdir '+workdir

;if n_elements(calib) eq 0 then calib = 0.756
if n_elements(calib) eq 0 then calib = 0.760
scalib = sigfig(calib,3)

reso_arcmin = 1.0
info = get_lps12_fieldinfo(field_idx) & 
npix = info.npix
f0   = info.lpf ; low-pass filter associated with downsampling


;-------------------------
; get MAPFILES
;-------------------------
mapfiles = get_lps12_runlist(field_idx, /xspec_maps)
if keyword_set(lr) then mapname = 'DMAP.MAP' else mapname = 'MAP.MAP'

;-------------------------
; get MCFILES
;-------------------------
;spawn,'ls /data/rkeisler/low_ell_sims/output/coadd_a_'+field_name+'*sav',mcfiles
mcfiles = [[[sim_dir+'coaddsim_'+field_name+'.dat']]]
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
;min_l = 550.
min_l = 250.
;max_l = 3050.
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
;simbeamfiles = get_lps12_beams(year, l_beam, bl_beam, /sim_run_03)
simbeamfiles = get_lps12_beams(year, l_beam, bl_beam, /sim_lmax4500)

;-------------------------
; get EXTRA_TF_STUFF
;-------------------------
;az_avg_lpf, l, field, f0=f0 ; temperature units


;-------------------------
; get SPECTRUMFILES
;-------------------------
;spectrumfiles = '/data/rkeisler/low_ell_sims/input/spectrum_file.txt'
spectrumfiles = '/data/rkeisler/low_ell_sims/input/dl_input_true_20101221_061503.txt'

;-------------------------
; if desired, get a KMASK
;   Obsolete?
;-------------------------
; if keyword_set(use_kmask) then begin
;     kmask = dblarr(nbig, nbig) + 1.
;     ll = make_fft_grid(reso_arcmin/60.*!dtor,nbig,nbig,fx=lx,fy=ly)
;     ll *= (2.*!pi)
;     lx *= (2.*!pi)
;     ly *= (2.*!pi)
;     whmask = where(abs(lx) ge 1100. and $
;                    abs(lx) lt 1500. and $
;                    abs(ly) lt 300.)
;     kmask[whmask]=0.
; endif


;-------------------------
; if desired, get a KWEIGHT
;-------------------------
if keyword_set(use_kweight) then begin
    kmask = get_lps12_kweight(field_idx)
    ;restore,kweight_dir+'weight_2d_'+field_name+'.sav'
    ;if n_elements(kmask) gt 0 then kmask *= weight_2d $
    ;else kmask = weight_2d
endif

end_to_end_powspec_lps12, mapfiles, mcfiles, mask_padded, reso_arcmin, f0,$
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
  field=field_name, $
  iskips=iskips, $
  nbands=nbands, nspectra=nspectra, $
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
  windowfunc=windowfunc, $
  all_data_spectra=all_data_spectra, $
  est1_data_cov=est1_data_cov, $
  est2_data_cov=est2_data_cov, $
  calib=calib

;if ~keyword_set(use_kmask) then kmask=-1 else kmask=1
kmask = 0
if n_elements(savename_in) eq 0 then begin
    savename = out_dir+'end_'+field_name+'_'+run
    if keyword_set(use_kmask) then savename = savename+'_kmask'
    if keyword_set(use_kweight) then savename = savename+'_kweight'
    if keyword_set(extra_masking) then savename = savename+'_extra_masking'
    if keyword_set(blob) then savename = savename+'_blob'
    if keyword_set(lr) then savename = savename+'_lr'
    savename = savename + '.sav'
endif else savename = savename_in
save,filename=savename

; remove scratch files
if keyword_set(save_mc) then begin
    dirname = 'mc_'+run
    if keyword_set(use_kmask) then dirname = dirname+'_kmask'
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






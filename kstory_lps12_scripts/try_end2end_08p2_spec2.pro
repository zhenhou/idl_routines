;;;
; NAME: try_end2end_08p2.pro
; PURPOSE:
;   Run end_to_end_powspec_lowell for lps12, pipeline test2, spec2
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
;  06/28/2012: (KTS) make this run_08: maps 0620, beams rev3.1, calib=1.0
;  07/30/2012: (KTS) make this run_08p2: Pipeline test
;  08/21/2012: (KTS) make this run_08p2_spec2: Pipeline test2, spec2
;;;

PRO try_end2end_08p2_spec2, field_idx, $
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
; Get transfer function, nbins, kernel, etc from normal run_08
;-------------------------
end_file = '/home/kstory/lps12/end2end/end_ra3h30dec-60_08_kweight.sav'
restore, end_file
; Undefine valiables that we want to fill later in call to unbiased_multiset_pspec
undefine, setdef_data
undefine, spectrum_data_raw
undefine, cov_data_raw
undefine, nmodes
undefine, mapname
undefine, all_data_spectra
undefine, convol_kern
undefine, nspectra
undefine, superkern
undefine, invkern
undefine, sim_superkern
undefine, sim_invkern
undefine, invkernmat
undefine, siminvkernmat
undefine, invkernmattr
undefine, siminvkernmattr

;-------------------------
; Setup
;-------------------------
use_kweight=1 ; force kweight
resume=0; force no-resume
field_idx=4 ; ra3h30dec-60 field

f = lps12_fieldstruct()
fst = f[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name
info = get_lps12_fieldinfo(field_idx)
nbig = info.nbig

run = '08p2_spec2'

; directories
sim_dir     = '/home/kstory/lps12/sims/' ; soft-linked to /data23/
workdir     = '/home/kstory/lps12/scratch/'+field_name+'/' ; soft-linked to /data23/
if keyword_set(use_kweight) then workdir     = '/home/kstory/lps12/scratch/'+field_name+'_'+run+'_kweight/'
if (n_elements(cut_name) ne 0) then workdir     = '/home/kstory/lps12/scratch/'+field_name+'_'+cut_name+'/'
mask_dir    = '/home/kstory/lps12/masks/masks_50mJy/'
kern_dir    = mask_dir
if (n_elements(out_dir) eq 0 ) then out_dir     = '/home/kstory/lps12/end2end/'

print, 'workdir = '+workdir
print, 'out_dir = '+out_dir
spawn, 'mkdir '+workdir

flag_calib = 1 ; only for output sav name
if n_elements(calib) eq 0 then begin
    calib = 1.0 ;0.825
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
;mapfiles = get_lps12_runlist(field_idx, /xspec_maps)
;if keyword_set(lr) then mapname = 'DMAP.MAP' else mapname = 'MAP.MAP'
mapfiles = file_search('/data23/kstory/lps12/sims/pipetest2/coaddsim_pipetest2_spec2_'+field_name+'.dat')
mapname = 'coadd'

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
mask_padded = get_lps12_mask(field_idx, /padded)

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
;spectrumfiles = '/home/kstory/lps12/cls_theory/Dls_theory.txt'

;-------------------------
; if desired, get a KWEIGHT
;-------------------------
if keyword_set(use_kweight) then begin
    kmask = get_lps12_kweight(field_idx)
endif

;-------------------------
; Calculate the autospectra
;-------------------------

unbiased_multiset_pspec, mapfiles, mask_padded, reso_arcmin, setdef=setdef_data, $
  spectrum=spectrum_data_raw, cov=cov_data_raw, $;resume=resume, $
  /cmbweighting, banddef=banddef, maxell=maxell, nmodes=nmodes, $
  persistdir=workdir, /auto, mapname=mapname, npix=npix,kmask=kmask, $
  allspectra=all_data_spectra, convol_kern=convol_kern
stop
nsets=(size(setdef_data))[2]
nspectra=(long(nsets)*(nsets+1))/2
help, all_data_spectra
print, 'nspectra = ', nspectra

; get beams
beams=dblarr(n_elements(ellkern), nspectra)
simbeams=dblarr(n_elements(ellkern), nspectra)

for i=0, nsets-1 do begin
    beam=read_ascii(beamfiles[i])
; fold the effect of the lpf used prior to down-sampling, 
; which is a really small effect anyway, like <0.5% in power at all ells
    lpf_ds = az_avg_lpf_lps12(beam.field1[0,*], field_name, f0=f0)
    beam_interp[*, i]=interpol(beam.field1[1, *]*lpf_ds, $
                               beam.field1[0, *], $
                               ellkern)    
    simbeam=read_ascii(simbeamfiles[i])
    simbeam_interp[*, i]=interpol(simbeam.field1[1, *], $
                               simbeam.field1[0, *], $
                               ellkern)    

endfor

k=0
for i=0, nsets-1 do begin
    for j=i, nsets-1 do begin
        idx=where((beam_interp[*, i] ge 0) * $
                  (beam_interp[*, j] ge 0))
        beams[idx, k]=$
          sqrt(beam_interp[idx, i]*beam_interp[idx, j])
        idx=where((beam_interp[*, i] ge 0) * $
                  (beam_interp[*, j] ge 0))
        simbeams[idx, k]=$
          sqrt(simbeam_interp[idx, i]*simbeam_interp[idx, j])
        k++
    endfor
endfor

;; Step 5: Rebin the Mode-mode coupling kernel adding the transfer
;; function, the beam-functions and the the pixelfunction
superkern=dblarr(nbands, nspectra, nbands, nspectra)
invkern=dblarr(nbands, nspectra, nbands, nspectra)
sim_superkern=dblarr(nbands, nspectra, nbands, nspectra)
sim_invkern=dblarr(nbands, nspectra, nbands, nspectra)
defaultskip=1
iskips = intarr(nspectra)
for i=0, nspectra-1 do begin 
    superkern[*, i, *, i]=rebin_coupling_matrix(kernel, ellkern, banddef, $
                                                transfer=transfer[*, i], $
                                                beam=beams[*, i])
    sim_superkern[*, i, *, i]=rebin_coupling_matrix(kernel, ellkern, banddef, $
                                                transfer=transfer[*, i], $
                                                beam=simbeams[*, i])

    wow = lindgen(nbands) + (1l*nbands)*nspectra*(lindgen(nbands))
    iskip = max([0,where(superkern(wow) eq 0)])+1
    iskips[i]=iskip
    ;; leave the first (usually bogus) bin out of the inversion
    invkern[iskip:*, i, iskip:*, i]=invert(/double,reform(superkern[iskip:*, i, iskip:*, i]))
    sim_invkern[iskip:*, i, iskip:*, i]=invert(/double,reform(sim_superkern[iskip:*, i, iskip:*, i]))
endfor

invkernmat=double(reform(invkern, nbands*nspectra, nbands*nspectra))
siminvkernmat=double(reform(sim_invkern, nbands*nspectra, nbands*nspectra))

invkernmattr=double(transpose(invkernmat))
siminvkernmattr=double(transpose(siminvkernmat))

;; Step 6: Apply the binned coupling kernel to the spectra
spectrum=reform(invkernmat##double(spectrum_data_raw), $
                nbands, nspectra)



;-------------------------
; Run end_to_end
;-------------------------

; end_to_end_powspec_lps12, mapfiles, mcfiles, mask_padded, reso_arcmin, f0,$
;   banddef=banddef, mapname=mapname, $
;   kernel=kernel, ellkern=ellkern, $
;   mcmapname=mcmapname, $
;   maxell=maxell, $
;   transferfunction=transfer, $
;   beamfiles=beamfiles, $
;   simbeamfiles=simbeamfiles,$
;   theoryspectrumfiles=spectrumfiles, $
;   pixfunc=pixfunc, $
;   spectrum=spectrum, covariance=cov, $
;   sample_covariance=sample_cov, $
;   measurement_covariance=meas_cov, $
;   measurement_chisq=meas_chisq, $
;   npix=npix, workdir=workdir, $
;   kmask=kmask,$
;   convol_kern=convol_kern, $;/skip_window_func, $
;   field=field_name, $
;   iskips=iskips, $
;   nbands=nbands, nspectra=nspectra, $
;   resume=resume, $
;   $
;   $
;   $
;   curlyw=curlyw, $
;   raw_spectrum_mc=spectrum_mc_raw, $
;   fine_raw_spectrum_mc=spectrum_mc_raw_fine, $
;   fine_banddef=banddef_fine, $
;   theory_cl=simspec_interp, $
;   cov_mc_raw=cov_mc_raw, $
;   setdef_mc=setdef_mc, $
;   raw_spectrum_data=spectrum_data_raw, $
;   cov_data_raw=cov_data_raw, $
;   setdef_data=setdef_data, $
;   nmodes=nmodes, $
;   fine_nmodes=nmodes_fine, $
;   transfer_iter=transfer_iter, $
;   binned_kern=superkern, $
;   sim_binned_kern=sim_superkern, $
;   invkern=invkern, $
;   sim_invkern=sim_invkern, $
;   beam_funcs=beam_interp, $
;   simbeam_funcs=simbeam_interp, $
;   windowfunc=windowfunc, $
;   sim_windowfunc=sim_windowfunc, $
;   all_data_spectra=all_data_spectra, $
;   all_mc_spectra=all_mc_spectra, $
;   est1_data_cov=est1_data_cov, $
;   est2_data_cov=est2_data_cov, $
;   calib=calib

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
print, 'TRY_END2END_08: save out put in file: ' + savename
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
; spawn,'rm -f '+workdir+'mc/*',soutput
; spawn,'rm -f '+workdir+'data/*',soutput2
spawn,'rm -f '+workdir,soutput2


if keyword_set(stopit) then stop
end






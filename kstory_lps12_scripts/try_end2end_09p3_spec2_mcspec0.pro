;;;
; NAME: try_end2end_08p2.pro
; PURPOSE:
;   Run end_to_end_powspec_lowell for lps12, pipeline test2, spec1
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
;  08/21/2012: (KTS) make this run_08p2_spec1: Pipeline test2, spec1
;  08/30/2012: (KTS) make this run_09p2_spec1: Pipeline test2, spec1, run_09, fix beams again
;  09/20/2012: (KTS) make this run_09p3_spec1: Pipeline test3, spec1
;  09/25/2012: (KTS) make this run_09p3_spec2_mcspec0: Pipeline test3, spec2
;;;

PRO try_end2end_09p3_spec2_mcspec0, field_idx, $
                    stopit=stopit, $
                    $;resume=resume, $ ; generally set this to true
                    $;use_kweight=use_kweight, $
                    calib=calib, $
                    cut_name=cut_name, $     ; this should be something like cut_name='randcut_90_seed1' or 'azrms_90'
                    savename=savename_in, $
                    save_mc=save_mc, $
                    out_dir=out_dir, $
                    lr=lr       ; set to use DMAP instead of MAP

;-------------------------
; Setup
;-------------------------
use_kweight=1 ; force kweight
resume=0
field_idx=4 ; ra3h30dec-60 field

f = lps12_fieldstruct()
fst = f[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name
info = get_lps12_fieldinfo(field_idx)
nbig = info.nbig

run = '09p3_spec2_mcspec0'

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
    calib = 1.0
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
mapfiles = file_search('/data23/kstory/lps12/sims/pipetest3/coaddsim_pipetest3_spec2_'+field_name+'.dat')
mapname = 'coadd'

; runlist cut option
if n_elements(cut_name) ne 0 then begin
    mapfiles = get_lps12_runlist(field_idx, /xspec_maps, cut_name=cut_name)
endif

;-------------------------
; get MCFILES
;-------------------------
;mcfiles = [[[sim_dir+'coaddsim_lmax8000_'+field_name+'.dat']]]
;mcmapname = 'coadd'
mcfiles = [[[sim_dir+'pipetest3/coaddsim_pipetest3_spec0_'+field_name+'.dat']]]
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
beamfiles = '/data/sptdat/beams/rev3.1/bl_2009_150.txt'
;simbeamfiles = get_lps12_beams(year, l_beam, bl_beam, /s09)
simbeamfiles = '/data/sptdat/beams/rev3.1/bl_2009_150.txt'
print, 'beamfiles = ', beamfiles
print, 'simbeamfiles = ', simbeamfiles

;-------------------------
; get SPECTRUMFILES
;-------------------------
spectrumfiles = '/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_0.txt'

;-------------------------
; if desired, get a KWEIGHT
;-------------------------
if keyword_set(use_kweight) then begin
    kmask = get_lps12_kweight(field_idx)
endif

;-------------------------
; Run end_to_end
;-------------------------

mask = mask_padded
field=field_name

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




;+ 
;
;Modification history:
;2010/06/29 - pass convol_kernel keyword, ES.
;-


;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 0: Check all the input arguments
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;; required arguments
if n_elements(mapfiles) eq 0 then begin
  message,  'Required argument \"mapfiles\" is undefined'
endif 

if n_elements(mcfiles) eq 0 then begin
  message,  'Required argument \"mcfiles\" is undefined'
endif 

if n_elements(mask) eq 0 then begin
  message,  'Required argument \"mask\" is undefined'
endif 

s=size(mask)

if (s[0] ne 2) then begin
    message, 'Map window must be 2-D'
endif

if (s[1] ne s[2]) then begin 
    message, 'Map window must be SQUARE'
endif

winsize=s[1]

if n_elements(reso_arcmin) eq 0 then begin
  message,  'Required argument \"reso_arcmin\" is undefined'
endif 

if n_elements(spectrumfiles) eq 0 then begin
  message,  'Required argument \"theoryspectrumfiles\" is undefined'
endif 

if n_elements(maxell) eq 0 then begin
  message,  'Required argument \"maxell\" is undefined'
endif 

if n_elements(beamfiles) eq 0 then begin
  message,  'Required argument \"beamfiles\" is undefined'
endif 
if n_elements(simbeamfiles) ne n_elements(beamfiles) then begin
    print,'Using same beam for real and sim data'
    simbeamfiles=beamfiles
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 1: Make the mode-mode coupling kernel.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

; do not make the kernel if it already exists and the resume keyword
; is set

; if not ( keyword_set(resume) and (n_elements(kernel) ne 0) ) then begin
;     kernel=coupling_kernel(mask, reso_arcmin, maxell=maxell, /changevar, $
;                            interp=1000, oversamp=8, /cheby, ellkern=ellkern, $
;                            curlyw=curlyw)
;     reso=double(reso_arcmin)/60*!dtor
;     kernsize=(size(kernel))[1]
;     u=(dindgen(kernsize)+0.5)#replicate(1./(reso*winsize)^4, kernsize)
;     kernel*=u    
; endif

;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 2: Evaluate the spectrum and covariance of the monte-carlo files
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

mcdir=workdir+'/mc'

; do not recompute the mc spectrum if the resume keyword is set AND
; the spectra and covariances have already been set

if not (keyword_set(resume) and (n_elements(spectrum_mc_raw) ne 0) $
        and (n_elements(spectrum_mc_raw_fine) ne 0) $
        and (n_elements(cov_mc_raw) ne 0)) then begin 
    print, "processing monte carlo maps"
    unbiased_multiset_pspec, mcfiles, mask, reso_arcmin, setdef=setdef_mc, $
      spectrum=spectrum_mc_raw, cov=cov_mc_raw, resume=resume, $
      /cmbweighting, banddef=banddef, maxell=maxell, nmodes=nmodes, $
      persistdir=mcdir, /auto, mapname=mcmapname, npix=npix,kmask=kmask, $
      allspectra=all_mc_spectra, convol_kern=convol_kern

; not sure why this is necessary...
;    if 0 then begin
;        spectrum_mc_raw *= (2.*!pi)
;        cov_mc_raw *= ((2.*!pi)^2.)
;        allspectra_mc *= (2.*!pi)
;    endif


    print, "processing monte carlo maps (finely)"
    ;; do this again sampling very fine for the transfer function
    ;; we can use the ffts in the last directory, since we are just rebinning
    unbiased_multiset_pspec, mcfiles, mask, reso_arcmin, setdef=setdef_mc, $
      spectrum=spectrum_mc_raw_fine, resume=1, $
      /cmbweighting, banddef=banddef_fine, maxell=maxell, $
      persistdir=mcdir, /auto, mapname=mcmapname, npix=npix, nmodes=nmodes_fine,$
      kmask=kmask, convol_kern=convol_kern

; not sure why this is necessary...
;    if 0 then begin
;        spectrum_mc_raw_fine *= (2.*!pi)
;    endif
    
endif
    
nbands=n_elements(banddef)
nbands_fine=n_elements(banddef_fine)
n_mc_spectra=n_elements(spectrum_mc_raw)/nbands
spectrum_mc_raw_reformed=reform(spectrum_mc_raw_fine, nbands_fine,$
                                n_mc_spectra)
cov_mc_raw_reformed=reform(cov_mc_raw, nbands, n_mc_spectra, nbands, $
                           n_mc_spectra)


;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 3: Evaluate the spectrum and covariance of the data files
;;
;;;;;;;;;;;;;;;;;;;;;;;;;

datadir=workdir+'/data'

if not (keyword_set(resume) and (n_elements(spectrum_data_raw) ne 0) $
        and (n_elements(cov_data_raw) ne 0)) then begin 
    print, 'processing real maps'
    unbiased_multiset_pspec, mapfiles, mask, reso_arcmin, setdef=setdef_data, $
      spectrum=spectrum_data_raw, cov=cov_data_raw, resume=resume,$
      /cmbweighting, banddef=banddef, maxell=maxell, $
      persistdir=datadir, /auto, mapname=mapname, npix=npix, chisq=meas_chisq,$
      kmask=kmask, allspectra=all_data_spectra, est1_cov=est1_data_cov, $
      est2_cov=est2_data_cov, convol_kern=convol_kern

; recalibrate the data spectrum and covariances.  note that this
; calibration factor depends on the normalization of the beams
; provided in beamfiles.
; CALIB is the number you multiply the SPT maps by.  For 150 GHz
; 2008/09, it should be ~0.76, but we'll have the default be 1.00
; (i.e., no correction) to keep things simple.
    if n_elements(calib) eq 0 then calib = 1.00
    
    spectrum_data_raw *= (calib^2.)
    cov_data_raw *= (calib^4.)
    all_data_spectra *= (calib^2.)
    est1_data_cov *= (calib^4.)
;    est2_data_cov *= (calib^4.) ; not for /auto
    
; not sure why this is necessary...
;    if 0 then begin
;        spectrum_data_raw *= (2.*!pi)
;        cov_data_raw *= ((2.*!pi)^2.)
;        all_data_spectra *= (2.*!pi)
;        est1_data_cov *= ((2.*!pi)^2.)
;        est2_data_cov *= ((2.*!pi)^2.)
;    endif

endif

n_data_spectra=n_elements(spectrum_data_raw)/nbands
spectrum_data_raw_reformed=reform(spectrum_data_raw, nbands, $
                                  n_data_spectra)
cov_data_raw_reformed=reform(cov_data_raw, nbands, $
                             n_data_spectra, nbands, n_data_spectra)

;; Check to see that the set definition for the data is the same as 
;; the set definition for the monte-carlos (ie are there the same
;; number of sets for each)?

if ((size(setdef_data))[2] ne (size(setdef_mc))[2]) then begin
    message, "monte carlo simulation set count does not match data set count"
endif

; if we get here the number of sets for the monte carlos and the
; real data are the same, which is good.  record the number of sets 
; in each.

nsets=(size(setdef_data))[2]
print, 'nsets = ', nsets ; DEBUGGING

; nspectra includes all of the cross spectra
nspectra=(long(nsets)*(nsets+1))/2

;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 4: Evaluate the transfer functions 
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

beams=dblarr(n_elements(ellkern), nspectra)
simbeams=dblarr(n_elements(ellkern), nspectra)
transfer=dblarr(n_elements(ellkern), nspectra)

if n_elements(beamfiles) ne nsets then begin
    message, 'Please specify one beamfile per set',n_elements(beamfiles),nsets
endif


beam_interp=dblarr(n_elements(ellkern), nsets)
simbeam_interp=dblarr(n_elements(ellkern), nsets)

for i=0, nsets-1 do begin
    beam=read_ascii(beamfiles[i])
; fold the effect of the lpf used prior to down-sampling, 
; which is a really small effect anyway, like <0.5% in power at all ells
    lpf_ds = az_avg_lpf_lps12(beam.field1[0,*], field, f0=f0)
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


;
; By default, we make one transfer function per set, and 
; cross spectra are scaled  by the geometric mean of the
; appropriate single set spectra. However it
; is an option to make a transfer function for cross spectrum 
;

beams_for_tf=simbeam_interp
ntfs=nsets ;  this might be changed below

; This is where we need to know the theory spectrum of the monte carlo sims
; The user has the option of inputting one file per set OR
;  a list of files or "components" per set. 
; Thirdly, if the user specifies one component per spectrum (including 
; cross spectra), then a separate transfer function will be made for
; the cross spectra, (rather than using the geometric mean of the
; the single frequency spectra)

s=size(spectrumfiles)

if (s[0] ge 3) then begin
    message, 'List of spectrumfiles should be at most 2 dimensional'
endif


if (nsets eq 1) then begin
    if(s[0] eq 0) then begin
        ncomponents=1
    endif
    if(s[0] eq 1) then begin 
        ncomponents=n_elements(spectrumfiles)
    endif
    if (s[0] eq 2) then begin
        if (s[2] gt 1) then begin
            message, "I'm confused by your list of spectrum files. For a single set give me a 1-D array"
        endif else begin
            ncomponents=s[1] 
        endelse
    endif
        
endif else begin
    if(s[0] eq 0) then begin 
        message, "Please input at least one theory spectrum per set"
    endif
    if(s[0] eq 1) then begin
        ncomponents=1
        if s[1] eq nsets then begin
            calc_cross_transfer_funcs=0
        endif else begin
            if s[1] eq (nsets*(nsets+1))/2 then begin
                calc_cross_transfer_funcs=1
                ntfs=nspectra
                ;; use the geometrically averaged beams for
                ;; computing transfer functions
                beams_for_tf=simbeams
                print, 'NOTICE: Making special transfer functions for cross spectra'
            endif else begin
                message, 'Please use either one spectrum per set, or one theory spectrum for each sets and each cross-spectrum'
            endelse
        endelse 
    endif
    if(s[0] eq 2 and s[2] ne nsets) then begin
        message, 'please input one list of files per set'
    endif
    if(s[0] eq 2 and s[2] eq nsets) then begin
        ncomponents=s[1]
    endif
endelse

bandcenters=([0, banddef]+banddef)/2
bandcenters_fine=([0, banddef_fine]+banddef_fine)/2

;; there is a chance that the monte-carlo spectrum as provided is 
;; sampled at different ells than the kernel is expecting 
simspec_interp=dblarr(n_elements(ellkern), ntfs)

niter=5
transfer_iter=dblarr(n_elements(ellkern), ntfs, niter+1)

fskyw2=mean(mask^2)

spectrumfiles=reform(spectrumfiles, ncomponents, ntfs)

;
; ok tricky detail, (that Martin imagines people will hate him for
; later):
; If you have n sets, unbiased_multiset_pspec takes n*(n-1)
; cross spectra, and intersperses them with the in-set auto-spectra:
; autospectrum set1xset1 is followed by set1xset2, set1xset3... set1xsetN, 
; which is then followed by set2xset2 and then set2xset3 etc.
;

spectrum_idx=0

for i=0, ntfs-1 do begin 
    ; the fine raw mc spectra should already be at the same binning as the 
    ; the kernel, but just in case we interpolate them to the same ell
    clmc=interpol(spectrum_mc_raw_reformed[*, spectrum_idx],$
                  bandcenters_fine, ellkern)*1d12; convert to uK^2
    if(keyword_set(calc_cross_transfer_funcs)) then begin
        spectrum_idx++          
    endif else begin
      ; increment spectrum_idx to account for
      ; nsets-1B-i cross spectra between auto-spectra
        spectrum_idx+=(nsets-i)
    endelse
    simspec_interp[*, i]=0
    simspec=0.
    for j=0, ncomponents-1 do begin
        if strlen(spectrumfiles[j, i]) eq 0 then continue
        if not file_test(spectrumfiles[j, i]) then begin
            message, 'theory spectrum file: "'+spectrumfiles[j, i]+'" not found'
        endif
        spectrumdata=read_ascii(spectrumfiles[j, i])
        simspec_interp[*, i]+=interpol(spectrumdata.field1[1, *], $
                                       spectrumdata.field1[0, *], $
                                       ellkern)
    endfor


;    simspec_interp *= (2.*!pi)
    ellfact=ellkern*(ellkern+1);/2./!pi;2PI hacked by RK
    transfer_iter[*, i, 0]=transfer_initial(ellkern, simspec_interp[*, i]/ellfact, $
                                            beams_for_tf[*, i], $
                                            fskyw2, clmc/ellfact)
    for j=1, niter do begin
        print, i
        transfer_iter[*, i, j]=transfer_iterate(ellkern,$
                                                simspec_interp[*, i]/ellfact,$
                                                beams_for_tf[*, i],$
                                                fskyw2, $
                                                clmc/ellfact, $
                                                reform(transfer_iter[*,i, j-1]),$
                                                kernel)
    endfor

endfor



k=0
for i=0, nsets-1 do begin
    for j=i, nsets-1 do begin
        idx=where((transfer_iter[*, i, niter] ge 0) * $
                  (transfer_iter[*, j, niter] ge 0))              
        if( keyword_set(calc_cross_transfer_funcs)) then begin
            
            transfer[idx, k]=transfer_iter[idx, k, niter]
            
        endif else begin 
            transfer[idx, k]=$
              sqrt(transfer_iter[idx, i, niter]$
                   *transfer_iter[idx, j, niter])
        endelse
        k++
    endfor
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 5: Rebin the Mode-mode coupling kernel adding the transfer
;; function, the beam-functions and the the pixelfunction
;;
;;;;;;;;;;;;;;;;;;;;;;;;;

; the large multiset binned kernel

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

;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 6: Apply the binned coupling kernel to the spectra
;;
;;;;;;;;;;;;;;;;;;;;;;;;;

spectrum=reform(invkernmat##double(spectrum_data_raw), $
                nbands, nspectra)

;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 7: Apply the binned coupling kernel to the covariances
;; 
;;;;;;;;;;;;;;;;;;;;;;;;;

sample_cov=reform(siminvkernmat##(double(cov_mc_raw)##siminvkernmattr),$
                  nbands, nspectra, nbands, nspectra)
meas_cov=reform(invkernmat##(double(cov_data_raw)##invkernmattr),$
                nbands, nspectra, nbands, nspectra)

;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 8: Sum the covariance to get the total covariance
;;
;;;;;;;;;;;;;;;;;;;;;;;;;

cov=meas_cov+sample_cov

if keyword_set(stopit) then stop

;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 9: Calculate window functions
;;
;;;;;;;;;;;;;;;;;;;;;;;;;
if not keyword_set(skip_window_func) then begin
    windowfunc = dblarr(4000-50-1,nbands*nspectra)
    sim_windowfunc = dblarr(4000-50-1,nbands*nspectra)
    for i=0, nspectra-1 do begin
        ; data
        iskip=iskips[i]
        transstruct = {ell:ellkern,$
                       kernel:kernel,$
                       transfer:transfer[*, i], $
                       bl:beams[*, i]}
        
        windowfunc[*,iskip+i*nbands:nbands-1+i*nbands] = $
          spt_window_func_lowell(150,banddef,transstruct=transstruct,nskip=iskip,/newdat)

        ; sims
        sim_transstruct = {ell:ellkern,$
                           kernel:kernel,$
                           transfer:transfer[*, i], $
                           bl:simbeams[*, i]}
        
        sim_windowfunc[*,iskip+i*nbands:nbands-1+i*nbands] = $
          spt_window_func_lowell(150,banddef,transstruct=sim_transstruct,nskip=iskip,/newdat)
    endfor

endif



;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Save output
;;
;;;;;;;;;;;;;;;;;;;;;;;;;


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
print, 'TRY_END2END_09: save out put in file: ' + savename
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
END






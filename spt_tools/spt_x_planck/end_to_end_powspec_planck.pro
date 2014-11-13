pro end_to_end_powspec_planck, mapfiles, mask, reso_arcmin, $
    banddef=banddef, mapname=mapname, $
    kernel=kernel, ellkern=ellkern, $
    maxell=maxell, $
    beamfiles=beamfiles, $
    pixfunc=pixfunc, $
    spectrum=spectrum, $
    npix=npix, workdir=workdir, $
    resume=resume, $ ;; The following optional outputs,  
                        $                ;; are for holding intermediate states
                        $                ;; in the calculation, for cross check purposes
                        $                ;; When starting from scratch all of these should be
                        $                ;; left undefined
    curlyw=curlyw,      $
    raw_spectrum_data=spectrum_data_raw, $
    cov_data_raw=cov_data_raw, $
    invkern=invkern, $
    windowfunc=windowfunc, $
    calib=calib

;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 0: Check all the input arguments
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;; required arguments
if n_elements(mapfiles) eq 0 then begin
    message,  'Required argument \"mapfiles\" is undefined'
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

if n_elements(maxell) eq 0 then begin
    message,  'Required argument \"maxell\" is undefined'
endif

if n_elements(beamfiles) eq 0 then begin
    message,  'Required argument \"beamfiles\" is undefined'
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 1: Make the mode-mode coupling kernel.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

; do not make the kernel if it already exists and the resume keyword
; is set

if not ( keyword_set(resume) and (n_elements(kernel) ne 0) ) then begin
    kernel=coupling_kernel(mask, reso_arcmin, maxell=maxell, /changevar, $
                           interp=1000, oversamp=8, /cheby, ellkern=ellkern, $
                           curlyw=curlyw)
    reso=double(reso_arcmin)/60*!dtor
    kernsize=(size(kernel))[1]
    u=(dindgen(kernsize)+0.5)#replicate(1./(reso*winsize)^4, kernsize)
    kernel*=u    
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 2: Evaluate the spectrum and covariance of the data files
;;
;;;;;;;;;;;;;;;;;;;;;;;;;

print, "Step 2: Evaluate the spectrum and covariance of the data files"
datadir=workdir+'/data'

if not (keyword_set(resume) and (n_elements(spectrum_data_raw) ne 0) $
        and (n_elements(cov_data_raw) ne 0)) then begin
    print, 'processing real maps'

    unbiased_multiset_pspec_test, mapfiles, mask, reso_arcmin, setdef=setdef_data, $
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
    ;est2_data_cov *= (calib^4.)
endif

nbands=n_elements(banddef)
;
;n_data_spectra=n_elements(spectrum_data_raw)/nbands
;spectrum_data_raw_reformed=reform(spectrum_data_raw, nbands, $
;                                  n_data_spectra)
;cov_data_raw_reformed=reform(cov_data_raw, nbands, $
;                             n_data_spectra, nbands, n_data_spectra)


;; Check to see that the set definition for the data is the same as 
;; the set definition for the monte-carlos (ie are there the same
;; number of sets for each)?

;if ((size(setdef_data))[2] ne (size(setdef_mc))[2]) then begin
;    message, "monte carlo simulation set count does not match data set count"
;endif

; if we get here the number of sets for the monte carlos and the
; real data are the same, which is good.  record the number of sets 
; in each.

nsets=(size(setdef_data))[2]
print, 'nsets = ', nsets ; DEBUGGING

; nspectra includes all of the cross spectra
nspectra=(long(nsets)*(nsets+1))/2


;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 3: Evaluate the transfer functions 
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

print, "Step 3: Evaluate the transfer functions"

beams=dblarr(n_elements(ellkern), nspectra)
;simbeams=dblarr(n_elements(ellkern), nspectra)
;transfer=dblarr(n_elements(ellkern), nspectra)

if n_elements(beamfiles) ne nsets then begin
    message, 'Please specify one beamfile per set',n_elements(beamfiles),nsets
endif


beam_interp=dblarr(n_elements(ellkern), nsets)
;simbeam_interp=dblarr(n_elements(ellkern), nsets)

for i=0, nsets-1 do begin
    beam=read_ascii(beamfiles[i])
; fold the effect of the lpf used prior to down-sampling, 
; which is a really small effect anyway, like <0.5% in power at all ells
;    lpf_ds = az_avg_lpf_lps12(beam.field1[0,*], field, f0=f0)
    beam_interp[*, i]=interpol(beam.field1[1, *], $
                               beam.field1[0, *], $
                               ellkern)
    ;beam_interp[*, i]=interpol(beam.field1[1, *]*lpf_ds, $
    ;                           beam.field1[0, *], $
    ;                           ellkern)    
    ;simbeam=read_ascii(simbeamfiles[i])
    ;simbeam_interp[*, i]=interpol(simbeam.field1[1, *], $
    ;                           simbeam.field1[0, *], $
    ;                           ellkern)    

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
        ;simbeams[idx, k]=$
        ;  sqrt(simbeam_interp[idx, i]*simbeam_interp[idx, j])
        k++
    endfor
endfor




;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 4: Rebin the Mode-mode coupling kernel adding the transfer
;; function, the beam-functions and the the pixelfunction
;;
;;;;;;;;;;;;;;;;;;;;;;;;;

; the large multiset binned kernel

print, "Step 4: Rebin the Mode-mode coupling kernel adding the transfer"
print, "function, the beam-functions and the the pixelfunction"

superkern=dblarr(nbands, nspectra, nbands, nspectra)
invkern=dblarr(nbands, nspectra, nbands, nspectra)
;sim_superkern=dblarr(nbands, nspectra, nbands, nspectra)
;sim_invkern=dblarr(nbands, nspectra, nbands, nspectra)
defaultskip=1
iskips = intarr(nspectra)
for i=0, nspectra-1 do begin 
    superkern[*, i, *, i]=rebin_coupling_matrix(kernel, ellkern, banddef, $
;                                                transfer=transfer[*, i], $
                                                beam=beams[*, i])
;    sim_superkern[*, i, *, i]=rebin_coupling_matrix(kernel, ellkern, banddef, $
;                                                transfer=transfer[*, i], $
;                                                beam=simbeams[*, i])

    wow = lindgen(nbands) + (1l*nbands)*nspectra*(lindgen(nbands))
    iskip = max([0,where(superkern(wow) eq 0)])+1
    iskips[i]=iskip
    ;; leave the first (usually bogus) bin out of the inversion
    invkern[iskip:*, i, iskip:*, i]=invert(/double,reform(superkern[iskip:*, i, iskip:*, i]))
;    sim_invkern[iskip:*, i, iskip:*, i]=invert(/double,reform(sim_superkern[iskip:*, i, iskip:*, i]))
endfor

invkernmat=double(reform(invkern, nbands*nspectra, nbands*nspectra))
;siminvkernmat=double(reform(sim_invkern, nbands*nspectra, nbands*nspectra))

invkernmattr=double(transpose(invkernmat))
;siminvkernmattr=double(transpose(siminvkernmat))

;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 5: Apply the binned coupling kernel to the spectra
;;
;;;;;;;;;;;;;;;;;;;;;;;;;

print, "Step 5: Apply the binned coupling kernel to the spectra"
spectrum=reform(invkernmat##double(spectrum_data_raw), $
                nbands, nspectra)


end

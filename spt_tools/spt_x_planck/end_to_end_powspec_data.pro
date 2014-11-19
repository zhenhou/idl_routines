function check_end2end_input, mapfiles, mask, reso_arcmin, $
         maxell=maxell, beamfiles=beamfiles

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

    return, 0

end


pro end_to_end_powspec_data, mapfiles, mask, reso_arcmin, $
    banddef=banddef, mapname=mapname, $
    kernel=kernel, ellkern=ellkern, maxell=maxell, $
    beamfiles=beamfiles, $
    auto=auto, $
    invkernmat=invkernmat, $
    spectrum=spectrum, $
    npix=npix, workdir=workdir, $
    apply_transfer=apply_transfer, $
    read_lps12_transfer=read_lps12_transfer, $
    transfer=transfer, $
    intfile_ident=intfile_ident, $
    delete_intfile=delete_intfile, $
    resume=resume

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;
    ;; Step 0: Check all the input arguments
    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;

    input_status = check_end2end_input(mapfiles, mask, reso_arcmin, $
    maxell=maxell, beamfiles=beamfiles)

    if (input_status ne 0) then begin
        print, "something wrong with the input variables"
        stop
    endif

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;
    ;; Step 1: Make the mode-mode coupling kernel.
    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ; do not make the kernel if it already exists and the resume keyword
    ; is set
    s = size(mask)
    winsize = s[1]

    if not (keyword_set(resume) and (n_elements(kernel) ne 0) ) then begin
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

    datadir=workdir+'/data'
    
    if not (keyword_set(resume) and (n_elements(spectrum_data_raw) ne 0) $
            and (n_elements(cov_data_raw) ne 0)) then begin
        unbiased_multiset_pspec_test, mapfiles, mask, reso_arcmin, setdef=setdef_data, $
        spectrum=spectrum_data_raw, cov=cov_data_raw, resume=resume,$
        /cmbweighting, banddef=banddef, maxell=maxell, $
        persistdir=datadir, auto=auto, mapname=mapname, npix=npix, chisq=meas_chisq,$
        kmask=kmask, allspectra=all_data_spectra, /no_cross_set, $
        intfile_ident=intfile_ident, delete_intfile=delete_intfile;, $
        ;est1_cov=est1_data_cov, $
        ;est2_cov=est2_data_cov, $
        ;convol_kern=convol_kern
    endif
    
    nbands=n_elements(banddef)
    nsets=(size(setdef_data))[2]
    
    ; nspectra includes all of the cross spectra
    ; nspectra=(long(nsets)*(nsets+1))/2
    nspectra = long(nsets)
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;
    ;; Step 3: the beam transfer functions in data
    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    beams=dblarr(n_elements(ellkern), nspectra)
    
    if n_elements(beamfiles) ne nsets then begin
        message, 'Please specify one beamfile per set',n_elements(beamfiles),nsets
    endif
    
    
    beam_interp=dblarr(n_elements(ellkern), nsets)
    
    for i=0, nsets-1 do begin
        beam=read_ascii(beamfiles[i])
        beam_interp[*, i]=interpol(beam.field1[1, *], $
                                   beam.field1[0, *], $
                                   ellkern)
    endfor
    
    ;; here I assume no_cross_set
    beams = beam_interp
    
    ;k=0
    ;for i=0, nsets-1 do begin
    ;    for j=i, nsets-1 do begin
    ;        idx = where((beam_interp[*, i] ge 0) * $
    ;                   (beam_interp[*, j] ge 0))
    ;        beams[idx, k] = $
    ;        sqrt(beam_interp[idx, i]*beam_interp[idx, j])
    ;        idx = where((beam_interp[*, i] ge 0) * $
    ;                   (beam_interp[*, j] ge 0))
    ;        k++
    ;    endfor
    ;endfor

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;
    ;; Step 3.1: the transfer functions
    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;

    if (keyword_set(apply_transfer)) then begin
        if (keyword_set())
    endif else begin

    endelse


    ;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;
    ;; Step 4: Rebin the Mode-mode coupling kernel adding the transfer
    ;; function, the beam-functions and the the pixelfunction
    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;

    if (n_elements(invkernmat) eq 0) then begin
        print, 'rebin and inverse the kernel'
        superkern=dblarr(nbands, nspectra, nbands, nspectra)
        invkern=dblarr(nbands, nspectra, nbands, nspectra)
        defaultskip=1
        iskips = intarr(nspectra)
        for i=0, nspectra-1 do begin
            superkern[*, i, *, i]=rebin_coupling_matrix(kernel, ellkern, banddef, beam=beams[*, i])

            wow = lindgen(nbands) + (1l*nbands)*nspectra*(lindgen(nbands))
            iskip = max([0,where(superkern(wow) eq 0)])+1
            iskips[i]=iskip
            invkern[iskip:*, i, iskip:*, i]=invert(/double,reform(superkern[iskip:*, i, iskip:*, i]))
        endfor

        invkernmat=double(reform(invkern, nbands*nspectra, nbands*nspectra))
        invkernmattr=double(transpose(invkernmat))
    endif


    ;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;
    ;; Step 5: Apply the binned coupling kernel to the spectra
    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;

    spectrum = reform(invkernmat##double(spectrum_data_raw), nbands, nspectra)
    
    return
end

pro get_cls, path, prefix, dls, out_prefix=out_prefix, camb_path=camb_path, add=add
    ;path = '../samples/'
    ;prefix = 'params_base_camspec_lowl_lowlike_const_omegab_zeq_thetas_thetad'
    
    if (not keyword_set(out_prefix)) then out_prefix=prefix

    filename = path+out_prefix+'.dls.sav'

    info = file_info(filename)
    
    if (info.exists) then begin
        restore, filename, /verbose
    endif else begin
        mcmc = read_chains(1,0,path,prefix,/single)
        lmax = 3000

        dls = dblarr(lmax+1,4,mcmc.nsamples)
    
        for i=1, mcmc.nsamples do begin
            get_sample_info, path, prefix, i-1, params_info, cmb_cls, mcmc=mcmc, camb_path=camb_path, add=add

            dls[0:lmax,0,i-1] = cmb_cls.TT[0:lmax]
            dls[0:lmax,1,i-1] = cmb_cls.EE[0:lmax]
            dls[0:lmax,2,i-1] = cmb_cls.TE[0:lmax]
            dls[0:lmax,3,i-1] = cmb_cls.scalcls_pp[0:lmax]
        endfor

        save, dls, filename=filename, /verbose
    endelse
    
    return
end

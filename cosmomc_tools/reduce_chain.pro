pro reduce_chain, path, prefix, output_path
    
    mcmc = read_chains(16, 1000, path, prefix)

    nparams = mcmc.nparams
    nsamples = mcmc.nsamples
    
    for ip=1L, nparams do begin
        pname = mcmc.paramnames[ip]
        len = strlen(pname)

        if (strmid(pname,0,1,/reverse) eq '*') then pname = 'derived_'+strmid(pname,0,len-1)

        filename = output_path+'/'+pname
        get_lun, unit
        openw, unit, filename
        for i=0L, nsamples-1L do begin
            printf, unit, format='(E16.7)', mcmc.chain_lists[ip,i]
        endfor
        free_lun, unit
    endfor

    filename = output_path+'/loglike'
    get_lun, unit
    openw, unit, filename
    for i=0L, nsamples-1L do begin
        printf, unit, format='(E16.7)', mcmc.chain_lists[0,i]
    endfor
    free_lun, unit

    filename = output_path+'/weight'
    get_lun, unit
    openw, unit, filename
    for i=0L, nsamples-1L do begin
        printf, unit, format='(E16.7)', mcmc.weight[i]
    endfor
    free_lun, unit

    range_file = path+'/'+prefix+'.ranges'
    res = file_info(range_file)
    if res.exists then begin
        spawn, ['cp',range_file,output_path+'/base_params.ranges'], /noshell
    endif
end

pro reunite_chain, path, params, output_path, output_prefix, num_chains=num_chains

    nparams = n_elements(params)
    
    readcol, path+'/loglike', loglike, format='D', nlines=nsamples

    chain_list = dblarr(nparams+1, nsamples)

    chain_list[0,*] = loglike

    param_exists = lonarr(nparams+1)
    param_exists[*] = 1L

    paramnames = params

    for ip=1, nparams do begin
        param_file = path+'/'+params[ip-1]
        res = file_info(param_file)

        if not res.exists then begin
            param_file = path+'/derived_'+params[ip-1]

            paramnames[ip-1] = params[ip-1]+'*'

            res1 = file_info(param_file)
            if not res1.exists then begin
                print, 'Parameter - '+params[ip-1]+', not found'
                param_exists[ip] = 0
            endif
        endif
        
        if (param_exists[ip] eq 1) then begin 
            readcol, param_file, p, format='D'
            chain_list[ip,*] = p
        endif
    endfor

    readcol, path+'/weight', weight, format='D'

    if (not keyword_set(num_chains)) then begin
        num_chains = 8
    endif

    nsamples_chain = long(nsamples/num_chains)

    ip_exists = where(param_exists eq 1)

    for ichain=1, num_chains do begin
        cidx = strcompress(string(ichain),/remove)
        file = output_path+'/'+output_prefix+'_'+cidx+'.txt'
        
        get_lun, unit
        openw, unit, file
        for i=(ichain-1)*nsamples_chain, ichain*nsamples_chain-1 do begin
            printf, unit, format='(E16.7,100E16.7)', weight[i], chain_list[ip_exists,i]
        endfor
        free_lun, unit
    endfor
    

    get_lun, unit
    openw, unit, output_path+'/'+output_prefix+'.paramnames'
    for ip=1, nparams do begin
        if (param_exists[ip] eq 1) then printf, unit, paramnames[ip-1]+'    '+paramnames[ip-1]
    endfor
    free_lun, unit

    range_file = path+'/base_params.ranges'
    res = file_info(range_file)

    ;if res.exists then spawn, ['cp', range_file, output_path+'/'+output_prefix+'.ranges'], /noshell
    
    a = ' '
    if res.exists then begin
        get_lun, unit_in
        get_lun, unit_out
        openr, unit_in, range_file
        openw, unit_out, output_path+'/'+output_prefix+'.ranges'

        while ~eof(unit_in) do begin
            readf, unit_in, a
            res = strsplit(a,/extract)
            
            for ip=1, nparams-1 do begin
                if (res[0] eq paramnames[ip-1]) then printf, unit_out, a
            endfor
        endwhile
        
        free_lun, unit_in
        free_lun, unit_out
    endif
end

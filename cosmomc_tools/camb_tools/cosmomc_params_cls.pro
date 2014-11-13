pro get_params, cmb_params=cmb_params, cos_params=cos_params
    
    cmb_params.Omegabh2  = cos_params.params_set[name2id('omegabh2',paramnames=cos_params.paramnames)]
    cmb_params.Omegach2  = cos_params.params_set[name2id('omegach2',paramnames=cos_params.paramnames)]
    cmb_params.Omeganuh2 = cos_params.params_set[name2id('omeganuh2*',paramnames=cos_params.paramnames)]
    cmb_params.H0        = cos_params.params_set[name2id('H0*',paramnames=cos_params.paramnames)]
    cmb_params.yp        = cos_params.params_set[name2id('yheused*',paramnames=cos_params.paramnames)]
    cmb_params.As        = cos_params.params_set[name2id('A*',paramnames=cos_params.paramnames)]
    cmb_params.ns        = cos_params.params_set[name2id('ns',paramnames=cos_params.paramnames)]
    cmb_params.tau       = cos_params.params_set[name2id('tau',paramnames=cos_params.paramnames)]
    
    cmb_params.Omegak = 0.00d0
    id = name2id('omegak',paramnames=cos_params.paramnames)
    if (id ne 0) then cmb_params.Omegak = cos_params.params_set[id]
    
    cmb_params.w = -1.00d0
    id = name2id('w',paramnames=cos_params.paramnames)
    if (id ne 0) then cmb_params.w = cos_params.params_set[id]

    cmb_params.neff = 3.04600d0
    id = name2id('nnu',paramnames=cos_params.paramnames)
    if (id ne 0) then cmb_params.neff = cos_params.params_set[id]

    cmb_params.Alens = 1.000d0
    id = name2id('Alens',paramnames=cos_params.paramnames)
    if (id ne 0) then cmb_params.Alens = cos_params.params_set[id]

    tags = tag_names(cmb_params)
    ipadd = where(tags eq strupcase('more_pnames'))
    if (ipadd[0] ne -1) then begin
        num_add_params = n_elements(cmb_params.more_pnames)

        for i_add=0, num_add_params-1 do begin
            ip = name2id(cmb_params.more_pnames[i_add],paramnames=cos_params.paramnames)
            if ip[0] ne 0 then begin
                cmb_params.more_params[i_add] = cos_params.params_set[ip]
            endif
        endfor
    endif

    ;print, cmb_params.Omegabh2
    ;print, cmb_params.Omegach2
    ;print, cmb_params.Omeganuh2
    ;print, cmb_params.H0
    ;print, cmb_params.yp
    ;print, cmb_params.As
    ;print, cmb_params.ns
    ;print, cmb_params.tau
    ;print, cmb_params.Omegak
    ;print, cmb_params.w
    ;print, cmb_params.neff

    return
end

function cosmomc_params_cls, lmax=lmax, cos_params=cos_params, old_camb=old_camb, camb_path=camb_path, pivot_k=pivot_k, more_pnames=more_pnames, more_params=more_params, add=add
    
    if (not keyword_set(lmax)) then begin
        lmax = 3300
    endif

    cmb_params = create_struct('lmax',lmax, 'Omegabh2',0.0d0, 'Omegach2',0.0d0, 'Omeganuh2',0.0d0, $
            'Omegak',0.0d0, 'H0',0.0d0, 'w',-1.00d0, 'yp',0.0d0, 'neff',3.046d0, 'As',0.0d0, $
            'ns',0.0d0, 'tau',0.0d0, 'Alens',1.00d0)
    
    if (keyword_set(more_pnames)) then begin
        num_more_params = n_elements(more_pnames)
        if (not keyword_set(more_params)) then more_params = dblarr(num_more_params)

        cmb_params = create_struct(cmb_params, 'more_params',more_params, 'more_pnames',more_pnames)
    endif

    

    get_params, cmb_params=cmb_params, cos_params=cos_params
    cmb_params.lmax = lmax

    rand = randomu(seed,/long)

    output_root = 'idlcamb_'+strcompress(string(rand),/remove)

    exe_camb_spt, cmb_params, output_root, cls, old_camb=old_camb, camb_path=camb_path, pivot_k=pivot_k, add=add
    
    return, cls
end

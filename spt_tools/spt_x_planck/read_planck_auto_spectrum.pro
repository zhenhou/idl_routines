pro read_planck_auto_spectrum, bandcenters, dls_planck_auto, savfile=savfile

    write_savfile = 0
    read_savfile  = 0

    if (keyword_set(savfile)) then begin
        res = file_info(savfile)
        if (res.exists) then read_savfile=1 else write_savfile=1
    endif

    if (not read_savfile) then begin
        workpath     = '/data23/hou/projects/spt_x_planck/powspec/planck'
        ptsrc_workpath = '/data23/hou/projects/spt_x_planck/powspec/planck_ptsrc'

        num_fields = 20
        f = lps12_fieldstruct()

        for ifield=0, num_fields-1 do begin
            field_idx  = ifield
            field_name = f[field_idx].name

            end_savfile = workpath+'/'+field_name+'/data/dls_auto_data_noise.sav'
            res = file_info(end_savfile)

            if (not res.exists) then begin
                print, end_savfile
                print, 'does not exist'
                stop
            endif

            restore, end_savfile

            if ifield eq 0 then begin
                num_bands = n_elements(bandcenters) 
                dls_planck_auto = dblarr(num_bands, 3, num_fields) 
                ;; 0 - data auto spectrum
                ;; 1 - noise spectrum from halfring diff
                ;; 2 - point sources leakage power
            endif

            dls_planck_auto[*,0,ifield] = dls_data*1.0d12
            dls_planck_auto[*,1,ifield] = dls_noise*1.0d12

            ptsrc_savfile = ptsrc_workpath+'/'+field_name+'/data/dls_auto_ptsrc_leakage.sav'
            res = file_info(ptsrc_savfile)
            if (not res.exists) then begin
                print, ptsrc_savfile
                print, 'does not exist'
                stop
            endif

            restore, ptsrc_savfile
            dls_planck_auto[*,2,ifield] = dls_ptsrc
        endfor
    endif else begin
        restore, savfile
    endelse

    if (write_savfile) then begin
        save, bandcenters, dls_planck_auto, filename=savfile
    endif

    return
end

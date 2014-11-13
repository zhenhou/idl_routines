pro read_spt_end2end_spectrum, bandcenters, dls_spt_data, savfile=savfile
    
    write_savfile = 0
    read_savfile  = 0

    if (keyword_set(savfile)) then begin
        res = file_info(savfile)
        if (res.exists) then read_savfile=1 else write_savfile=1
    endif
    
    if (not read_savfile) then begin
        num_fields = 20
        end_path = '/home/kstory/lps12/end2end/'
        fields_ptsrc_corr = ['ra0h50dec-50','ra1hdec-42.5','ra1hdec-60']

        f = lps12_fieldstruct()

        for ifield=0, num_fields-1 do begin
            field_name = f[ifield].name
            
            ip = where(fields_ptsrc_corr eq field_name)
            if ip[0] eq -1 then $
            end_savfile = end_path+'end_'+field_name+'_09_kweight.sav' $
            else $
            end_savfile = end_path+'end_'+field_name+'_ptscor_kweight.sav'
            
            print, 'read -'
            print, end_savfile
            restore, end_savfile

            if (ifield eq 0) then begin
                num_bands = n_elements(banddef) - 1
                dls_spt_data = dblarr(num_bands, num_fields)
            endif

            dls_spt_data[*,ifield] = spectrum[1:*]*1.0d12
        endfor

        nbands=n_elements(banddef)-1
        bandcenters=(banddef[0:nbands-1]+banddef[1:nbands])/2
    endif else begin
        restore, savfile
    endelse

    if (write_savfile) then begin
        save, bandcenters, dls_spt_data, filename=savfile
    endif

    return
end

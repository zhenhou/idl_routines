pro sptsz_field_coord, ifield, theta, phi, coord=coord, map_struct=map_struct, freq=freq, nx=nx, ny=ny, ra0=ra0, dec0=dec0, reso_arcmin=reso_arcmin, proj=proj
    
    fields_info = spt_fieldsinfo()

    if (keyword_set(nx)) then begin
        npixels = [nx, ny]
        radec0 = [ra0, dec0]
    endif else begin
        map_path   = fields_info[ifield].xspec_map_dir
        if (keyword_set(freq)) then begin
            fits_files = file_search(map_path, '*_'+strcompress(string(freq),/remove)+'_*.fits')
        endif else begin
            fits_files = file_search(map_path, '*.fits')
        endelse

        res = read_spt_fits(fits_files[0])
        info = res.mapinfo

        nx = info.nsidex
        ny = info.nsidey
        ra0 = info.ra0
        dec0 = info.dec0
        npixels = [nx, ny]
        radec0  = [ra0, dec0]
        reso_arcmin = info.reso_arcmin
        ;proj = info.projection

        map_struct = expand_fits_struct(res)
    endelse

    pix2ang_anyproj, npixels, radec0, reso_arcmin, ra, dec, proj=proj

    if (keyword_set(coord)) then begin 
        if (coord eq 'G') then begin
            euler, ra, dec, gl, gb, 1
            theta = (90.-gb)*!dtor
            phi = gl*!dtor
        endif else begin
            print, "Now only support coord=G"
            stop
        endelse
    endif else begin
        print, "Please do coord=G for now."
    endelse
    
    return
end

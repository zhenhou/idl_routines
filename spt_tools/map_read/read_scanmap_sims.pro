function read_scanmap_sims, filename

    resolve_routine, 'proj_planck_sptsz'

    hdr = spt_read_binary_header(filename)
    offset = long64(4*8)
    npix = hdr.nx*long64(hdr.ny)
    print,hdr.nx, hdr.ny

    pting  = spt_read_binary_pting(filename)
    weight = spt_read_binary_wt(filename)
    npp = n_elements(pting)
    map_tmp      = fltarr(npix)
    map_original = fltarr(hdr.nx, hdr.ny, hdr.nmap)
    map_weighted = fltarr(hdr.nx, hdr.ny, hdr.nmap)

    map_read = fltarr(npp)
    
    map_tmp[pting] = weight
    array_1dto2d, map_tmp, w2d, nx=hdr.nx, ny=hdr.ny
    
    get_lun, lun
    openr, lun, filename
    point_lun, lun, offset
    for imap=0, hdr.nmap-1 do begin
        readu, lun, map_read
        map_tmp[pting] = map_read

        array_1dto2d, map_tmp, m2d, nx=hdr.nx, ny=hdr.ny
        map_original[*,*,imap] = m2d
        map_weighted[*,*,imap] = m2d * w2d
    endfor
    close, lun
    free_lun, lun

    m = create_struct('maps',map_original, 'maps_weighted',map_weighted)
    
    return, m
end

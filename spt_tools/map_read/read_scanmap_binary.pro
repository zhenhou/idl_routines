function read_scanmap_binary, filename

    hdr = spt_read_binary_header(filename)
    offset = long64(4*8)
    npix = hdr.nx*long64(hdr.ny)

    pting  = spt_read_binary_pting(filename)
    weight = spt_read_binary_wt(filename)
    npp = n_elements(pting)

    map_read = fltarr(npp)
    maps = fltarr(npp,hdr.nmap)

    get_lun, lun
    openr, lun, filename
    point_lun, lun, offset
    for imap=0, hdr.nmap-1 do begin
        readu, lun, map_read
        maps[*,imap] = map_read
    endfor
    close, lun
    free_lun, lun

    m = create_struct('maps',double(maps), 'pting',pting, 'weight',double(weight))

    return, m

end


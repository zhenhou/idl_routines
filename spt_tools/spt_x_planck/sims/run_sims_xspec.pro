pro run_sims_xspec, field_idx, sim_map_root, num_sims, $
    mapname=mapname, $
    workpath=workpath, beamfile=beamfile, $
    bandcenters=bandcenters, dls_sims=dls_sims, $
    dls_sav_root=dls_sav_root, $
    delete_intfile=delete_intfile, sims_type=sims_type, $
    resume=resume

    f = lps12_fieldstruct()
    fst = f[field_idx]
    field_name = fst.name

    reso_arcmin = 1.0
    info = get_lps12_fieldinfo(field_idx)
    npix = info.npix

    if (not keyword_set(sims_type)) then sims_type = 'halfmission'

    ;------------- 
    ; directories 
    ;-------------
    host        = getenv('HOSTNAME')
    if host eq 'spt' then home = '/home/hou/'
    if host eq 'midway' then home = '/home/zhenhou/'

    mask_dir    = home+'data/projects/spt_x_planck/lps12/masks_50mJy/'
    kern_dir    = mask_dir
    workdir     = workpath+'/'+field_name+'/'
    res         = file_info(workdir)
    if (not res.exists) then spawn, ['mkdir', '-p', workdir], /noshell

    ;-------------------------
    ; get MASK and KERNEL
    ;   restores 'kernel' and 'mask_padded'
    ;-------------------------
    maskfile = mask_dir+'mask_'+field_name+'.sav'
    kernfile = kern_dir+'kern_'+field_name+'.sav'
    restore, kernfile
    restore, maskfile

    ;-------------------------
    ; get BANDDEF
    ;-------------------------
    delta_l = 50.
    min_l = 50.
    max_l = 3000.
    nl = floor((max_l-min_l)/delta_l)
    lhi = dindgen(nl)*delta_l + min_l
    banddef = lhi
    ; this can possibly be much lower:
    maxell = max(banddef)*1.2

    nbands=n_elements(banddef)-1
    bandcenters=(banddef[0:nbands-1]+banddef[1:nbands])/2
    
    if (not keyword_set(mapname)) then mapname = 'MAP.MAP'
    if (not keyword_set(dls_sav_root)) then dls_sav_root='dls_'+sims_type+'_sims_xspec'

    sim_map_dir = workdir+'fits_maps/'
    s = file_info(sim_map_dir)
    if (not s.exists) then file_mkdir, sim_map_dir

    sim_dls_dir = workdir+'sim_dls/'
    s = file_info(sim_dls_dir)
    if (not s.exists) then file_mkdir, sim_dls_dir

    for isim=0, num_sims-1 do begin
        sidx = strcompress(string(isim),/remove)
        intfile_ident = 'sim'+sidx

        dls_savfile = sim_dls_dir+dls_sav_root+'_'+intfile_ident+'.sav'
        res = file_info(dls_savfile)

        if res.exists then begin
            restore, dls_savfile
        endif else begin
            ;-------------
            ; sim maps 
            ;-------------

            sim_map_fits1 = sim_map_dir+sim_map_root[0]+'_'+intfile_ident+'.fits'
            sim_map_fits2 = sim_map_dir+sim_map_root[1]+'_'+intfile_ident+'.fits'

            sim_mapfiles = [sim_map_fits1, sim_map_fits2]

            end_to_end_powspec_data, sim_mapfiles, mask_padded, reso_arcmin, $
            banddef=banddef, mapname=mapname, $
            kernel=kernel, ellkern=ellkern, $
            maxell=maxell, $
            beamfiles=beamfile, $
            spectrum=spectrum, $
            invkernmat=invkernmat, $
            npix=npix, workdir=workdir, $
            intfile_ident=intfile_ident, $
            delete_intfile=delete_intfile, $
            resume=resume

            dls = spectrum[1:*,0]
            save, bandcenters, dls, filename=dls_savfile

            print, 'isim = '+sidx+' finished.'
            ;;;;;file_delete, sim_mapfiles
        endelse
        
        if (isim eq 0) then begin
            num_bands = n_elements(dls)
            dls_sims = dblarr(num_bands, num_sims)
        endif

        dls_sims[*,isim] = dls[*]
    endfor
end

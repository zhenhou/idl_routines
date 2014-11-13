pro run_noise_sims_cross, field_idx, covmap_fits_files, num_sims, $
    workpath=workpath, beamfile=beamfile, $
    bandcenters=bandcenters, dls_noise_sims=dls_noise_sims, $
    intfile_ident=intfile_ident, $
    resume=resume

    f = lps12_fieldstruct()
    fst = f[field_idx]
    field_name = fst.name

    reso_arcmin = 1.0
    info = get_lps12_fieldinfo(field_idx)
    npix = info.npix

    ;------------- 
    ; directories 
    ;-------------
    mask_dir    = '/data23/hou/projects/spt_x_planck/lps12/masks_50mJy/'
    kern_dir    = mask_dir
    workdir     = workpath+'/'+field_name+'/'
    res = file_info(workdir)
    if (not res.exists) then spawn, ['mkdir', '-p', workdir], /noshell

    ;-------------------------
    ; get MASK and KERNEL
    ;   restores 'kernel' and 'mask_padded'
    ;-------------------------
    kernfile = kern_dir+'kern_'+field_name+'.sav'
    restore, kernfile
    mask_padded = get_lps12_mask(field_idx, /padded)

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
    
    ;seed = -100L*field_idx
    mapname = 'MAP.MAP'

    sim_map_dir = workdir+'sim_maps/'
    s = file_info(sim_map_dir)
    if (not s.exists) then file_mkdir, sim_map_dir

    sim_dls_dir = workdir+'sim_dls/'
    s = file_info(sim_dls_dir)
    if (not s.exists) then file_mkdir, sim_dls_dir

    for isim=0, num_sims-1 do begin
        sidx = strcompress(string(isim),/remove)
        intfile_ident = 'sim'+sidx

        dls_savfile = sim_dls_dir+'dls_halfring_noise_cross_'+intfile_ident+'.sav'
        res = file_info(dls_savfile)

        if res.exists then begin
            restore, dls_savfile
        endif else begin
            ;-------------
            ; sim maps 
            ;-------------

            sim_map_fits1 = sim_map_dir+'simmap_halfring1_'+intfile_ident+'.fits'
            sim_map_fits2 = sim_map_dir+'simmap_halfring2_'+intfile_ident+'.fits'

            ;gauss_noise_sim_covmap, covmap_fits_files[0], noise_map_scale=1.0d6, seed=seed, output_fits_file=sim_map_fits1
            ;gauss_noise_sim_covmap, covmap_fits_files[1], noise_map_scale=1.0d6, seed=seed, output_fits_file=sim_map_fits2

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
            resume=resume

            dls_noise = spectrum[1:*,0]
            save, bandcenters, dls_noise, filename=dls_savfile

            ;file_delete, sim_mapfiles
        endelse
        
        if (isim eq 0) then begin
            num_bands = n_elements(dls_noise)
            dls_noise_sims = dblarr(num_bands, num_sims)
        endif

        dls_noise_sims[*,isim] = dls_noise[*]
    endfor
end

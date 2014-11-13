pro run_end2end_planck_field, field_idx, $
mapfile_path=mapfile_path, mapfile_root=mapfile_root, $
mapname=mapname, $
beamfile=beamfile, $
invkernmat=invkernmat, $
auto=auto, $
workpath=workpath, $
bandcenters=bandcenters, spectrum=spectrum, resume=resume

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

    ;-------------
    ; planck maps 
    ;-------------
    ;num_maps = n_elements(mapfile_root)
    ;mapfiles = strarr(num_maps/2,2)
    ptsrc_maps_path = mapfile_path+'/'+field_name+'/'
    mapfiles = ptsrc_maps_path+mapfile_root+'.fits'
    if (not keyword_set(mapname)) then mapname = 'MAP.MAP'

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

    ;-------------------------
    ; get BEAMFILES
    ;-------------------------
    beamfiles = beamfile
    ;beamfiles = '/data23/hou/planck_data/2013/beams/hfi_beam_143x143_nominal_wpix_R1.10.txt'

    ;-------------------------
    ; Run end_to_end
    ;-------------------------
    end_to_end_powspec_data, mapfiles, mask_padded, reso_arcmin, $
    banddef=banddef, mapname=mapname, $
    kernel=kernel, ellkern=ellkern, $
    maxell=maxell, $
    beamfiles=beamfiles, $
    auto=auto, $
    invkernmat=invkernmat, $
    spectrum=spectrum, $
    npix=npix, workdir=workdir, $
    resume=resume
    ;curlyw=curlyw, $
    ;raw_spectrum_data=spectrum_data_raw, $
    ;cov_data_raw=cov_data_raw, $
    ;calib=calib
end

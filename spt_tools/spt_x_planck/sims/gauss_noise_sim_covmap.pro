pro gauss_noise_sim_covmap, covmap_fits_file, noise_map, noise_map_scale=noise_map_scale, seed=seed, output_fits_file=output_fits_file
    
    res = file_info(covmap_fits_file)
    if (not res.exists) then begin
        print, covmap_fits_file
        print, "does not exist."
        stop
    endif

    m = read_spt_fits(covmap_fits_file)
    nx0 = n_elements(m.map.map[*,0])
    ny0 = n_elements(m.map.map[0,*])

    nx = m.mapinfo.nsidex
    ny = m.mapinfo.nsidey

    if (nx0 ne nx or ny0 ne ny) then begin
        print, covmap_fits_map
        print, 'internal dimension error.'
        stop
    endif

    cov = double(m.map.map)
    
    if (keyword_set(seed)) then seed0=seed
    rand_gauss = randomn(seed0, nx, ny, /double)
    seed = seed0
    
    if (keyword_set(noise_map_scale)) then scale=noise_map_scale else scale=1.0d0
    noise_map = sqrt(cov) * rand_gauss * double(scale)

    if (keyword_set(output_fits_file)) then begin
        m.map.map = float(noise_map)

        write_processed_fits, m, output_fits_file
        print, 'written - ',output_fits_file
    endif
    
    return
end

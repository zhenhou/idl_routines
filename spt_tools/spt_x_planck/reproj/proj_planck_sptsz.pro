pro array_2dto1d, arr2d, arr1d
    
    nx = n_elements(arr2d[*,0])
    ny = n_elements(arr2d[0,*])

    n = nx * ny
    arr1d = dblarr(n)
    arr1d[0:n-1] = arr2d[0:n-1]

    return
end


pro array_1dto2d, arr1d, arr2d, nx=nx, ny=ny
    dim = [nx, ny]

    n = nx * ny
    if (n_elements(arr1d) ne n) then begin
        print, "dimensions do not match"
        print, n_elements(arr1d), n
        stop
    endif

    loc = lindgen(n)
    ind = ARRAY_INDICES(dim, loc, /DIMENSIONS)

    ipx = ind[0,*]
    ipy = ind[1,*]

    arr2d = dblarr(nx,ny)

    arr2d[ipx,ipy] = arr1d[*]

    return
end


function proj_struct, fields_info, hpx_map_file, max_order, spt_output_path=spt_output_path, spt_output_root=spt_output_root, $
         planck_output_path=planck_output_path, planck_output_root=planck_output_root, $
         jackhalf=jackhalf, half_roots=half_roots, projection=projection, reso_arcmin=reso_arcmin

    num_fields = 20

    tmp = {theta_file: ' ', $
           phi_file: ' ', $
           prj_file: ' ', $
           hpx_map_file: ' ', $
           hpx_nside: 2048L, $
           lmax: 3000L, $
           max_order: 0, $
           projection: 5, $
           reso_arcmin: 1.0, $
           jackhalf: 0, $
           half_files: [' ', ' '] $
    }

    proj = replicate(tmp, num_fields)

    res = getsize_fits(hpx_map_file, nside=nside)

    for i_field=0, num_fields-1 do begin
        file_mkdir, spt_output_path+'/'+fields_info[i_field].name
        file_mkdir, planck_output_path+'/'+fields_info[i_field].name

        proj[i_field].theta_file = spt_output_path+'/'+fields_info[i_field].name+'/'+spt_output_root+'_theta'
        proj[i_field].phi_file = spt_output_path+'/'+fields_info[i_field].name+'/'+spt_output_root+'_phi'
        proj[i_field].prj_file = planck_output_path+'/'+fields_info[i_field].name+'/'+planck_output_root+'_prj'
        proj[i_field].hpx_map_file = hpx_map_file
        proj[i_field].hpx_nside = nside
        proj[i_field].lmax = 2L*nside
        proj[i_field].max_order = max_order
        proj[i_field].projection = projection
        proj[i_field].reso_arcmin = reso_arcmin

        if keyword_set(jackhalf) and jackhalf eq 1 then begin
            proj[i_field].jackhalf = 1
            if (n_elements(half_roots) ne 2) then begin
                print, "half_roots should have 2 elements"
                stop
            endif
            proj[i_field].half_files = planck_output_path+'/'+fields_info[i_field].name+'/'+half_roots+'_prj'
        endif
    endfor
    
    return, proj
end


pro create_proj_ini, ini_file, proj

    num_fields = n_elements(proj[*])
    
    get_lun, unit
    openw, unit, ini_file
    printf, unit, format='(1024A)', 'num_fields = ', strcompress(string(n_elements(proj[*])),/remove)
    printf, unit, format='(1024A)', 'hpx_map_nside = ', strcompress(string(proj[0].hpx_nside),/remove)
    printf, unit, format='(1024A)', 'hpx_map_file  = ', proj[0].hpx_map_file
    printf, unit, format='(1024A)', 'has_pol = F'
    printf, unit, format='(1024A)', 'expand_lmax = ',strcompress(string(proj[0].lmax),/remove)
    printf, unit, format='(1024A)', 'max_expand_order = ',strcompress(string(proj[0].max_order),/remove)
    printf, unit, '## io for each field ##'

    for i_field=0, num_fields-1 do begin
        fstr = strcompress(string(i_field),/remove)
        printf, format='(1024A)', unit, 'gal_theta_file'+fstr+' = '+proj[i_field].theta_file
        printf, format='(1024A)', unit, 'gal_phi_file'+fstr+'   = '+proj[i_field].phi_file
        printf, format='(1024A)', unit, 'out_prj_file'+fstr+'   = '+proj[i_field].prj_file
    endfor
    free_lun, unit
end


pro write_prj_fits, proj, fields_info, spt_freq=spt_freq

    num_fields = n_elements(proj[*])
    for i_field=0, num_fields-1 do begin
        prj_file = proj[i_field].prj_file
        prj_info = file_info(prj_file)
        
        ;if (getenv('HOSTNAME') eq 'spt') then begin
        ;    map_path = fields_info[i_field].xspec_map_dir
        ;    if (keyword_set(spt_freq)) then begin
        ;        fits_files = file_search(map_path, '*_'+strcompress(string(spt_freq),/remove)+'_*.fits')
        ;    endif else begin
        ;        fits_files = file_search(map_path, '*.fits')
        ;    endelse
        ;    res = read_spt_fits(fits_files[0])
        ;endif
    
        if (getenv('HOSTNAME') eq 'midway') then begin
            fits_file = '/home/zhenhou/scratch-data/projects/spt_x_planck/planck_2013/reproj/'+fields_info[i_field].name+'/hfi_SkyMap_143_nominal_ringfull_maxOrder4_prj.fits'
            res = read_spt_fits(fits_file)
        endif

        map_struct = expand_fits_struct(res)
        
        if (prj_info.exists) then begin
            num_pixels = prj_info.size / 8
            prj_data = dblarr(num_pixels)
            
            get_lun, unit
            openr, unit, prj_file
            readu, unit, prj_data
            free_lun, unit

            nx = fields_info[i_field].nx_clt
            ny = fields_info[i_field].ny_clt

            prj_2d = reform(prj_data, nx, ny)
            map_struct.mapinfo.nsidex = nx
            map_struct.mapinfo.nsidey = ny
            map_struct.mapinfo.projection = proj[i_field].projection
            
            map_struct.mapinfo.ra0  = fields_info[i_field].ra0_clt
            map_struct.mapinfo.dec0 = fields_info[i_field].dec0_clt
            map_struct.mapinfo.reso_arcmin = proj[i_field].reso_arcmin

            ;map_struct.map.map = prj_2d
            map_str = create_struct('MAP',prj_2d)

            ;if proj[i_field].jackhalf then begin
            ;    prj_half1 = dblarr(num_pixels)
            ;    prj_half2 = dblarr(num_pixels)
            ;    half1_file = proj[i_field].half_files[0]
            ;    get_lun, unit
            ;    openr, unit, half1_file
            ;    readu, unit, prj_half1
            ;    free_lun, unit

            ;    half2_file = proj[i_field].half_files[1]
            ;    get_lun, unit
            ;    openr, unit, half2_file
            ;    readu, unit, prj_half2
            ;    free_lun, unit

            ;    dhalf = prj_half1 - prj_half2
            ;    array_1dto2d, dhalf, dhalf_2d, $
            ;    nx=map_struct.mapinfo.nsidex, ny=map_struct.mapinfo.nsidey

            ;    map_struct.dmap.map = dhalf_2d
            ;endif else begin
            ;    map_struct.dmap.map[*] = 0.0
            ;endelse

            map_write = create_struct('MAP',map_str, 'MAPINFO',map_struct.mapinfo)

            fits_file = prj_file+'.fits'
            write_processed_fits, map_write, fits_file
            print, 'written - ',fits_file
        endif else begin
            print, 'prj_file: '+prj_file
            print, 'does not exist. skip'
        endelse
    endfor
end


pro proj_planck_sptsz, max_order, planck_map_file=planck_map_file, rewrite_thetaphi=rewrite_thetaphi, $
    spt_output_path=spt_output_path, spt_output_root=spt_output_root, spt_freq=spt_freq, $
    planck_output_path=planck_output_path, planck_output_root=planck_output_root, $
    jackhalf=jackhalf, half_roots=half_roots, projection=projection, reso_arcmin=reso_arcmin

    home = getenv('HOME')

    data_path = home+'/data/'

    if (not keyword_set(spt_output_path)) then spt_output_path = data_path+'projects/spt_x_planck/reproj/'
    if (not keyword_set(planck_output_path)) then planck_output_path = data_path+'projects/spt_x_planck/reproj/'
    if (not keyword_set(spt_freq)) then spt_freq = 150

    if (not keyword_set(planck_map_file)) then $
    planck_map_file=data_path+'planck_data/2014/all_sky_maps/single_field_maps/HFI_SkyMap_RING_143_2048_R1.10_nominal.fits'
    
    fields = spt_fieldsinfo()
    proj   = proj_struct(fields, planck_map_file, max_order, spt_output_path=spt_output_path, spt_output_root=spt_output_root, $
             planck_output_path=planck_output_path, planck_output_root=planck_output_root, $
             jackhalf=jackhalf, half_roots=half_roots, projection=projection, reso_arcmin=reso_arcmin)

    ;; make ini file for fortran program ;;
    ini_file = 'ini_zone/'+planck_output_root+'_PRJTO_'+spt_output_root+'.ini'
    create_proj_ini, ini_file, proj

    num_fields = n_elements(proj[*])

    ;; first run, create dirs 

    ;for i=0, num_fields-1 do begin
    ;    spawn, ['mkdir','-p',output_path+fields.name], /noshell
    ;endfor
    
    prj_exists = 0L
    for i_field=0, num_fields-1 do begin
        theta_file = proj[i_field].theta_file
        phi_file   = proj[i_field].phi_file
        prj_file   = proj[i_field].prj_file
        
        theta_info = file_info(theta_file)
        phi_info   = file_info(phi_file)
        prj_info   = file_info(prj_file)
        if ((not theta_info.exists) or (not phi_info.exists) or keyword_set(rewrite_thetaphi)) then begin
            
            sptsz_field_coord, i_field, theta, phi, coord='G', freq=spt_freq, $
            nx=fields[i_field].nx_clt, ny=fields[i_field].ny_clt, $
            ra0=fields[i_field].ra0_clt, dec0=fields[i_field].dec0_clt, reso_arcmin=reso_arcmin, proj=proj[i_field].projection

            array_2dto1d, theta, theta_1d
            array_2dto1d, phi,   phi_1d

            get_lun, unit
            openw, unit, theta_file
            writeu, unit, theta_1d
            free_lun, unit

            get_lun, unit
            openw, unit, phi_file
            writeu, unit, phi_1d
            free_lun, unit
        endif

        print, "theta, phi written for "+fields[i_field].name

        prj_exists += long(prj_info.exists)
    endfor
    
    hpx_taylor = home+'/Projects/projects/spt_x_planck/proc_maps/reproj/healpix_taylor/hpx_taylor'
    if (prj_exists eq num_fields) then begin
        print, "all prj files ready. skip hpx_taylor."
    endif else begin
        spawn, [hpx_taylor, ini_file], /noshell
    endelse
	

    write_prj_fits, proj, fields, spt_freq=spt_freq
end



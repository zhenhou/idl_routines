pro make_planck_powspec_01, bandcenters, Dls_ave

    f = lps12_fieldstruct()
    
    num_fields = 19

    field_index = [0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]

    Dls_all = fltarr(57,num_fields)
    Dls_ave = fltarr(57)

    weight = dblarr(num_fields)

    for ifield=0, num_fields-1 do begin
        path = '/data/hou/projects/spt_x_planck/powspec/planck/'+f[field_index[ifield]].name+'/data/'
        savfile = path+'Dls_01_autoraw_noise.sav'

        res = file_info(savfile)
        if res.exists then begin
            restore, savfile
        endif else begin
            try_end2end_planck_01, field_index[ifield], banddef=banddef, spectrum=Dls_raw0, /resume
            try_end2end_planck_01, field_index[ifield], banddef=banddef, spectrum=Dls_noise0, /resume, /jackhalf
            
            Dls_raw = Dls_raw0[1:*]
            Dls_noise = Dls_noise0[1:*]

            nbands=n_elements(banddef)-1
            bandcenters=(banddef[0:nbands-1]+banddef[1:nbands])/2

            save, bandcenters, Dls_raw, Dls_noise, filename=savfile
        endelse

        Dls_all[*,ifield] = Dls_raw - 0.25*Dls_noise
        weight[ifield] = double(f[ifield].area)
    endfor

    nm = total(weight)
    weight /= nm

    for i=0, 56 do Dls_ave[i] = total(Dls_all[i,*]*weight[*])

    return
end

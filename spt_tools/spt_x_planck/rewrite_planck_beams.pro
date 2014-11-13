pro rewrite_planck_beams, convol_wpix=convol_wpix, nominal_only=nominal_only
    
    fits_file = '/data23/hou/planck_data/2013/RIMO/HFI_RIMO_R1.10.fits'
    output_path = '/data23/hou/planck_data/2013/beams/'
    
    freqs = [100, 143, 217, 353, 545, 857]
    
    if (keyword_set(wpix)) then begin
        wpix = fltarr(8001,1)
        wpix[*] = 1.00e0
    endif else begin
        wpix = healpixwindow(2048)
    endelse

    extno_start = 8
    
    extno_ct = 0
    for ib1=0, 5 do begin
        for ib2=ib1, 5 do begin
            str = strcompress(string(freqs[ib1]),/remove)+'x'+strcompress(string(freqs[ib2]),/remove)

            str_add = str

            if (keyword_set(nominal_only)) then str_add+='_nominal'
            
            res=mrdfits(fits_file, extno_start+extno_ct, hdr)
            n = n_elements(hdr)
            pos = -1
            for i=0, n-1 do begin
                pos = strpos(hdr[i], str)
                if pos ne -1 then break
            endfor
            
            if pos eq -1 then begin
                print, "something wrong while reading"
                print, hdr
                stop
            endif
            
            if (keyword_set(convol_wpix)) then $
            output_file = output_path+'hfi_beam_'+str_add+'_wpix_R1.10.txt' $
            else $
            output_file = output_path+'hfi_beam_'+str_add+'_R1.10.txt'
            openw, 5, output_file
            for il=0L, 4000L do begin
                if (keyword_set(nominal_only)) then begin
                    printf, 5, format='(I6,E16.7)', il, (res.nominal)[il]*wpix[il,0]
                endif else begin
                    printf, 5, format='(I6,6E16.7)', il, (res.nominal)[il]*wpix[il,0], (res.EIGEN_1)[il], (res.EIGEN_2)[il], (res.EIGEN_3)[il], (res.EIGEN_4)[il], (res.EIGEN_5)[il]
                endelse
            endfor
            close, 5

            extno_ct += 1
        endfor
    endfor
end

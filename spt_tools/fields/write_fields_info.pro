pro write_fields_info

    f = lps12_fieldstruct()
    num_fields = 20

    openw, 5, 'sptsz_fields_info.txt'
    printf, 5, format='(A-20,A10,A10,A10,A10,A10)', '#field_name', 'RA0', 'DEC0', 'NSIDEX', 'NSIDEY', 'RESO_ARCMIN'
    for i=0, num_fields-1 do begin
        fits_file='/project/kicp/zhenhou/projects/spt_x_planck/reproj/'+f[i].name+'/hfi_143_nominal_maxOrder0_prj.fits'
        res = read_spt_fits(fits_file)
        printf, 5, format='(A-20,F10.4,F10.4,I10,I10,F10.6)', f[i].name, f[i].RA0, f[i].DEC0, res.mapinfo.nsidex, res.mapinfo.nsidey,res.mapinfo.reso_arcmin
    endfor
    close, 5
end

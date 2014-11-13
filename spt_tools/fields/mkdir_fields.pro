pro mkdir_fields, path, num_fields=num_fields

    res = file_info(path)
    if (not res.exists) then begin
        print, 'do "mkdir '+path+'" first!'
        stop
    endif

    f = lps12_fieldstruct()

    if (keyword_set(num_fields)) then nfields = num_fields else nfields = 20

    for ifield=0, nfields-1 do begin
        field_name = f[ifield].name
        spawn, ['mkdir',path+'/'+field_name], /noshell
    endfor

end

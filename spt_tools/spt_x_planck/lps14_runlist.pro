pro lps14_runlist, ifield, list, runlist=runlist

    runlist_path = '/home/zhenhou/scratch-data/projects/sptsz_lowl/runlists/'
    f = lps12_fieldstruct()
    field_name = f[ifield].name

    ;print, field_name

    runlist = runlist_path+'runlist_lps14_'+field_name+'_0602.txt'
    nlines = file_lines(runlist)

    tmp = strarr(nlines)
    openr, 5, runlist
    readf, 5, tmp
    close, 5

    ncolumns = n_elements(strsplit(tmp[0],/extract))
    list_tmp = strarr(nlines, ncolumns)

    ;print, 'nlines = '+strcompress(string(nlines),/remove)+', ncolumns = '+strcompress(string(ncolumns),/remove)
    
    for i=0, nlines-1 do begin
        list_tmp[i,*] = strsplit(tmp[i],/extract)
    endfor

    list = strarr(nlines*ncolumns)
    for ic=0, ncolumns-1 do begin
        list[ic*nlines:(ic+1)*nlines-1] = list_tmp[0:nlines-1,ic]
    endfor

    return
end

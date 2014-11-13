function coadd_signalonly_dir,dir,fstub,nowt=nowt

list = file_search(dir+fstub+'*',count=nl)

if nl le 0 then stop

info =  spt_read_binary_header(list[0])

if info.isdbl eq 1 then stop ; not yet handled
;stop

coadds=fltarr(info.nx,info.ny,info.nmap)
wts=fltarr(info.nx,info.ny)
tmap=fltarr(info.nx,info.ny)
twt=tmap

nmappix = long64(n_elements(tmap))

for i=0,nl-1 do begin
    twt = spt_read_binary_wt(list[i])
    if n_elements(twt) lt nmappix then $
      pting = spt_read_binary_pting(list[i]) $
    else $
      pting=lindgen(nmappix)

    if keyword_set(nowt) then begin
        ind = where(twt gt 0)
        twt(ind)=1.0
    endif
    wts(pting)+=twt

    for j=0,info.nmap-1 do begin
        tmap=spt_read_binary_mapi(list[i],j)
        coadds[pting+j*nmappix]+=tmap*twt
    endfor
    

endfor

ind = where(wts gt 0)
nn=info.nx*long64(info.ny)
for j=0,info.nmap-1 do coadds[ind+nn*j]/=wts(ind)
return,{maps:coadds,wt:wts}

end

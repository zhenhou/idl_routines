function spt_read_binary_header,file

ndims = 8

dims = lonarr(8)


openr,lun,file,/get_lun
readu,lun,dims
free_lun,lun
;stop
isdbl = dims[3]
nmap=dims[2]
nx=dims[0]
ny=dims[1]
havejack=dims[4]
version=dims[5]
nmappix=dims[6]
npix=dims[7]
return,{nx:nx,ny:ny,nmap:nmap,isdbl:isdbl,havejack:havejack,version:version,npix:npix,nmappix:nmappix}

end

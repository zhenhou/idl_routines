function spt_read_binary_mapi,file,i,forceversion=forceversion

if n_elements(forceversion) ne 1 then forceversion = -1

ndims = 8
if forceversion eq 0 then $
  ndims = 4

dims = lonarr(8)

openr,lun,file,/get_lun
readu,lun,dims
if dims[5] lt 1 then ndims = 4;
if dims[5] lt 4 then ndims = 6;


if i ge dims[2] then begin
free_lun,lun
return,-1
endif

isdbl = dims[3]
;isdmap = dims[4]
;print,dims[0:2]
;version = dims[5]

sizeflt = long64(4)
sizedbl = long64(8)
sizelon = long64(4)
sizel64 = long64(8)
sizeint = long64(2)

offset = long64(sizelon*ndims)

mapsize=dims[6]
version = dims[5]
if forceversion eq 0 then version = 0

if version lt 4 then mapsize = dims[0]*dims[1]
fmapsize=dims[0]*dims[1]
if fmapsize eq mapsize then begin
    if isdbl then begin
        offset += sizedbl*mapsize*i
        map = dblarr(dims[0],dims[1])
    endif else begin
        offset += sizeflt*mapsize*i
        map = fltarr(dims[0],dims[1])
    endelse
endif else begin
    if isdbl then begin
        offset += sizedbl*mapsize*i
        map = dblarr(mapsize)
    endif else begin
        offset += sizeflt*mapsize*i
        map = fltarr(mapsize)
    endelse
endelse
point_lun,lun,offset

readu,lun,map
free_lun,lun
return,map
end

function spt_read_binary_wt,file,forceversion=forceversion

if n_elements(forceversion) ne 1 then forceversion = -1

ndims = 8
if forceversion eq 0 then ndims = 4

dims = lonarr(8)
openr,lun,file,/get_lun
readu,lun,dims

;readu,lun,dims
if dims[5] lt 1 then ndims = 4
if dims[5] lt 4 then ndims = 6  ;     
if forceversion eq 0 then ndims = 4

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
        offset += sizedbl*mapsize*dims[2]
        wt = dblarr(dims[0],dims[1])
    endif else begin
        offset += sizeflt*mapsize*dims[2]
        wt = fltarr(dims[0],dims[1])
    endelse
endif else begin
    if isdbl then begin
        offset += sizedbl*mapsize*dims[2]
        wt = dblarr(mapsize)
    endif else begin
        offset += sizeflt*mapsize*dims[2]
        wt = fltarr(mapsize)
    endelse
endelse

point_lun,lun,offset
    
readu,lun,wt
free_lun,lun
return,wt
end

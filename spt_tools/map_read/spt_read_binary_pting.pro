function spt_read_binary_pting,file,forceversion=forceversion,n1=n1,n2=n2

if n_elements(forceversion) ne 1 then forceversion = -1

dims = lonarr(8)

openr,lun,file,/get_lun
readu,lun,dims
if dims[5] lt 4 or forceversion eq 0 or dims[7] eq 0 then begin
    free_lun,lun
    return,-1
endif

n1=dims[0]
n2=dims[1]
npixarr=dims[7]
ndims=8

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

if isdbl then begin
    offset += sizedbl*mapsize*(dims[2]+1)
    
endif else begin
    offset += sizeflt*mapsize*(dims[2]+1)
endelse
tmp = lonarr(npixarr)
point_lun,lun,offset
readu,lun,tmp
free_lun,lun
if dims[6] eq npixarr then return,tmp
pting = lonarr(dims[6])
k=0l
npair = npixarr/2
for i=0l,npair-1 do begin
    for j=tmp[2*i],tmp[2*i+1] do $
      pting[k++]=j
endfor

return,pting
end

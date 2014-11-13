function spt_check_offset,file,poly

sizeflt = long64(4)
sizedbl = long64(8)
sizelon = long64(4)
sizel64 = long64(8)
sizeint = long64(2)


points = fltarr(50)
openr,lun,file,/get_lun
dims = lonarr(3)
readu,lun,dims
n1 = dims[0]
n2 = dims[1]
nsim=dims[2]
offset = long64(sizelon*long64(n1)*long64(n2)*long64(nsim+1))
point_lun,lun,offset

readu,lun,points
;print,points[4:10]
ind = where(points eq 19)

free_lun,lun
return,ind[0]-2
end

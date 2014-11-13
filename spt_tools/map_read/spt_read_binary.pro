function spt_read_binary,file
;print,'check before using!'
;stop
;

dims = lonarr(6)
openr,lun,file,/get_lun
readu,lun,dims

isdbl = dims[3]
isdmap = dims[4]
print,dims[0:2]
version = dims[5]

;nobs = lon64arr(dims[0],dims[1])
if isdbl then begin
    map =dblarr(dims[0],dims[1],dims[2])
    wt = dblarr(dims[0],dims[1])
 ;   if usedmap then    dmap = map
endif else begin
    map =fltarr(dims[0],dims[1],dims[2])
    wt = fltarr(dims[0],dims[1])
;    if usedmap then    dmap = map
endelse
readu,lun,map
readu,lun,wt


params = fltarr(7)
readu,lun,params
print,params
fpdim = lon64arr(2)
readu,lun,fpdim
print,fpdim
bolowt = dblarr(fpdim[0]);
bolocal = fltarr(fpdim[0]);
readu,lun,bolowt
readu,lun,bolocal
print,'c'
flags = lonarr(fpdim[0],fpdim[1]);
readu,lun,flags
free_lun,lun
print,'ret'

return,{map:map,wt:wt,nbolo:fpdim[0],nscan:fpdim[1],bolowt:bolowt,bolocal:bolocal,flags:flags,params:params,desc:dims}

end

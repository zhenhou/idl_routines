function spt_check_version,file


openr,lun,file,/get_lun
dims = lonarr(6)
readu,lun,dims
free_lun,lun
return,dims[5]
end

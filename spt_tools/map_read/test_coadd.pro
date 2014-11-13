pro test_coadd

files = ['~/data/map_ra5h30dec-55_150_20080229_112837_isims-0-74.dat',$
         '~/data/map_ra5h30dec-55_150_20080229_081332_isims-0-74.dat']

hdr = spt_read_binary_header(files[0])

;first get the franklin map
maps = fltarr(hdr.nx,hdr.ny,2)
openr,lun,/get_lun,'~/data/test.dat'
readu,lun,maps
free_lun,lun

;next read in the input maps
p1=spt_read_binary_pting(files[0])
wt1=spt_read_binary_wt(files[0])
uwt1 = wt1*0
uwt1(where(wt1 gt 0))=1.0
wt2=spt_read_binary_wt(files[1])
uwt2 = wt1*0
uwt2(where(wt2 gt 0))=1.0

m1a = spt_read_binary_mapi(files[0],0)
m1b = spt_read_binary_mapi(files[0],1)
m2a = spt_read_binary_mapi(files[1],0)
m2b = spt_read_binary_mapi(files[1],1)

ma = m1a+m2a
mb = m1b+m2b
wt = uwt1+uwt2
ind = where(wt gt 0)
ma(ind) /= wt(ind)
mb(ind) /= wt(ind)

omaps = maps*0.0
omaps[p1]=ma
omaps[p1+hdr.nx*hdr.ny]=mb
stop
end

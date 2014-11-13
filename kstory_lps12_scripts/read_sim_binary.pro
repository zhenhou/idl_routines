;;;
; script to read binary sim file
; 20120511
;;;
FUNCTION read_sim_binary, fname


hdr = spt_read_binary_header(fname)
nsim=long64(hdr.nmap)
offset = long64(4*8)
npix = hdr.nx*long64(hdr.ny)
tmap = fltarr(npix,hdr.nmap)

pting = spt_read_binary_pting(fname)
wt = spt_read_binary_wt(fname)
;twt[pting]+=wt
npp=n_elements(pting)
scr = fltarr(npp)
get_lun,lun
openr,lun,fname
point_lun,lun,offset
for j=0,nsim-1 do begin
    readu,lun,scr
    tmap[pting,j]+=scr*wt
endfor
close,lun        


RETURN, reform(tmap, hdr.nx, hdr.ny, nsim)
END

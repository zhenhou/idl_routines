;;;
; function to read in a .bin file of maps
; 06/07/2012: (KTS) created
; 03/07/2013: (KTS) add nowt option
;;;
FUNCTION read_bin_map, fname, nowt=nowt


hdr = spt_read_binary_header(fname)
nsim=long64(hdr.nmap)
offset = long64(4*8)
npix = hdr.nx*long64(hdr.ny)
tmap = fltarr(npix,hdr.nmap)

pting = spt_read_binary_pting(fname)
wt = spt_read_binary_wt(fname)
if keyword_set(nowt) then wt = 1.
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
close,/all


RETURN, reform(tmap, hdr.nx, hdr.ny, nsim)
END

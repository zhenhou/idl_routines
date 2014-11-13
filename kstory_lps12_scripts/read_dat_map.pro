;;;
; function to read in a .dat file of maps
; 04/25/2012: (KTS) created
;;;

FUNCTION read_dat_map, fname, npix_x, npix_y, nmaps

maps = fltarr(npix_x,npix_y, nmaps)

get_lun,u
openr,u,fname
readu,u,maps
free_lun,u

RETURN, maps
END

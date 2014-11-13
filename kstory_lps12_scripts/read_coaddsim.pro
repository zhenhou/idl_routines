;;;
; script for oliver
; 20120511
;;;
FUNCTION read_coaddsim, idx, fname

f = lps12_fieldstruct()
field = f[idx].name

info = get_lps12_fieldinfo(idx) 
nx = info.npix[0] & ny = info.npix[1]
n_sim_coadds = 100

;if ~keyword_set(fname) then fname = bdir+'coaddsim_'+stub+'_'+field+'.dat'
maps = read_dat_map(fname, nx, ny, n_sim_coadds)

RETURN, maps
END

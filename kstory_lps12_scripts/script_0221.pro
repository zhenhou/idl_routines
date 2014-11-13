;;;
; NAME: script_0221
; PURPOSE:
;   investigate apod masks
;
; INPUTS:
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/21/2012: (KTS) Created
;;;


;------------------------
; test noise psds
pro run_0221

test_coadd = '/data17/rkeisler/ps09/post_29Apr2008_coadd_ra5h30dec-55_150.fits'
d = read_spt_fits(test_coadd)

sky = d.map.map
noise = d.dmap.map

n = long(n_elements(sky[*,0]))
n2 = n*n

ss = shift(abs(fft(sky,1))^2, n/2.,n/2.)
ns = shift(abs(fft(noise,1))^2, n/2.,n/2.)
;ns = shift(abs(fft(noise*apod_mask,1))^2, n/2.,n/2.)


; Make ell-grid
reso_rad = d.mapinfo.reso_arcmin/60.*!dtor
lgrid = shift(make_fft_grid(reso_rad/(2.*!pi),n,n,fx=lx, fy=ly), n/2.,n/2.)
lx = shift(lx, n/2.,0) & ly = shift(ly, 0, n/2.)

; Define binning of final spectrum
lbinsize = 50.
lmax = 10000.
nbin = floor(lmax/lbinsize)
lbin = findgen(nbin)*lbinsize + 0.
cbin = lbin + lbin[1]/2. ; use bin centers for calculation

;------------------------------
; Calculate Cl
;------------------------------
print, '   *** calculate cl'

cl   = dblarr(nbin) ; spectrum
cl_n = dblarr(nbin) ; noise spectrum

for il=1, nbin-1 do begin
    wh = where( lgrid ge lbin[il] and $
                lgrid lt (lbin[il] + lbinsize), nwh)

    b = where( lgrid ge lbin[il] and lgrid lt (lbin[il] + lbinsize) )

    cl[il]    = total(ss[b])
    cl_n[il]  = total(ns[b])
endfor

wset, 1
plot, cbin, cl, xrange=[0.,lmax],yrange=[10.,10000.],/ylog,psym=3

wset, 2
plot, cbin, cl_n, xrange=[0.,lmax],yrange=[10.,10000.],/ylog,psym=3

end


;;;
; NAME: make_wtrms_coadds
; PURPOSE:
;   Calculate dl at 3000 for coadded dmaps
;
; INPUTS:
;    field_name,     name of field
;    cut_idx_list,   array of cuts to apply, from cut_set()
;
; MODIFICATION HISTORY:
;  08/21/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro calc_dmap_ps_0821, field_name, cut_idx_list
compile_opt IDL2, HIDDEN

for cc=0, n_elements(cut_idx_list)-1 do begin
    calc_dmap_ps_single_field, field_name, cut_idx_list[cc]
endfor

end

;...................................................................
; Calculate the Dmap PS for a single coadd file
;
pro calc_dmap_ps_single_field, field_name, cut_idx
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'
cutSetNames = cut_set_names()

print, "Read fits file"
field_idx=2
field_name = all_fields[field_idx]
cut_name = cutSetNames[cut_idx]

coadd_fname = '/data/kstory/projects/lps12/runlists/coadds/coadd_wtrms_'+cut_name+'_'+field_name+'.fits'

print, 'Reading coadd file: '+coadd_fname
data = read_spt_fits(coadd_fname)
dmap = data.dmap.map

nsidex = data.mapinfo.nsidex
nsidey = data.mapinfo.nsidey

;------------------------
; Apply Window to middle thrid of the map
print, "Make Hanning Window"

; Get inner third of the map
meshgrid, nsidex, nsidey, xx, yy
xmax = round(nsidex*(2./3)) & xmin = round(nsidex*(1./3))
n_sm = xmax-xmin
;cmap = intarr(nsidex, nsidey)
;ymax = round(nsidey*(2./3)) & ymin = round(nsidey*(1./3))
;wh = where( (xx lt xmax) and (xx gt xmin) and (yy lt xmax) and (yy gt xmin) , nwh)
;cmap[wh] = 1

h = hanning(n_sm)#(fltarr(n_sm)+1.d0) & h*=transpose(h); Hanning window
tmp = DBLARR(nsidex, nsidey)
for ii=0, n_sm-1 do begin
    tmp[0:n_sm-1,ii] = h[*,ii]
endfor
han = shift(tmp, n_sm, n_sm)

dmap_h = dmap*han ; Hanning-windowed data

;------------------------
; Make ell-grid
reso_rad = data.mapinfo.reso_arcmin/60.*!dtor
lgrid = shift(make_fft_grid(reso_rad/(2.*!pi),nsidex,nsidey,fx=lx, fy=ly), nsidex/2.,nsidey/2.)
lx = shift(lx, nsidex/2.,0) & ly = shift(ly, 0, nsidey/2.)

; Define binning of final spectrum
lbinsize = 50.
lmax = 10000.
nbin = floor(lmax/lbinsize)
lbin = findgen(nbin)*lbinsize + 0.
cbin = lbin + lbin[1]/2.

; Make k-space weight to ignore low kx modes
kweight = dblarr(nsidex,nsidey)+1.
kx_cut = 600.
kweight[where(abs(lx) lt kx_cut)] = 0.

;------------------------
; FFT
print, 'Calculate FFT'
ss = shift( abs( fft(dmap_h,1))^2., nsidex/2., nsidey/2.)

sss = ss * kweight

;------------------------
; Calculate Cl
print, 'Calculate Cl'
cl   = dblarr(nbin)

for il=1, nbin-1 do begin
    wh = where( lgrid ge lbin[il] and $
                lgrid lt (lbin[il] + lbinsize) )
    cl[il]    = total(sss[wh]) / total(kweight[wh])
endfor

; Convert to dl
dl    = cl*cbin*(cbin+1)/(2*!pi)


;------------------------
; Plots
;------------------------
wset, 3
plot, cbin,dl, xrange=[0.,lmax],yrange=[10.,10000.],/ylog,psym=3


;------------------------
; print out cl at 700, 3000
print, "*** Field: "+field_name
print, "    cut_name: "+cut_name
print, "    Dl(700)  = ", dl[ where(lbin eq 700) ]
print, "    Dl(3000) = ", dl[ where(lbin eq 3000) ]

stop
;save, file_name='dmap_ps_0821.sav'

end

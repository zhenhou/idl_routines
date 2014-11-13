;;;
; NAME: calc_dmap_ps
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
pro calc_dmap_ps, field_name, cut_idx_list
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
cut_name = cutSetNames[cut_idx]

print, "Read fits file"

; Need special case for ra21hdec-60 because it was made before 
;   I always made first_half, second_half coadds
;if (field_name eq 'ra21hdec-60' ) then begin
if (0) then begin

    coadd_fname = '/data/kstory/projects/lps12/runlists/coadds/'+field_name+'/coadd_wtrms_'+cut_name+'_'+field_name+'.fits'

    print, 'Reading coadd file: '+coadd_fname
    data1 = read_spt_fits(coadd_fname)
    dmap = data1.dmap.map

; All other fields
endif else begin
    coadd_fname1 = '/data/kstory/projects/lps12/runlists/coadds/'+field_name+'/coadd_wtrms_'+cutSetNames[cut_idx]+'_'+field_name+'_first_half.fits'
    coadd_fname2 = '/data/kstory/projects/lps12/runlists/coadds/'+field_name+'/coadd_wtrms_'+cutSetNames[cut_idx]+'_'+field_name+'_second_half.fits'

    print, 'Reading coadd file: '+coadd_fname1
    data1 = read_spt_fits(coadd_fname1)
    data2 = read_spt_fits(coadd_fname2)

    dmap = data1.map.map - data2.map.map
endelse
    
tv_spt_map, dmap, winnum=20

nsidex = data1.mapinfo.nsidex
nsidey = data1.mapinfo.nsidey

;------------------------
; Apply Window to middle thrid of the map
print, "Make Hanning Window"

; Get inner third of the map
meshgrid, nsidex, nsidey, xx, yy
if (nsidex lt nsidey) then begin
    xmax = round(nsidex*(7./8)) & xmin = round(nsidex*(1./8))
endif else begin
    xmax = round(nsidey*(7./8)) & xmin = round(nsidey*(1./8))
endelse
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
han = shift(tmp, nsidex/2 - n_sm/2, nsidey/2 - n_sm/2)

dmap_h = dmap*han ; Hanning-windowed data

tv_spt_map, dmap_h, winnum=21

;------------------------
; Make ell-grid
reso_rad = data1.mapinfo.reso_arcmin/60.*!dtor
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

tv_spt_map, sss, /norms, winnum=22

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
window, 3
plot, cbin,dl, xrange=[0.,lmax],yrange=[100.,100000.],/ylog,psym=3


;------------------------
; print out cl at 700, 3000
print, "*** Field: "+field_name
print, "    cut_name: "+cut_name
print, "    Dl(3000) = ", dl[ where(lbin eq 3000) ], " -------- Dl(2000) = ", dl[ where(lbin eq 2000) ], " -------- Dl(700) = ", dl[ where(lbin eq 700)]

;save, file_name='dmap_ps_0821.sav'
end

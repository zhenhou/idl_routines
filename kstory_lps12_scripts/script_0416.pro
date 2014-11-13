;;;
; NAME: script_0416
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) plot_kmap_for_wiki, plot kmaps

; MODIFICATION HISTORY:
;  04/16/2012: (KTS) Created
;;;

;...................................................................
; plot kmaps
PRO plot_kmap_for_wiki, idx, save_plots=save_plots, $
               stopit=stopit
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()

coadd_dir = '/data/kstory/projects/lps12/maps/20120305/coadds/'
field_name = f[idx].name
file = coadd_dir + 'coadd_'+field_name+'_50mJy.fits'
d = krf(file)

mymap  = d.map.map * 1e6
mydmap = d.dmap.map * 1e6

; fft's
s = size(mymap)
nx = s[1] & ny = s[2]
myfft  = shift(alog(abs(fft(mymap))), nx/2, ny/2)
mydfft = shift(alog(abs(fft(mydmap))), nx/2, ny/2)

; make ellgrid
reso_rad = 1.0/60.*!dtor
lgrid = shift(make_fft_grid(reso_rad/(2.*!pi),nx,ny,fx=lx, fy=ly), nx/2.,ny/2.)
lx = shift(lx, nx/2.,0) & ly = shift(ly, 0, ny/2.)
lxres = lx[2] - lx[1]
lyres = ly[2] - ly[1]

; plots
wh = where( (abs(lx) lt 3500) and (abs(ly) lt 3500) )

whx = where( abs(lx[*,0]) lt 3500 )
why = where( abs(ly[0,*]) lt 3500 )
nx1 = n_elements(whx)
ny1 = n_elements(why)

print, float(ny1)/nx1 * 620
;stop

ss  = myfft[whx[0]:whx[nx1-1], why[0]:why[ny1-1]]
ssd = mydfft[whx[0]:whx[nx1-1], why[0]:why[ny1-1]]

tv_spt_map, ss, resolution=lxres, xtitle='Lx', ytitle='wrong scale', scale=2, /norms, /forcesize, title=field_name+', Kmap', winnum=5, colt=39
tv_spt_map, ssd, resolution=lxres, xtitle='Lx', ytitle='wrong scale', scale=2, /norms, /forcesize, title=field_name+', K_dmap', winnum=6, colt=39


if keyword_set(save_plots) then begin
    wset, 5
    figdir = '/home/kstory/lps12/figs/'
    filename = figdir+'kmap_0416_'+field_name
    err = tvread(/png, filename=filename, /nodialog)

    wset, 6
    figdir = '/home/kstory/lps12/figs/'
    filename = figdir+'kdmap_0416_'+field_name
    err = tvread(/png, filename=filename, /nodialog)
endif
; plot

if keyword_set(stopit) then stop
END

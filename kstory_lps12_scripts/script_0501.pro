;;;
; NAME: script_0501
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) wiki_1,  make plots to check noise
; 2) old_kweight, run old kweight code and plot it.
; 3) process_field, hold-over from 0430
;
; MODIFICATION HISTORY:
;  05/01/2012: (KTS) Created
;;;

;...................................................................
; Check inputs to kweight calculation / noise psd normalization.
PRO wiki_1

field_name='ra5h30dec-55_2008'
nbig = 4320
l = make_fft_grid(1./60.*!dtor,nbig,nbig)*2.*!pi

; cl_th
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
cl_th = interpol(cl_uK2, l_vec, l)

; noise
restore, '/data/kstory/projects/lps12/noise_psds/noise_psd_'+field_name+'.sav'
; psd is in units of K-rad
calib = 0.76 ;
psd_uK = psd * 1d6 * calib
noise = ( psd_uK )^2

; get the TF for this field
restore, '/data/kstory/projects/lps12/twod_tfs/tf_'+field_name+'.sav'
tf = TF_W_BEAM

; find where tf is not NaN, and tf ne 0 (can't divide by 0)
whgood = where( (finite(tf) ne 0) and (tf ne 0), nwh, complement=whbad)
if (nwh ne 0) then begin
    ntmp = dblarr(nbig, nbig)
    ntmp[whgood] = noise[whgood] / tf[whgood]
endif else begin
    ntmp = noise/tf
endelse
noise_tf = ntmp

; get coadded DMAP
d = krf('/home/kstory/lps12/maps/20120420/coadds/coadd_ra5h30dec-55_2008_50mJy.fits')
nsmall = 960
reso = 1.0 & reso_rad = reso/60.*!dtor
dmap = dblarr(nbig, nbig)
dmap[0:nsmall-1, 0:nsmall-1] = d.dmap.map
mask = get_lps12_mask(0, /padded)
dpower = abs(fft(dmap*mask))^2.
dpower *= 1d12 / mean(mask^2.) * nbig^2. * reso_rad^2. * calib^2.
dpower_uK = sqrt(dpower)

; psd from cluster finding
cluster_file = '/data15/sptdat/run2_2010/ra5h30dec-55/coadds/psd_ra5h30dec-55_150.sav'
restore, cluster_file
l_cl = make_fft_grid(0.25/60.*!dtor,3120,3120)*2.*!pi
power_cl = (psd * 1e6 * calib)^2 ; (uK-rad)^2

; Plot 1
window, 5, xsize=700, ysize=500
vec = indgen(1000)
plot, l[vec,vec], cl_th[vec,vec], /ylog, xtitle='ell', ytitle='Cl [uK-rad]^2', title='ra5h30dec-55_2008', xra=[0,4000]
oplot, l[vec,vec], noise[vec,vec], color=!red
oplot, l[vec,vec], dpower[vec,vec], color=!blue
oplot, l[vec,vec], noise_tf[vec,vec], color=!orange
oplot, l_cl[vec,vec], power_cl[vec,vec], color=!green
legend_str = ['Cl_th', 'coadded DMAP', 'noise_psd', 'noise_TF', 'cluster PSD']
legend, legend_str, linestyle=0, colors=[!white, !blue, !red, !orange, !green], position=[2500, 0.9e2]
;err = tvread(/png, filename='/home/kstory/public_html/notebook/spt_lps12/noise_psd_0501', /nodialog)

stop
END


; Check old code
PRO old_kweight
field = 'ra5h30dec-55'
nbig = 2160 & reso = 1.0
l = make_fft_grid(reso/60.*!dtor,nbig,nbig)*2.*!pi

; get mask
    restore,'/data/rkeisler/ps09/mask_'+field+'_20101009_055240.sav'
; get bigmask (zero-padded mask)
    bigmask = fltarr(nbig,nbig)
    bigmask[0:n_elements(mask[*,0])-1, 0:n_elements(mask[0,*])-1] = mask

; get data coadd (this is the MAP (left+right), and DMAP
; (left-right)).  we'll use the DMAP to estimate the noise power.
    restore,'/data/rkeisler/ps09/1.39/coadd_'+field+'.sav'
    bigdmap = fltarr(nbig,nbig)
    bigdmap[0:n_elements(dmap[*,0])-1, 0:n_elements(dmap[0,*])-1] = dmap
    big_dmap_power = abs(fft(bigdmap*bigmask))^2.

    noise = big_dmap_power

readcol,'/data/rkeisler/low_ell_sims/input/dl_input_true_20101221_061503.txt',l_vec,dl_uK2
cl_uK2 = dl_uK2*2.*!pi/l_vec/(l_vec+1.)
cl_uK2[0] = cl_uK2[1]       ; get rid of infinity from dividing by l=0
cl_K2 = cl_uK2*(1d-12)
cl_th = interpol(cl_K2, l_vec, l)

; get the TF for this field
    restore, '/data/rkeisler/ps09/twod_tfs/tf_pad_'+field+'.sav'

    tf = BIG_TF_W_BEAM
    whgood = where( (finite(tf) ne 0) and (tf ne 0), nwh) 
    ntmp = dblarr(nbig, nbig)
    ntmp[whgood] = noise[whgood] / tf[whgood]

!p.multi=[0,1,2]
window, 5, xsize=600, ysize=750
plot, l[0:1000], cl_th[0:1000], /ylog, color=!red, xtitle='ell', ytitle='OLD: (uK-rad)^2'
oplot, l[0:1000], ntmp[0:1000]
plot, l[0:1000], ntmp[0:1000], /ylog
!p.multi=0

stop
END



; Process a field
PRO process_field, idx
coadd_maps_0420, idx
make_mask_lps12, field_idx=idx
make_coupling_kernel, idx
;make_noise_psd_lps12, idx
;try_end2end, idx ??
;make_twod_tfs, idx ??
;make_twod_kweights_lps12, idx ??
END



;...................................................................
; Plot all inputs to the kweight calculation.
PRO kweight_plots

save_files = 0
save_plots = 0
; Set the FWHM for Gaussian smoothing:
l_fwhm=1000.
calib = 0.76 ; use 2008 value
field_idx=0

;------------------------
; Setup
;------------------------
f = lps12_fieldstruct()
field_name = f[field_idx].name
info = get_lps12_fieldinfo(field_idx)
nbig = info.nbig

reso = 1. ; arcminutes per pixel in data maps

; get the l-grid
l = make_fft_grid(reso/60.*!dtor,nbig,nbig)*2.*!pi
lbinsize = l[1] - l[0]
nbin = round(5000./lbinsize)
lbin = findgen(nbin)*lbinsize + 0.


lbin = findgen(nbin)*lbinsize + 0.

; precompute 2d indices for the l-bins
for j=0,nbin-2 do begin
    ex=execute('wh'+strtrim(floor(j),2)+'=where(l ge lbin[j] and l lt lbin[j+1],nwh)')
    ex=execute('nwh'+strtrim(floor(j),2)+'=nwh')
endfor

; get the average power spectrum from the simulated, noise-free observations.
;  Relevant variables: s#.power
;restore,'/data/rkeisler/ps09/create_twodim_tfs.sav'
;restore, '/home/kstory/lps12/scripts/make_twod_tfs.sav'

;------------------------
; Get Inputs
;------------------------

; Get Theoretical Cl's
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
cl_th = interpol(cl_uK2, l_vec, l)

; get noise psd
restore, '/data/kstory/projects/lps12/noise_psds/noise_psd_'+field_name+'.sav'
; psd is in units of K-rad
psd_uK = psd * 1d6
noise = ( psd_uK )^2


; get the TF for this field
restore, '/data/kstory/projects/lps12/twod_tfs/tf_'+field_name+'.sav'
tf = TF_W_BEAM

;----------------------
; calculate the weight_2d
;----------------------
print, 'Gathered all inputs.  Now make the weight_2d'
; find where tf is not NaN, and tf ne 0 (can't divide by 0)
whgood = where( (finite(tf) ne 0) and (tf ne 0), nwh) 
if (nwh ne 0) then begin
    ntmp = dblarr(nbig, nbig)
    ntmp[whgood] = noise[whgood] / tf[whgood]
    weight_2d = ( gauss_smooth(cl_th,l_fwhm,10.) + gauss_smooth(ntmp,l_fwhm,10.) )^(-2.)
endif else begin
    weight_2d = ( gauss_smooth(cl_th,l_fwhm,10.) + gauss_smooth(noise/tf,l_fwhm,10.) )^(-2.)
endelse

w1 = weight_2d ; unnormalized

; renormalize the weight within each l-ring, (l-bin), i.e. divide by
; the maximum weight within that bin.
for j=0,nbin-2 do begin
    ex=execute('wh=wh'+strtrim(floor(j),2))
    ex=execute('nwh=nwh'+strtrim(floor(j),2))
    if nwh gt 0 then begin
        weight_2d[wh] = weight_2d[wh] / max(weight_2d[wh])
    endif
endfor

;----------------------
; view this field's weight    
;----------------------
fig_dir = '/home/kstory/public_html/notebook/spt_lps12/'

tv_spt_map,shift(weight_2d,nbig/2.,nbig/2.), /norms, min=0,max=1,reso=5,xr=[-3000,3000],yr=[-3000,3000],title='weight_2d=1 / ( GS{cl_th} + GS{noise/tf} )^2, '+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
err = tvread(/png, filename=fig_dir+'kweight_ra5h30dec-55_2008_0501', /nodialog)

; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,xr=[-3000,3000],yr=[-3000,3000],title='un-normalized LOG weight_2d=1 / ( GS{cl_th} + GS{noise/tf} )^2, '+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
err = tvread(/png, filename=fig_dir+'kweight_unnorm_ra5h30dec-55_2008_0501', /nodialog)

; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,title='un-normalized LOG weight_2d=1 / ( GS{cl_th} + GS{noise/tf} )^2, '+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
err = tvread(/png, filename=fig_dir+'kweight_unnorm_zoom_ra5h30dec-55_2008_0501', /nodialog)

; inputs plot
window, 5, xsize=700, ysize=500
vec = indgen(1000)
plot, l[vec,vec], cl_th[vec,vec], /ylog, xtitle='ell', ytitle='Cl [uK-rad]^2', title='ra5h30dec-55_2008', xra=[0,4000]
oplot, l[vec,vec], noise[vec,vec], color=!red
oplot, l[vec,vec], ntmp[vec,vec], color=!orange
legend_str = ['Cl_th', 'noise_psd', 'noise_TF']
legend, legend_str, linestyle=0, colors=[!white, !red, !orange], position=[2500, 0.9e2]
err = tvread(/png, filename=fig_dir+'kweight_inputs_0501', /nodialog)

stop
END

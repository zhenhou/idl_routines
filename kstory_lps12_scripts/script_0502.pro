;;;
; NAME: script_0502
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) process_fields
; 2) check_kweight
; 3) check_sims
;
; MODIFICATION HISTORY:
;  05/02/2012: (KTS) Created
;;;

PRO make_maps_lps12, idx
f = lps12_fieldstruct()
field = f[idx].name

; make the maps
map_script = '/data/kstory/projects/lps12/maps/20120420/make_maps_'+field+'_150.txt'
log_file   = '/data/kstory/projects/lps12/maps/20120420/run0_'+field+'.out'
spawn, '/home/kstory/sptcdevel/mapping/mpibatch.x '+map_script+' > & '+log_file

END

; Process a field
PRO process_field, idx
print, 'Make Maps, idx = ', idx
make_maps_lps12, idx

print, 'coadd maps, idx = ', idx
coadd_maps_0420, idx

print, 'make mask, idx = ', idx
make_mask_lps12, field_idx=idx

print, 'make noise psd, idx = ', idx
make_noise_psd_lps12, idx

;make_coupling_kernel, idx
;try_end2end, idx ??
;make_twod_tfs, idx ??
;make_twod_kweights_lps12, idx ??
END

; List to process
;[2, 9, 12, 13, 14, 15, 17, 18, 19]

;;; Process all new fields
PRO process_all
idx_list = [2, 9, 12, 13, 14, 15, 17, 18, 19]
nlist = n_elements(idx_list)
for jj=0, nlist-1 do begin
    idx = idx_list[jj]

    print, 'process_field, ', idx
    process_field, idx

endfor
END

PRO maps_1314
make_maps_lps12, 13
make_maps_lps12,14
END


PRO process_9121314
idx_list = [9, 12, 13, 14]
nlist = n_elements(idx_list)
for jj=0, nlist-1 do begin
    idx = idx_list[jj]

    print, 'process_field, ', idx
    process_field, idx

endfor
END

;;; Process all new fields
PRO process_15171819
idx_list = [15, 17, 18, 19]
nlist = n_elements(idx_list)
for jj=0, nlist-1 do begin
    idx = idx_list[jj]

    print, 'process_field, ', idx
    process_field, idx

endfor
END

; coupling kernels
PRO make_coupling_kernels
idx_list = 0
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Inputs to kweight
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO check_kweight

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
err = tvread(/png, filename=fig_dir+'kweight_ra5h30dec-55_2008_0502', /nodialog)

; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,xr=[-3000,3000],yr=[-3000,3000],title='un-normalized LOG weight_2d'+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
err = tvread(/png, filename=fig_dir+'kweight_unnorm_ra5h30dec-55_2008_0502', /nodialog)

; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,title='un-normalized LOG weight_2d'+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
err = tvread(/png, filename=fig_dir+'kweight_unnorm_zoom_ra5h30dec-55_2008_0502', /nodialog)

; inputs plot
window, 5, xsize=700, ysize=500
vec = indgen(1000)
plot, l[vec,vec], cl_th[vec,vec], /ylog, xtitle='ell', ytitle='Cl [uK-rad]^2', title='ra5h30dec-55_2008', xra=[0,4000]
oplot, l[vec,vec], noise[vec,vec], color=!red
oplot, l[vec,vec], ntmp[vec,vec], color=!orange
legend_str = ['Cl_th', 'noise_psd', 'noise_TF']
legend, legend_str, linestyle=0, colors=[!white, !red, !orange], position=[2500, 0.9e2]
err = tvread(/png, filename=fig_dir+'kweight_inputs_0502', /nodialog)

stop
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Check sims
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;...................................................................
; Check sims
PRO check_sims

; setup
reso = 1.0
reso_rad = reso/60.*!dtor
field_name = 'ra5h30dec-50_2008'
fig_dir = '/home/kstory/public_html/notebook/spt_lps12/'

;;; get RK's data
restore, '/home/kstory/lps12/twod_tfs/tf_rk/create_twod_tfs.sav'
nbig_rk=2160L
l_rk = make_fft_grid(reso/60.*!dtor/2./!pi,nbig_rk,nbig_rk,fx=lx,fy=ly)
s_rk = s0
power_rk = s_rk.bigpower*1d12/s_rk.bigmask_factor*(nbig_rk*nbig_rk)*(reso_rad^2.)

;tv_spt_map,shift(power_rk,nbig_rk/2.,nbig_rk/2.),reso=(l_rk[1]-l_rk[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='RK, '+field_name, xra=[-4000, 4000], yra=[-4000, 4000], min=1d-4, max=1.2e-2, winnum=3
;err = tvread(/png, filename=fig_dir+'sim_power_rk_0502', /nodialog)

; get KS's data
restore, 'make_twod_tfs.sav'
nbig_ks = 4320L
l_ks = make_fft_grid(reso/60.*!dtor/2./!pi,nbig_ks,nbig_ks,fx=lx,fy=ly)
s_ks = s0
power_ks = s_ks.power*1d12/s_ks.mask_factor*(nbig_ks*nbig_ks)*(reso_rad^2.)


;tv_spt_map,shift(power_ks,nbig_ks/2.,nbig_ks/2.),reso=(l_ks[1]-l_ks[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='KS, '+field_name, xra=[-4000, 4000], yra=[-4000, 4000], min=1d-4, max=1.2e-2, winnum=5
;err = tvread(/png, filename=fig_dir+'sim_power_ks_0502', /nodialog)

;----------------------------
; cut plots: coadded and Single obs
;----------------------------

;;;;;;;;; RK
field = 'ra5h30dec-55'
simdir = '/data/rkeisler/low_ell_sims/output/'
spawn,'ls '+simdir+'coadd_a_'+field+'*.sav',list
nlist=n_elements(list)
restore,'/data/rkeisler/ps09/mask_'+field+'_20101009_055240.sav' ; get mask

nx=960 & ny=960
ca = kscolor_array()
window, 8, xsize=1100,ysize=500
!p.multi=[0,2,1]

plot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), power_rk[0:nbig_rk-1, 0], /ylog, xrange=[0,4000], title='RK, average sim power', xtitle='ell', ytitle='Power, cut in x-dir'

for jj=0, 4 do begin
    restore, list[jj]           ; pick random one
    bigcoadd = fltarr(nbig_rk,nbig_rk)
    bigcoadd[0:nx-1, 0:ny-1] = coadd
    bigmask = fltarr(nbig_rk,nbig_rk)
    bigmask[0:n_elements(mask[*,0])-1, 0:n_elements(mask[0,*])-1] = mask
    bigmask_factor = mean(bigmask^2.)

    rk_pow_1 = (abs(fft(bigcoadd*bigmask))^2.)
    power_rk_1 = rk_pow_1*1d12/bigmask_factor*(nbig_rk*nbig_rk)*(reso_rad^2.)

    if jj eq 0 then plot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), power_rk_1[0:nbig_rk-1, 0], /ylog, yrange=[1e-8,1e0], xrange=[0,4000], title='RK, single sim power', xtitle='ell', ytitle='Power, cut in x-dir'
    oplot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), power_rk_1[0:nbig_rk-1, 0], color=ca[jj]
endfor
!p.multi=0
err = tvread(/png, filename=fig_dir+'sim_power_cut_rk_0502', /nodialog)


;;;;;;;;; KS
f = lps12_fieldstruct()
idx = 0
field_name = f[idx].name
info = get_lps12_fieldinfo(idx) 
nx = info.npix[0] & ny = info.npix[1]
sim_dir  = '/home/kstory/lps12/sims/'
mask_dir = '/home/kstory/lps12/masks/masks_50mJy/'
mcfiles = sim_dir+'coaddsim_'+field_name+'.dat'

mask_padded = get_lps12_mask(0, /padded)
mask_factor = mean(mask_padded^2.)

maps = read_dat_map(mcfiles, nx, ny, 100)

ca = kscolor_array()
!p.multi=[0,2,1]
window, 6, xsize=1100,ysize=500
plot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks[0:nbig_ks-1, 0], /ylog, xrange=[0,4000], title='KS, average sim power', xtitle='ell', ytitle='Power, cut in x-dir'

for jj=0, 4 do begin
    map_padded = pad_array(maps[*,*,10+jj], nbig_ks)
    power_ks_1 = dblarr(nbig_ks, nbig_ks)
    power_ks_1[*,*] += (abs(fft(map_padded*mask_padded))^2.)
    power_ks_1 *= 1d12/mask_factor*(nbig_ks*nbig_ks)*(reso_rad^2.)

    ; re-bin
    p2 = rebin(power_ks_1, 2160, 2160)

;     if jj eq 0 then plot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks_1[0:nbig_ks-1, 0], /ylog, xrange=[0,4000], title='KS, single sim power', xtitle='ell', ytitle='Power, cut in x-dir'
;     oplot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks_1[0:nbig_ks-1, 0], color=ca[jj]

    if jj eq 0 then plot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), p2[0:nbig_rk-1, 0], /ylog, yrange=[1e-8, 1e0], xrange=[0,4000], title='KS, single sim power', xtitle='ell', ytitle='Power, cut in x-dir'
    oplot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), p2[0:nbig_rk-1, 0], color=ca[jj]

endfor
!p.multi=0
err = tvread(/png, filename=fig_dir+'sim_power_cut_ks_0502', /nodialog)


stop
END



;-----------------------------------------------------------
;; temporary
PRO sss
reso = 1.0
reso_rad = reso/60.*!dtor
nbig_ks = 4320L
l_ks = make_fft_grid(reso/60.*!dtor/2./!pi,nbig_ks,nbig_ks,fx=lx,fy=ly)

; coadd
restore, 'make_twod_tfs.sav'
nbig_ks = 4320L
l_ks = make_fft_grid(reso/60.*!dtor/2./!pi,nbig_ks,nbig_ks,fx=lx,fy=ly)
s_ks = s0
power_ks = s_ks.power*1d12/s_ks.mask_factor*(nbig_ks*nbig_ks)*(reso_rad^2.)

; re-bin
p2 = rebin(power_ks, 2160, 2160)

;window, 6 & plot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks[0:nbig_ks-1, 0], xrange=[0,4000], title='KS, average sim power', xtitle='ell', ytitle='Power, cut in x-dir', /ylog
window, 10 & plot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), p2[0:nbig_rk-1, 0], xrange=[0,4000], title='KS, average sim power', xtitle='ell', ytitle='Power, cut in x-dir', /ylog
;err = tvread(/png, filename=fig_dir+'sim_power_cut_ks_0502', /nodialog)

; single obs
f = lps12_fieldstruct()
idx = 0
field_name = f[idx].name
info = get_lps12_fieldinfo(idx) 
nx = info.npix[0] & ny = info.npix[1]
sim_dir  = '/home/kstory/lps12/sims/'
mask_dir = '/home/kstory/lps12/masks/masks_50mJy/'
mcfiles = sim_dir+'coaddsim_'+field_name+'.dat'

mask_padded = get_lps12_mask(0, /padded)
mask_factor = mean(mask_padded^2.)

maps = read_dat_map(mcfiles, nx, ny, 100)

ca = kscolor_array()
!p.multi=[0,2,1]
window, 6, xsize=1100,ysize=500
plot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks[0:nbig_ks-1, 0], xrange=[0,4000], title='KS, average sim power', xtitle='ell', ytitle='Power, cut in x-dir'

for jj=0, 4 do begin
    map_padded = pad_array(maps[*,*,10+jj], nbig_ks)
    power_ks_1 = dblarr(nbig_ks, nbig_ks)
    power_ks_1[*,*] += (abs(fft(map_padded*mask_padded))^2.)
    power_ks_1 *= 1d12/mask_factor*(nbig_ks*nbig_ks)*(reso_rad^2.)

    if jj eq 0 then plot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks_1[0:nbig_ks-1, 0], xrange=[0,4000], title='KS, single sim power', xtitle='ell', ytitle='Power, cut in x-dir'
    oplot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks_1[0:nbig_ks-1, 0], color=ca[jj]
endfor
!p.multi=0

stop
END

;;;
; NAME: script_0503
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) process_fields
; 2) check_kweight
; 3) check_sims
;
; MODIFICATION HISTORY:
;  05/03/2012: (KTS) Created
;;;

PRO make_maps_lps12, idx
f = lps12_fieldstruct()
field = f[idx].name

; make the maps
map_script = '/data/kstory/projects/lps12/maps/20120420/make_maps_'+field+'_150.txt'
log_file   = '/data/kstory/projects/lps12/maps/20120420/run2_'+field+'.out'
spawn, '/home/kstory/sptcdevel/mapping/mpibatch.x '+map_script+' > & '+log_file

END

; first set
PRO maps_26
print, 'Make maps, idx = ', 2
make_maps_lps12, 2
print, 'Make maps, idx = ', 6
make_maps_lps12, 6
END

PRO maps_78
print, 'Make maps, idx = ', 7
make_maps_lps12, 7
print, 'Make maps, idx = ', 8
make_maps_lps12, 8
END

; second set
PRO maps_11121314
print, 'Make maps, idx = ', 11
make_maps_lps12, 11
print, 'Make maps, idx = ', 12
make_maps_lps12, 12
print, 'Make maps, idx = ', 13
make_maps_lps12, 13
print, 'Make maps, idx = ', 14
make_maps_lps12, 14
END

PRO maps_1516171819
print, 'Make maps, idx = ', 15
make_maps_lps12, 15
print, 'Make maps, idx = ', 16
make_maps_lps12, 16
print, 'Make maps, idx = ', 17
make_maps_lps12, 17
print, 'Make maps, idx = ', 18
make_maps_lps12, 18
print, 'Make maps, idx = ', 19
make_maps_lps12, 19
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
; end2end
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO end_0
try_end2end, 0, run='03', /use_kweight, /resume
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; tfs
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;plot,ell[*,0],tf[*,0],xr=[0,5000]



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Check sims
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO ss1
check_sims, 0
check_sims, 1
check_sims, 3
check_sims, 4
;check_sims, 5
END

;...................................................................
; Check sims
PRO check_sims, idx
;idx = 3
rk_idx = idx le 1 ? idx : idx-1
print, 'rk_idx = ', rk_idx
f = lps12_fieldstruct()

; setup
reso = 1.0
reso_rad = reso/60.*!dtor
field_name = f[idx].dir_name
fig_dir = '/home/kstory/public_html/notebook/spt_lps12/'

;;; get RK's data
; 100 coadded sims
restore, '/home/kstory/lps12/twod_tfs/tf_rk/create_twod_tfs_100s.sav'
nbig_rk=2160L
l_rk = make_fft_grid(reso/60.*!dtor/2./!pi,nbig_rk,nbig_rk,fx=lx,fy=ly)
lvec_rk = indgen(nbig_rk)*(l_rk[1]-l_rk[0])
ex = execute('s_rk1 = s'+strtrim(string(rk_idx),2))
power_rk1 = s_rk1.bigpower*1d12/s_rk1.bigmask_factor*(nbig_rk*nbig_rk)*(reso_rad^2.)

; 300 coadded sims
restore, '/home/kstory/lps12/twod_tfs/tf_rk/create_twod_tfs.sav'
ex = execute('s_rk3 = s'+strtrim(string(rk_idx),2))
power_rk3 = s_rk3.bigpower*1d12/s_rk3.bigmask_factor*(nbig_rk*nbig_rk)*(reso_rad^2.)

; window, 2, xsize=1100, ysize=500
; !p.multi=[0,2,1]
; plot, lvec_rk, power_rk[0:nbig_rk-1, 0], /ylog, xrange=[0,4000], title='RK, 100 <psim>'+field_name, xtitle='ell', ytitle='Power, cut in x-dir'
; plot, lvec_rk, power_rk[0:nbig_rk-1, 0], xrange=[0,4000], title='RK, 100 <psim>'+field_name, xtitle='ell', ytitle='Power, cut in x-dir'
; !p.multi=0
; err = tvread(/png, filename=fig_dir+'sim_100_power_rk_'+field_name+'_0503', /nodialog)

;tv_spt_map,shift(power_rk,nbig_rk/2.,nbig_rk/2.),reso=(l_rk[1]-l_rk[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='RK, '+field_name, xra=[-4000, 4000], yra=[-4000, 4000], min=1d-4, max=1.2e-2, winnum=3

; get KS's data
restore, 'make_twod_tfs.sav'
nbig_ks = 4320L
l_ks = make_fft_grid(reso/60.*!dtor/2./!pi,nbig_ks,nbig_ks,fx=lx,fy=ly)
ex = execute('s_ks = s'+strtrim(string(idx),2))
power_ks = s_ks.power*1d12/s_ks.mask_factor*(nbig_ks*nbig_ks)*(reso_rad^2.)
p2 = rebin(power_ks, 2160, 2160)

; PLOT

window, 12, xsize=1100, ysize=500
!p.multi=[0,2,1]
plot, lvec_rk, p2[0:nbig_rk-1, 0], /ylog, xrange=[0,4000], title='lps12, <psim>'+field_name, xtitle='ell', ytitle='Power, cut in x-dir'
oplot, lvec_rk, power_rk3[0:nbig_rk-1, 0], color=!green
oplot, lvec_rk, power_rk1[0:nbig_rk-1, 0], color=!orange
oplot, lvec_rk, p2[0:nbig_rk-1, 0], color=!red

plot, lvec_rk, p2[0:nbig_rk-1, 0], xrange=[0,4000], title='lps12, <psim>'+field_name, xtitle='ell', ytitle='Power, cut in x-dir'
oplot, lvec_rk, power_rk3[0:nbig_rk-1, 0], color=!green
oplot, lvec_rk, power_rk1[0:nbig_rk-1, 0], color=!orange
oplot, lvec_rk, p2[0:nbig_rk-1, 0], color=!red

legend_str = ['lps12', 'k11-100', 'k11-300']
legend, legend_str, linestyle=0, colors=[!red, !orange, !green], position=[2500, 0.9e-2]
!p.multi=0
err = tvread(/png, filename=fig_dir+'sim_power_lps12Vk11_'+field_name+'_0503', /nodialog)

;tv_spt_map,shift(power_ks,nbig_ks/2.,nbig_ks/2.),reso=(l_ks[1]-l_ks[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='KS, '+field_name, xra=[-4000, 4000], yra=[-4000, 4000], min=1d-4, max=1.2e-2, winnum=5
;err = tvread(/png, filename=fig_dir+'sim_power_ks_0503', /nodialog)
;stop
END

; ;----------------------------
; ; cut plots: coadded and Single obs
; ;----------------------------

; ;;;;;;;;; RK
; field = 'ra5h30dec-55'
; simdir = '/data/rkeisler/low_ell_sims/output/'
; spawn,'ls '+simdir+'coadd_a_'+field+'*.sav',list
; nlist=n_elements(list)
; restore,'/data/rkeisler/ps09/mask_'+field+'_20101009_055240.sav' ; get mask

; nx=960 & ny=960
; ca = kscolor_array()
; window, 8, xsize=1100,ysize=500
; !p.multi=[0,2,1]

; plot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), power_rk[0:nbig_rk-1, 0], /ylog, xrange=[0,4000], title='RK, average sim power', xtitle='ell', ytitle='Power, cut in x-dir'

; for jj=0, 4 do begin
;     restore, list[jj]           ; pick random one
;     bigcoadd = fltarr(nbig_rk,nbig_rk)
;     bigcoadd[0:nx-1, 0:ny-1] = coadd
;     bigmask = fltarr(nbig_rk,nbig_rk)
;     bigmask[0:n_elements(mask[*,0])-1, 0:n_elements(mask[0,*])-1] = mask
;     bigmask_factor = mean(bigmask^2.)

;     rk_pow_1 = (abs(fft(bigcoadd*bigmask))^2.)
;     power_rk_1 = rk_pow_1*1d12/bigmask_factor*(nbig_rk*nbig_rk)*(reso_rad^2.)

;     if jj eq 0 then plot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), power_rk_1[0:nbig_rk-1, 0], /ylog, yrange=[1e-8,1e0], xrange=[0,4000], title='RK, single sim power', xtitle='ell', ytitle='Power, cut in x-dir'
;     oplot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), power_rk_1[0:nbig_rk-1, 0], color=ca[jj]
; endfor
; !p.multi=0
; ;err = tvread(/png, filename=fig_dir+'sim_power_cut_rk_0503', /nodialog)


; ;;;;;;;;; KS
; f = lps12_fieldstruct()
; idx = 0
; field_name = f[idx].name
; info = get_lps12_fieldinfo(idx) 
; nx = info.npix[0] & ny = info.npix[1]
; sim_dir  = '/home/kstory/lps12/sims/'
; mask_dir = '/home/kstory/lps12/masks/masks_50mJy/'
; mcfiles = sim_dir+'coaddsim_'+field_name+'.dat'

; mask_padded = get_lps12_mask(0, /padded)
; mask_factor = mean(mask_padded^2.)

; maps = read_dat_map(mcfiles, nx, ny, 100)

; ca = kscolor_array()
; !p.multi=[0,2,1]
; window, 6, xsize=1100,ysize=500
; plot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks[0:nbig_ks-1, 0], /ylog, xrange=[0,4000], title='KS, average sim power', xtitle='ell', ytitle='Power, cut in x-dir'

; for jj=0, 4 do begin
;     map_padded = pad_array(maps[*,*,10+jj], nbig_ks)
;     power_ks_1 = dblarr(nbig_ks, nbig_ks)
;     power_ks_1[*,*] += (abs(fft(map_padded*mask_padded))^2.)
;     power_ks_1 *= 1d12/mask_factor*(nbig_ks*nbig_ks)*(reso_rad^2.)

;     ; re-bin
;     p2 = rebin(power_ks_1, 2160, 2160)

; ;     if jj eq 0 then plot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks_1[0:nbig_ks-1, 0], /ylog, xrange=[0,4000], title='KS, single sim power', xtitle='ell', ytitle='Power, cut in x-dir'
; ;     oplot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks_1[0:nbig_ks-1, 0], color=ca[jj]

;     if jj eq 0 then plot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), p2[0:nbig_rk-1, 0], /ylog, yrange=[1e-8, 1e0], xrange=[0,4000], title='KS, single sim power', xtitle='ell', ytitle='Power, cut in x-dir'
;     oplot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), p2[0:nbig_rk-1, 0], color=ca[jj]

; endfor
; !p.multi=0
; ;err = tvread(/png, filename=fig_dir+'sim_power_cut_ks_0503', /nodialog)


; stop
; END



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
;err = tvread(/png, filename=fig_dir+'sim_power_cut_ks_0503', /nodialog)

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

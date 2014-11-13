;;;
; NAME: script_0430
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) check_noise
;
; MODIFICATION HISTORY:
;  04/30/2012: (KTS) Created
;;;

;...................................................................
; Check noise psds
PRO sss

field_name='ra5h30dec-55_2008'
l = make_fft_grid(1./60.*!dtor,4320,4320)*2.*!pi

; cl_th
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
cl_th = interpol(cl_uK2, l_vec, l)

; noise
restore, '/data/kstory/projects/lps12/noise_psds/noise_psd_'+field_name+'.sav'
; psd is in units of K-rad
psd_uK = psd * 1d6
noise = ( psd_uK )^2

; get the TF for this field
restore, '/data/kstory/projects/lps12/twod_tfs/tf_'+field_name+'.sav'
tf = TF_W_BEAM

!p.multi=[0,1,2] 
plot, l[0:1000], noise[0:1000]
plot, l[0:1000], cl_th[0:1000], yra=[0,1e-4]
!p.multi=0

stop
END


;...................................................................
; Check noise psds
PRO check_noise1

reso = 1.0 ; arcmin
reso_rad = reso/60.*!dtor

; get new noise
restore, '/home/kstory/lps12/noise_psds/noise_psd_ra5h30dec-55_2008.sav'
npix = 4320
; psd should be in units of K-rad.
psd_uKarcmin = psd * 1d6 *60./!dtor
kpower = ( psd_uKarcmin )^2
;kpower = ( psd )^2. * 1d12
;((noise in uK-arcmin)/60.*!dtor)^2.

ks_l = make_fft_grid(reso/60.*!dtor/2./!pi,npix,npix,fx=ks_lx,fy=ks_ly)

;---------------------------
; plots
;---------------------------
kpower_s = shift(kpower, npix/2., npix/2.)
ks_lvec = indgen(npix/2.)*(ks_l[1]-ks_l[0])

;tv_spt_map, Alog(kpower_s), /norms, winnum=4
tv_spt_map,Alog(kpower_s),winnum=4,/norms,xra=[-4000,4000],yra=[-4000,4000],reso=(ks_l[1]-ks_l[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='KS power, ra5h30dec-55'
;err = tvread(/png, filename='/home/kstory/public_html/notebook/spt_lps12/noise_ks_ra5h30dec-55_0430', /nodialog)

!p.multi=[0,2,1]
window, 5, xsize=1000, ysize=500
plot, ks_lvec, psd_uKarcmin[0:npix/2.], xra=[0,4000], title='KS'
plot, ks_lvec, psd_uKarcmin[0:npix/2.], xra=[2000,4000], title='KS'
!p.multi=0

window, 6
;plot, ks_lvec, kpower[0:npix/2.]*ks_lvec*(ks_lvec+1)/(2*!pi), xra=[0,4000], xtitle='ell', ytitle='Dl'
plot, ks_lvec, psd_uKarcmin[0:npix/2.], xra=[2000,4000], title='KS'

stop
END

;...................................................................
; Check noise psds
PRO check_noise

reso = 1.0 ; arcmin
reso_rad = reso/60.*!dtor

; get new noise
restore, '/home/kstory/lps12/noise_psds/noise_psd_ra5h30dec-55_2008.sav'
ks_nbig = 4320
kpower = ( psd * 1./0.25 * !dtor/60. )^2. * 1d12
;((noise in uK-arcmin)/60.*!dtor)^2.

ks_l = make_fft_grid(reso/60.*!dtor/2./!pi,ks_nbig,ks_nbig,fx=ks_lx,fy=ks_ly)

; get old noise
rk_nbig = 2160
restore,'/data/rkeisler/ps09/mask_ra5h30dec-55_20101009_055240.sav'
bigmask = fltarr(rk_nbig,rk_nbig)
bigmask[0:n_elements(mask[*,0])-1, 0:n_elements(mask[0,*])-1] = mask

restore, '/data/rkeisler/ps09/1.39/coadd_ra5h30dec-55.sav'
bigdmap = fltarr(rk_nbig,rk_nbig)
bigdmap[0:n_elements(dmap[*,0])-1, 0:n_elements(dmap[0,*])-1] = dmap
big_dmap_power = abs(fft(bigdmap*bigmask))^2.

maskfac = mean(bigmask^2.)
rpower = big_dmap_power*1d12/maskfac*(rk_nbig*rk_nbig)*(reso_rad^2.)

rk_l = make_fft_grid(reso/60.*!dtor/2./!pi,rk_nbig,rk_nbig,fx=rk_lx,fy=rk_ly)


;---------------------------
; plots
;---------------------------
kpower_s = shift(kpower, ks_nbig/2., ks_nbig/2.)
ks_lvec = indgen(ks_nbig/2.)*(ks_l[1]-ks_l[0])
rpower_s = shift(big_dmap_power, rk_nbig/2., rk_nbig/2.)
rk_lvec = indgen(rk_nbig/2.)*(rk_l[1]-rk_l[0])

;tv_spt_map, Alog(kpower_s), /norms, winnum=4
tv_spt_map,Alog(kpower_s),winnum=4,/norms,xra=[-4000,4000],yra=[-4000,4000],reso=(ks_l[1]-ks_l[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='KS power, ra5h30dec-55'
;err = tvread(/png, filename='/home/kstory/public_html/notebook/spt_lps12/noise_ks_ra5h30dec-55_0430', /nodialog)

!p.multi=[0,2,1]
window, 5, xsize=1000, ysize=500
plot, ks_lvec, kpower[0:ks_nbig/2.], xra=[0,4000], title='KS'
plot, ks_lvec, kpower[0:ks_nbig/2.], xra=[2000,4000], title='KS'
!p.multi=0

tv_spt_map,Alog(rpower_s),winnum=14,/norms,xra=[-4000,4000],yra=[-4000,4000],reso=(rk_l[1]-rk_l[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='RK power, ra5h30dec-55'
;err = tvread(/png, filename='/home/kstory/public_html/notebook/spt_lps12/noise_rk_ra5h30dec-55_0430', /nodialog)

window, 15
plot, rk_lvec, rpower[0:rk_nbig/2.], xra=[0,4000], title='RK'


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

PRO process_0428
;print, "********* process field 10 ****************"
;process_field, 10
END

PRO maps_8
; make the maps
map_script = '/data/kstory/projects/lps12/maps/20120420/make_maps_ra2h30dec-50_150.txt'
log_file   = '/data/kstory/projects/lps12/maps/20120420/run0_ra2h30dec-50.out'
spawn, '/home/kstory/sptcdevel/mapping/mpibatch.x '+map_script+' > & '+log_file

END


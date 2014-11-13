;;;
; NAME: script_0427
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) run end2end on 10, 11, 16
;
; MODIFICATION HISTORY:
;  04/27/2012: (KTS) Created
;;;

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


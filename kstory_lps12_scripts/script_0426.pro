;;;
; NAME: script_0426
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) run lr jacks
; 2) run end2end on 3,4,5
; 3) end_6101116, run end2end on 6,10,11,16
; 4) make_noise, make noise psds
; 5) check_tf, compare TF's agains rk
; 6) plot_for_zhen, make TF plot for zhen
;
; MODIFICATION HISTORY:
;  04/26/2012: (KTS) Created
;;;

;...................................................................
; Run jackknives on fields 10, 11, 16
PRO jack_lr ; done
lps12_jack, 6, 'lr'
END

PRO end_345 ; done
try_end2end, 3, run='02', /resume
try_end2end, 4, run='02', /resume
try_end2end, 5, run='02', /resume
END

PRO end_6101116 ; done
try_end2end, 6, run='02', /resume
try_end2end, 10, run='02', /resume
try_end2end, 11, run='02', /resume
try_end2end, 16, run='02', /resume
END

;...................................................................
; Make noise psds
PRO make_noise
make_noise_psd_lps12, 1
make_noise_psd_lps12, 3
make_noise_psd_lps12, 4
make_noise_psd_lps12, 5
make_noise_psd_lps12, 6
END

PRO make_noise2
make_noise_psd_lps12, 10
make_noise_psd_lps12, 11
make_noise_psd_lps12, 16
END


;...................................................................
; Run jackknives on fields 10, 11, 16
PRO check_tf

; setup
reso = 1.0
reso_rad = reso/60.*!dtor
field_name = 'ra5m30dec-50_2008'
fig_dir = '/home/kstory/public_html/notebook/spt_lps12/'

;;; get RK's data
restore, '/home/kstory/lps12/twod_tfs/tf_rk/create_twodim_tfs.sav'
;nx=960 & ny=960
nbig_rk=2160L
l_rk = make_fft_grid(reso/60.*!dtor/2./!pi,nbig_rk,nbig_rk,fx=lx,fy=ly)
s_rk = s0
power_rk = s_rk.bigpower*1d12/s_rk.bigmask_factor*(nbig_rk*nbig_rk)*(reso_rad^2.)

tv_spt_map,shift(power_rk,nbig_rk/2.,nbig_rk/2.),reso=(l_rk[1]-l_rk[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='RK, '+field_name, xra=[-4000, 4000], yra=[-4000, 4000], min=1d-4, max=1.2e-2, winnum=3
err = tvread(/png, filename=fig_dir+'sim_power_rk_0426', /nodialog)

window, 4 & plot, indgen(nbig_rk)*(l_rk[1]-l_rk[0]), power_rk[0:nbig_rk-1, 0], xrange=[0,4000], title='RK, average sim power', xtitle='ell', ytitle='Power, cut in x-dir'
err = tvread(/png, filename=fig_dir+'sim_power_cut_rk_0426', /nodialog)

; get KS's data
restore, 'make_twodim_tfs.sav'
nbig_ks = 4320L
l_ks = make_fft_grid(reso/60.*!dtor/2./!pi,nbig_ks,nbig_ks,fx=lx,fy=ly)
s_ks = s0
power_ks = s_ks.power*1d12/s_ks.mask_factor*(nbig_ks*nbig_ks)*(reso_rad^2.)


tv_spt_map,shift(power_ks,nbig_ks/2.,nbig_ks/2.),reso=(l_ks[1]-l_ks[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='KS, '+field_name, xra=[-4000, 4000], yra=[-4000, 4000], min=1d-4, max=1.2e-2, winnum=5
err = tvread(/png, filename=fig_dir+'sim_power_ks_0426', /nodialog)

window, 6 & plot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks[0:nbig_ks-1, 0], xrange=[0,4000], title='KS, average sim power', xtitle='ell', ytitle='Power, cut in x-dir'
err = tvread(/png, filename=fig_dir+'sim_power_cut_ks_0426', /nodialog)

stop
END

;...................................................................
; Run jackknives on fields 10, 11, 16
PRO plot_for_zhen

; setup
reso = 1.0
reso_rad = reso/60.*!dtor
field_name = 'ra5m30dec-50_2008'
fig_dir = '/home/kstory/public_html/notebook/spt_lps12/'


; get KS's data
restore, 'make_twodim_tfs.sav'
nbig_ks = 4320L
l_ks = make_fft_grid(reso/60.*!dtor/2./!pi,nbig_ks,nbig_ks,fx=lx,fy=ly)
s_ks = s0
power_ks = s_ks.power*1d12/s_ks.mask_factor*(nbig_ks*nbig_ks)*(reso_rad^2.)

window, 8 & plot, indgen(nbig_ks)*(l_ks[1]-l_ks[0]), power_ks[0:nbig_ks-1, 0], xrange=[800,3000], /ylog, title='KS, average sim power', xtitle='ell', ytitle='log(Power), cut in x-dir'
err = tvread(/png, filename=fig_dir+'sim_power_cut_log_ks_0426', /nodialog)

END

;;;
; NAME: script_0621
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Make 2d kweight for reso=0.5, npix=4320
;
; MODIFICATION HISTORY:
;  06/21/2012: (KTS) Created
;;;




;;;;;;;;;;;;;;;;;;;;;;
; Make 2d kweight 
;;;;;;;;;;;;;;;;;;;;;;

PRO mk_weight
kk = get_lps12_kweight(4)

kk_tmp = shift(congrid(kk, 2160,2160), 1080,1080)

kk_new = dblarr(4320,4320) + 1.

kk_new[1080:3239, 1080:3239] = kk_tmp
kk_new = shift(kk_new, 2160,2160)
kk_new[where(kk_new eq 5)] = 1.

stop
save, weight_2d, filename='/home/kstory/lps12/end2end/run_07p1/crtmp_weight_2d_07p1.sav'
END






;;;;;;;;;;;;;;;;;;;;;;
; Plot ps error bars
;;;;;;;;;;;;;;;;;;;;;;

PRO plot_ps

; plotting setup
if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
whplot = indgen(45) + 10

;run_07
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_144542_kweight.sav'
dl_07 = dl_all
diag07 = diag_nobeam

;run_05
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120614_192137_kweight.sav' ; Should be final
dl_05  = dl_all
diag05 = diag_nobeam

; run_07; 0809 only
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120619_195218_kweight_0809.sav'
dl_0809     = dl_all
diag_0809 = diag_nobeam
dd_e2e = dblarr(58)
for i=0, 57 do dd_e2e[i] = sqrt(cov_all_nobeam_e2e[i,i])

restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_205039_kweight_0809.sav' ; from comb2_0809
;dl2_0809     = dl_all_init
diag2_0809 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_rk = diag_nobeam

; cl_theory
;readcol,'/home/kstory/lps12/theory_spectrum/Cls_100sims_ave_from_combined_map.txt',l_vec,cl_uK2
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)

;----------------
; Plots
;----------------
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

; ; ps plot
; window, 0
; plot,l[whplot],dl_07[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
;   xtitle=xtitle,ytitle=ytitle,chars=chars, title='Full 2500 deg^2, run_07'
; oplot,l_vec,dl_th,color=!red
; errplot,l[whplot],(dl_07-diag07)[whplot]*1d12,(dl_07+diag07)[whplot]*1d12
; ; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0619')

; ; compare 2500 against K11
; window, 1
; plot, l[whplot], ((err_rk )/diag07)[whplot], xtitle=xtitle, ytitle='diag_k11/diag07', title='Error bars, k11 v.s. run_07'
; oplot, [0,10000], [sqrt(2544./790),sqrt(2544./790)], linestyle=1
; ;err=tvread(/png,/nodialog,filename=fdir+'err_k11v07_0619')

; compare 0809 against K11
window, 2
plot, l[whplot], ((err_rk )/diag_0809)[whplot], xtitle=xtitle, ytitle='diag_k11/diag_0809', title='Error bars, comb_0809'
oplot, [0,10000], [1,1], linestyle=1
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v0809_0619')

;window, 3
;plot, l[whplot], ((err_rk )/diag2_0809)[whplot], xtitle=xtitle, ytitle='diag_k11/comb_0809', title='Error bars, comb_0809'

; window, 3
; plot, l[whplot], (dl_05/dl_07)[whplot], xtitle=xtitle, ytitle='dl_05 / dl_07'

; window, 4
; plot, l[whplot], (diag05/diag07)[whplot], xtitle=xtitle, ytitle='dl_05 / dl_07'


stop
END



;;;;;;;;;;;;;;;;;;;;
; Investigate ra3h30dec-60 errors
;;;;;;;;;;;;;;;;;;;;
PRO err_4
; get RK
restore, '/home/rkeisler/ps09/1.37/end_ra3h30dec-60_1.37_kweight.sav'
rk_diag_meas = dblarr(58)
rk_diag_sv   = dblarr(58)
rk_diag      = dblarr(58)
rk_dl        = spectrum
for i=0, 57 do begin
    rk_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    rk_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    rk_diag[i]      = sqrt(cov[i,0,i,0])
endfor
rk_tf = transfer
rk_ellkern = ellkern

; get KS
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07_kweight.sav'
ks_diag_meas = dblarr(58)
ks_diag_sv   = dblarr(58)
ks_diag      = dblarr(58)
ks_dl        = spectrum
for i=0, 57 do begin
    ks_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks_diag[i]      = sqrt(cov[i,0,i,0])
    ;ks_diag[i]      = sqrt( ks_diag_meas[i]^2. + ks_diag_sv[i]^2.)
endfor
ks_tf = transfer
ks_ellkern = ellkern

; get KS run_07p1 (pad to 2160)
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07p1_kweight.sav'
ks07p1_diag_meas = dblarr(58)
ks07p1_diag_sv   = dblarr(58)
ks07p1_diag      = dblarr(58)
ks07p1_dl        = spectrum
for i=0, 57 do begin
    ks07p1_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks07p1_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks07p1_diag[i]      = sqrt(cov[i,0,i,0])
endfor

; get KS run_07p2 (use end2end_lowell)
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07p2_kweight.sav'
ks07p2_diag_meas = dblarr(58)
ks07p2_diag_sv   = dblarr(58)
ks07p2_diag      = dblarr(58)
ks07p2_dl        = spectrum
for i=0, 57 do begin
    ks07p2_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks07p2_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks07p2_diag[i]      = sqrt(cov[i,0,i,0])
endfor

; Try this:
; ks_diag_meas *= 0.825
; ks_diag_sv *= sqrt(0.825)
; ;ks_diag *= sqrt(0.825)
; ks_diag = sqrt(ks_diag_meas^2. + ks_diag_sv^2.)

; run_07; 0809 only
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120619_195218_kweight_0809.sav'
diag_0809 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
k11_diag = diag_nobeam

;---------------------------
; PLOTS
;---------------------------
whplot = indgen(45) + 10
fdir = '/home/kstory/public_html/notebook/spt_lps12/'
l = banddef-25.

window, 0, xsize=800, ysize=500
plot, l[whplot], (rk_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, ra3h30dec-60', $ ; tot
  ytitle='err_rk / err_ks', xtitle='ell'
oplot, l[whplot], (rk_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (rk_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
oplot, l[whplot], (k11_diag/diag_0809)[whplot], color=!red, linestyle=3       ; 0809
;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0621')

; window, 1, xsize=800, ysize=900
; !p.multi=[0,1,2]
; plot, l[whplot], (ks07p1_diag/ks_diag)[whplot], yr=[.95,1.05],/yst,title='Errors, changing mask padding', $ ; tot
;   ytitle='err_2160 / err_4320', xtitle='ell'
; oplot, l[whplot], (ks07p1_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
; oplot, l[whplot], (ks07p1_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
; oplot, l[whplot], (k11_diag/diag_0809)[whplot], color=!red, linestyle=3       ; 0809

; window, 1
; plot, l[whplot], (ks07p1_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, ra3h30dec-60', $ ; tot
;   ytitle='err_2160 / err_4320', xtitle='ell'
; oplot, l[whplot], (ks07p1_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
; oplot, l[whplot], (ks07p1_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
; oplot, l[whplot], (k11_diag/diag_0809)[whplot], color=!red, linestyle=3       ; 0809
; legend, ['Tot', 'SV', 'NV', 'K11 v 0809'], color=[!white, !green, !skyblue, !red], linestyle=[0,0,0,3], pos=[600, 1.2]
; ;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0620')
; !p.multi=0

window, 2
plot, l[whplot], (ks07p2_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, 7p2, ra3h30dec-60', $ ; tot
  ytitle='err_2160 / err_4320', xtitle='ell'
oplot, l[whplot], (ks07p2_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (ks07p2_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
oplot, l[whplot], (k11_diag/diag_0809)[whplot], color=!red, linestyle=3       ; 0809
legend, ['Tot', 'SV', 'NV', 'K11 v 0809'], color=[!white, !green, !skyblue, !red], linestyle=[0,0,0,3], pos=[600, 1.2]
;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0620')
!p.multi=0


; window, 2
; plot, ks_ellkern, ks_tf, xr=[0,4500], title='TFs, ra3h30dec-60'
; oplot, rk_ellkern, rk_tf, color=!red
; legend, ['KS', 'RK'], color=[!white, !red], linestyle=[0,0], pos=[3000,0.4]
; ;err=tvread(/png,/nodialog,filename=fdir+'tf_field4_0620')

; window, 3
; rk_tf2 = interpol(rk_tf, 1860)
; plot, ks_ellkern, rk_tf2/ks_tf, xr=[650,3000], /xst, yr=[0.9,1.2],/yst, $
;   title='TFs, ra3h30dec-60', xtitle=ell, ytitle='TF_RK / tf_KS'
;plot, ks_ellkern, ks_tf, xr=[0,4500], title='TFs, ra3h30dec-60'
;oplot, ks_ellkern, rk_tf2, color=!red
;legend, ['KS', 'RK'], color=[!white, !red], linestyle=[0,0], pos=[3000,0.4]
;err=tvread(/png,/nodialog,filename=fdir+'tf_frac_field4_0620')

stop
END





PRO err_4p2
; get KS
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07_kweight.sav'
ks_diag_meas = dblarr(58)
ks_diag_sv   = dblarr(58)
ks_diag      = dblarr(58)
ks_dl        = spectrum
ks_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks_diag[i]      = sqrt(cov[i,0,i,0])
    ;ks_diag[i]      = sqrt( ks_diag_meas[i]^2. + ks_diag_sv[i]^2.)
endfor
ks_tf = transfer
ks_ellkern = ellkern

; get KS run_07p2 (use end2end_lowell)
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07p2_kweight.sav'
ks07p2_diag_meas = dblarr(58)
ks07p2_diag_sv   = dblarr(58)
ks07p2_diag      = dblarr(58)
ks07p2_dl        = spectrum
ks07p2_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks07p2_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks07p2_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks07p2_diag[i]      = sqrt(cov[i,0,i,0])
endfor

; get KS run_07p3 (use end2end original)
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07p3_kweight.sav'
ks07p3_diag_meas = dblarr(58)
ks07p3_diag_sv   = dblarr(58)
ks07p3_diag      = dblarr(58)
ks07p3_dl        = spectrum * 0.825^2.
ks07p3_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks07p3_diag_meas[i] = sqrt(meas_cov[i,0,i,0]) * (0.825^4.)
    ks07p3_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ;ks07p3_diag[i]      = sqrt(cov[i,0,i,0])
    ks07p3_diag[i]      = sqrt(meas_cov[i,0,i,0]*(0.825^4.) + sample_cov[i,0,i,0])
endfor

; get KS run_07p4 (use end2end_lps12, re-make kernel)
restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_07p4_kweight.sav'
ks07p4_diag_meas = dblarr(58)
ks07p4_diag_sv   = dblarr(58)
ks07p4_diag      = dblarr(58)
ks07p4_dl        = spectrum
ks07p4_cov_mc_raw = cov_mc_raw
for i=0, 57 do begin
    ks07p4_diag_meas[i] = sqrt(meas_cov[i,0,i,0])
    ks07p4_diag_sv[i]   = sqrt(sample_cov[i,0,i,0])
    ks07p4_diag[i]      = sqrt(cov[i,0,i,0])
endfor

;---------------------------
; PLOTS
;---------------------------
whplot = indgen(45) + 10
fdir = '/home/kstory/public_html/notebook/spt_lps12/'
l = banddef-25.

window, 2
plot, l[whplot], (ks07p2_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, 7p2, ra3h30dec-60', $ ; tot
  ytitle='err_2160 / err_4320', xtitle='ell'
oplot, l[whplot], (ks07p2_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (ks07p2_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
;oplot, l[whplot], (k11_diag/diag_0809)[whplot], color=!red, linestyle=3       ; 0809
legend, ['Tot', 'SV', 'NV', 'K11 v 0809'], color=[!white, !green, !skyblue, !red], linestyle=[0,0,0,3], pos=[600, 1.2]
;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0620')

window, 3
plot, l[whplot], (ks07p3_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, 7p3, ra3h30dec-60', $ ; tot
  ytitle='err_2160 / err_4320', xtitle='ell'
oplot, l[whplot], (ks07p3_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (ks07p3_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0620')

window, 4
plot, l[whplot], (ks07p4_diag/ks_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, 7p4, ra3h30dec-60', $ ; tot
  ytitle='err_2160 / err_4320', xtitle='ell'
oplot, l[whplot], (ks07p4_diag_meas/ks_diag_meas)[whplot], color=!skyblue         ; meas
oplot, l[whplot], (ks07p4_diag_sv/ks_diag_sv)[whplot], color=!green               ; sv
;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0620')

; window, 4
; plot, l[whplot], (ks07p3_diag/ks07p2_diag)[whplot], yr=[.7,1.3],/yst,title='Errors, 7p3, ra3h30dec-60', $ ; tot
;   ytitle='err_2160 / err_4320', xtitle='ell'
; oplot, l[whplot], (ks07p3_diag_meas/ks07p2_diag_meas)[whplot], color=!skyblue         ; meas
; oplot, l[whplot], (ks07p3_diag_sv/ks07p2_diag_sv)[whplot], color=!green               ; sv
; ;err=tvread(/png,/nodialog,filename=fdir+'err_field4_0620')

stop
END

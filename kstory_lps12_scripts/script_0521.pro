;;;
; NAME: script_0521
; PURPOSE:
;   General script for today
;
; NOTES:
;
;
; MODIFICATION HISTORY:
;  05/21/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; compare error bars for individual fields
PRO comp_err;, idx
bdir = '/home/kstory/lps12/end2end/'

; get ks_new
restore, bdir+'end_ra4h10dec-50_03.sav'
dtot_ks_new = dblarr(nl)   ; raw cov diags
ds_ks_new   = dblarr(nl)   ; sample variance
dm_ks_new   = dblarr(nl)   ; noise variance

dls_ks_new = spectrum & cov_ks_new = reform(cov)
ellkern_ks_new = ellkern & nb = n_elements(ellkern_ks_new)
tf_ks_new = transfer
beam_ks_new = beam_interp
for j=0,nl-1 do begin 
    dtot_ks_new[j] = sqrt(cov_ks_new[j,j])
    ds_ks_new[j]  = sqrt(sample_cov[j,0,j,0])
    dm_ks_new[j]  = sqrt(meas_cov[j,0,j,0])
endfor

; get ks_badBeam
restore, bdir+'run_03_badBeam/end_ra4h10dec-50_03.sav'
dtot_ks_badBeam = dblarr(nl)   ; raw cov diags
ds_ks_badBeam   = dblarr(nl)   ; sample variance
dm_ks_badBeam   = dblarr(nl)   ; noise variance

dls_ks_badBeam = spectrum & cov_ks_badBeam = reform(cov)
ellkern_ks_badBeam = ellkern & nb = n_elements(ellkern_ks_badBeam)
tf_ks_badBeam = transfer
beam_ks_badBeam = beam_interp
for j=0,nl-1 do begin 
    dtot_ks_badBeam[j] = sqrt(cov_ks_badBeam[j,j])
    ds_ks_badBeam[j]  = sqrt(sample_cov[j,0,j,0])
    dm_ks_badBeam[j]  = sqrt(meas_cov[j,0,j,0])
endfor

xtot = (dtot_ks_new - dtot_ks_badBeam)/dtot_ks_new
xs = (ds_ks_new - ds_ks_badBeam)/ds_ks_new
xm = (dm_ks_new - dm_ks_badBeam)/dm_ks_new

x_tf = (tf_ks_new - tf_ks_badBeam)/tf_ks_new
x_beam = (beam_ks_new - beam_ks_badBeam)/beam_ks_new

;;;;;;;;;;;;;;;;;;;;;
; plots
;;;;;;;;;;;;;;;;;;;;;
figdir = '/home/kstory/public_html/notebook/spt_lps12/'
wh  = where(banddef gt 600)
whe = where(ellkern ge 650 and ellkern lt 3000)

window, 0
plot, banddef[wh], xtot[wh], title='fractional change in errors, ra4h10dec-50', xtitle='ell', ytitle='(dDl_new - dDl_old) / dDl_new'
;err = tvread(/png,/nodialog,filename=figdir+'frac_err_ra4h10dec-50_0521')

window, 1
plot, banddef[wh], dtot_ks_new[wh], /ylog, yr=[1e-12, 1e-9], title='errors, ra4h10dec-50', xtitle='ell', ytitle='dDl'
oplot, banddef[wh], dtot_ks_new[wh], linestyle=3
oplot, banddef[wh], ds_ks_new[wh], color=!red
oplot, banddef[wh], dm_ks_new[wh], color=!yellow
legend, ['dDl_tot', 'dDl_sv', 'dDl_meas'], colors=[!white, !red, !yellow], position=[2300, 1e-9], linestyle=[3,0,0]
;err = tvread(/png,/nodialog,filename=figdir+'err_ra4h10dec-50_0521')

window, 2
plot, ellkern_ks_new[whe], x_tf[whe], xr=[0,3000], yr=[-0.04, 0.04], title='TF, errors, ra4h10dec-50', ytitle='(tf_new - tf_old)/tf_new', xtitle='ell'
oplot, ellkern_ks_new[whe], x_tf[whe], color=!purple
oplot, banddef[wh], xtot[wh], linestyle=3
legend, ['TF comparison', 'dDl comparison'], colors=[!purple, !white], position=[1000, 0.1], linestyle=[0,3]
;err = tvread(/png,/nodialog,filename=figdir+'tf_err_ra4h10dec-50_0521')

stop
END


; look at differences in beams
PRO plot_beams
bf = get_lps12_beams(2009, l09, b09)
bf = get_lps12_beams(2010, l10, b10)
wh= where(l09 gt 600 and l09 lt 3000)

figdir = '/home/kstory/public_html/notebook/spt_lps12/'

window, 1
plot, l09, b09, xr=[0,3200], yr=[0.75, 1.1], xtitle='ell', ytitle='Bl', title='beams, 2009 v 2010-prelim'
oplot, l10, b10, color=!red, linestyle=3
oplot, [650,650], [0, 10], linestyle=2
oplot, [3000,3000], [0, 10], linestyle=2
legend, ['2009', '2010-prelim'], colors=[!white, !red], position=[1500, 1.05], linestyle=[0,3]
err = tvread(/png,/nodialog,filename=figdir+'beams_09_10_0521')

window, 2
plot, l09[wh], (b09[wh]-b10[wh])/b09[wh], xtitle='ell', ytitle='(b09-b10)/b09', title='Beam difference, 2009 2010'
err = tvread(/png,/nodialog,filename=figdir+'frac_beams_09_10_0521')

stop
END




;----------------------------------
;;; combined error bars
PRO comp_err_comb
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120520_001439_kweight.sav' & err_ks = diag_nobeam
;restore, '/home/kstory/lps12/end2end/run_03/obsolete/combined_spectrum_20120521_203118_kweight.sav'
;restore, '/home/kstory/lps12/end2end/run_03/obsolete/combined_spectrum_20120521_212127_kweight.sav' ; looks best
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120521_221658_kweight.sav'
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120521_223459_kweight.sav'
restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120521_224750_kweight.sav'
err_nobeam_ks = diag_nobeam & err_ks = diag

restore, '/home/kstory/lps12/end2end/run_03_badBeam/run_03/combined_spectrum_20120511_200438_kweight.sav'
err_nobeam_ks_old = diag_nobeam & err_ks_old = diag

restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_nobeam_rk = diag_nobeam & err_rk = diag

l_bdf = l
x = (err_nobeam_rk - err_nobeam_ks)/err_nobeam_rk
v = 6+indgen(58-6)
res = linfit( l_bdf[v], x[v] )

; PLOTS
window, 1
plot, l_bdf[v], sqrt(err_nobeam_rk[v]/err_nobeam_ks[v]), title='change in error error', xtitle='ell', ytitle='err_k11 / err_lps12'

stop
END




;;; calculate the sample variance error bars for lps12 vs k11 based on
;;; the TF
PRO lx_cut

; Get the actual difference in error bars between lps12 and k11 based on TF
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120516_190917_kweight_0809.sav' & err_ks = diag_nobeam
restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120517_223348_0809.sav' & err_ks = diag_nobeam
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav' & err_rk = diag_nobeam
l_bdf = l
x = (err_rk - err_ks)/err_rk
v = 6+indgen(58-6)
;v = 2+indgen(58-2)
res = linfit( l_bdf[v], x[v] )


;-------------------------
; Setup
;-------------------------
reso = 1.0
nbig_ks = 4320L
nbig_rk = 2160L
f = lps12_fieldstruct()
tf_dir = '/home/kstory/lps12/twod_tfs/'
kweight_dir = '/home/kstory/lps12/twod_kweights/'

;-------------------------
; get BANDDEF
;-------------------------
delta_l = 50.
min_l = 250.
max_l = 3150.
nl = floor((max_l-min_l)/delta_l)
lhi = dindgen(nl)*delta_l + min_l
banddef = lhi
; this can possibly be much lower:
maxell = max(banddef)*1.5

; get the l-grid
l_ks = make_fft_grid(reso/60.*!dtor,nbig_ks,nbig_ks)*2.*!pi
l_rk = make_fft_grid(reso/60.*!dtor,nbig_rk,nbig_rk)*2.*!pi

;; PLOT
window, 1
ksc = kscolor_array()
plot, l_bdf[v], x[v], title='sample variance error', xtitle='ell', ytitle='(err_k11 - err_lps12)/err_k11'

; loop over all 08 09 fields
idx_list = [0,1,3,4,5]
;idx_list = [0]
for ii=0, 4 do begin
    idx = idx_list[ii]

    ; get TF_ks
    restore, tf_dir+'tf_'+f[idx].name+'.sav'
    tf_ks = TF
    ;kweight_ks = get_lps12_kweight(idx)
    kweight_ks = intarr(nbig_ks, nbig_ks) + 1

    ; get tf_rk
    restore, '/home/rkeisler/ps09/twod_tfs/tf_pad_'+f[idx].dir_name+'.sav'
    tf_rk = BIG_TF
    restore, '/home/rkeisler/ps09/kweight/weight_2d_'+f[idx].dir_name+'_20100117.sav'
    kweight_rk = WEIGHT_2D


    x_ks  = fltarr(nl)
    nx_ks = fltarr(nl)
    
    x_rk  = fltarr(nl)
    nx_rk = fltarr(nl)

    ; get indicies of annuli, corresponding to banddef
    for ib=0, nl-2 do begin
        if ib mod 10 eq 2 then print, 'ib = ', ib
        
        ; get the proxy for sample variance
        wh_ks = where( (l_ks gt banddef[ib]) and (l_ks lt banddef[ib+1]) , nwh_ks)
        x_ks[ib]  = 1./sqrt( total( TF_KS[wh_ks] * KWEIGHT_KS[wh_ks] ) )
        nx_ks[ib] = nwh_ks
        
        ; get the proxy for sample variance
        wh_rk = where( (l_rk gt banddef[ib]) and (l_rk lt banddef[ib+1]) , nwh_rk)
        x_rk[ib]  = 1./sqrt( total( TF_RK[wh_rk] * KWEIGHT_RK[wh_rk]  ) )
        nx_rk[ib] = nwh_rk
        
    endfor

    ; add to plot
    oplot, banddef[v], (x_rk[v] - x_ks[v]*2)/x_rk[v], color=ksc[ii]
        
endfor

;plot, banddef[v], (x_ks[v]*2 - x_rk[v])/x_rk[v], title='sample variance error, '+f[idx].name, xtitle='ell', ytitle='(err_lps12 - err_k11)/err_kll'
;plot, banddef[v], x_rk[v], 
;oplot, banddef[v], x_ks[v]*2, color=!red

bdir = '/home/kstory/public_html/notebook/spt_lps12/'
;err=tvread(/png,/nodialog,filename=bdir+'sv_0517')
stop

END



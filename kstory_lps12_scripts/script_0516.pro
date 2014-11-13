;;;
; NAME: script_0516
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) comp_error_0809, compare end2end error bars between lps12 and
;    k11, 2008 and 2009 fields
; 2) lx_cut, look at the effect of the lx cut on sample variance error bars
; 3) run azrms_70 and azrms_90 jackknives for both data and sims
;
;
; MODIFICATION HISTORY:
;  05/16/2012: (KTS) Created
;;;

;;; compare error bars for 08 09 fields against k11
PRO comp_error_0809
restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120516_190917_kweight_0809.sav' & err_ks = diag_nobeam
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav' & err_rk = diag_nobeam
x = (err_rk - err_ks)/err_rk
v = 6+indgen(52)
res = linfit( l[v], x[v] )

plot, l[v], x[v]
oplot, l[v], res[0] + l[v]*res[1], linestyle=3
stop

END

;;; calculate the sample variance error bars for lps12 vs k11 based on
;;; the TF
PRO lx_cut

; Get the actual difference in error bars between lps12 and k11 based on TF
restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120516_190917_kweight_0809.sav' & err_ks = diag_nobeam
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav' & err_rk = diag_nobeam
l_bdf = l
x = (err_rk - err_ks)/err_rk
;v = 6+indgen(58-6)
v = 2+indgen(58-2)
res = linfit( l_bdf[v], x[v] )


;-------------------------
; Setup
;-------------------------
reso = 1.0
nbig_ks = 4320L
nbig_rk = 2160L
f = lps12_fieldstruct()
tf_dir = '/home/kstory/lps12/twod_tfs/'

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
;idx_list = [0,1,3,4,5]
idx_list = [0]
for ii=0, 4 do begin
    idx = idx_list[ii]

    ; get TF_ks
    restore, tf_dir+'tf_'+f[idx].name+'.sav'
    tf_ks = TF

    ; get tf_rk
    restore, '/home/rkeisler/ps09/twod_tfs/tf_pad_'+f[idx].dir_name+'.sav'
    tf_rk = BIG_TF

    x_ks  = fltarr(nl)
    nx_ks = fltarr(nl)
    
    x_rk  = fltarr(nl)
    nx_rk = fltarr(nl)

    ; get indicies of annuli, corresponding to banddef
    for ib=0, nl-2 do begin
        if ib mod 10 eq 2 then print, 'ib = ', ib
        
        ; get the proxy for sample variance
        wh_ks = where( (l_ks gt banddef[ib]) and (l_ks lt banddef[ib+1]) , nwh_ks)
        x_ks[ib]  = 1./sqrt( total( TF_KS[wh_ks] ) )
        nx_ks[ib] = nwh_ks
        
        ; get the proxy for sample variance
        wh_rk = where( (l_rk gt banddef[ib]) and (l_rk lt banddef[ib+1]) , nwh_rk)
        x_rk[ib]  = 1./sqrt( total( TF_RK[wh_rk] ) )
        nx_rk[ib] = nwh_rk
        
        ; add to plot
        oplot, banddef[v], (x_rk[v] - x_ks[v]*2)/x_rk[v], color=ksc[ii]
        
    endfor
endfor

;plot, banddef[v], (x_ks[v]*2 - x_rk[v])/x_rk[v], title='sample variance error, '+f[idx].name, xtitle='ell', ytitle='(err_lps12 - err_k11)/err_kll'
;plot, banddef[v], x_rk[v], 
;oplot, banddef[v], x_ks[v]*2, color=!red

bdir = '/home/kstory/public_html/notebook/spt_lps12/'
;err=tvread(/png,/nodialog,filename=bdir+'sv_'+f[idx].name+'_0516')
stop

END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Processing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Make azrms_70 files
PRO goodfiles_azrms_70
for ii=0, 19 do make_azrms_jack_list_keeppct, ii, 70
END

; Make azrms_90 files
PRO goodfiles_azrms_90
for ii=0, 19 do make_azrms_jack_list_keeppct, ii, 90
END


; run jackknives
PRO azrms_70_jack
for idx=0, 19 do begin
    print, 'azrms_70 jack, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_70'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO azrms_90_jack
for idx=0, 19 do begin
    print, 'azrms_90 jack, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_90'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


;;; sim jackknives
PRO azrms_70_jack_sim
for idx=0, 19 do begin
    print, 'coadd jack sims azrms_70, field ', idx
    t0 = systime(0, /seconds)
    coadd_lps12_jack_sims, idx, 'azrms_70'
    print, 'That coadd took ', (systime(0, /seconds) - t0)/60., ' minutes.'

    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_70', /sim
    print, 'That jack sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO azrms_90_jack_sim
for idx=0, 19 do begin
    print, 'coadd jack sims azrms_90, field ', idx
    t0 = systime(0, /seconds)
    coadd_lps12_jack_sims, idx, 'azrms_90'
    print, 'That coadd took ', (systime(0, /seconds) - t0)/60., ' minutes.'

    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_90', /sim
    print, 'That jack sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

;;;
; NAME: script_0530
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) compare error bars
;
; MODIFICATION HISTORY:
;  05/30/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; compare errors between run_03 and run_04
PRO comp_er, idx

; setup
f = lps12_fieldstruct()
end_dir = '/home/kstory/lps12/end2end/'

field = f[idx].name
f03 = end_dir+'run_03/end_'+field+'_03_kweight.sav'
f04 = end_dir+'end_'+field+'_04_kweight.sav'

nbins_off_diag = 5
reso=1.0

; get stuff from f03
print,idx, ',  '+f03
restore,f03
cov_03 = reform(cov)
l = banddef-25.
nl = n_elements(l)
area_03 = total(mask_padded)*(reso/60.)^2.

cov_cond_03 = condition_cov_matrix(cov_03, nbins_off_diag, $
                                corr=this_corr, uncondcorr=this_uncorr, $
                                banddef=banddef, $
                                ellmin=ellmin, ellmax=ellmax, $
                                noaverage=noaverage,$
                                knowncorr=knowncorr)
these_diags2_03 = dblarr(nl)
for j=0,nl-1 do these_diags2_03[j] = sqrt(cov_cond_03[j,j])
diags2_03 = these_diags2_03

; get stuff from f04
print,idx, ',  '+f04
restore,f04
cov_04 = reform(cov)
l = banddef-25.
nl = n_elements(l)
area_04 = total(mask_padded)*(reso/60.)^2.

cov_cond_04 = condition_cov_matrix(cov_04, nbins_off_diag, $
                                corr=this_corr, uncondcorr=this_uncorr, $
                                banddef=banddef, $
                                ellmin=ellmin, ellmax=ellmax, $
                                noaverage=noaverage,$
                                knowncorr=knowncorr)
these_diags2_04 = dblarr(nl)
for j=0,nl-1 do these_diags2_04[j] = sqrt(cov_cond_04[j,j])
diags2_04 = these_diags2_04


; plots
whplot=indgen(46)+10
plot, l[whplot], (diags2_03/diags2_04)[whplot], yr=[0.85, 1.]
print, 'expected: ', sqrt(area_04 / area_03)
stop
END






;;; compare error bars for individual fields, from TF and kweight
PRO tf_err, idx
; setup
f = lps12_fieldstruct()
field = f[idx].name
tf_dir = '/home/kstory/lps12/twod_tfs/'
kweight_dir = '/home/kstory/lps12/twod_kweights/'
nbig = 4320L
l_ks = make_fft_grid(1.0/60.*!dtor,nbig,nbig)*2.*!pi
v = 6+indgen(58-6)

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

; get TF_04
restore, tf_dir+'tf_'+f[idx].name+'.sav'
tf_04 = TF
kweight_04 = get_lps12_kweight(idx)

; get TF_03
restore, tf_dir+'run_03/tf_'+f[idx].name+'.sav'
tf_03 = TF
restore, kweight_dir+'run_03/weight_2d_'+f[idx].name+'.sav'
kweight_03 = weight_2d

x_04  = fltarr(nl)
y_04  = fltarr(nl)
nx_04 = fltarr(nl)

x_03  = fltarr(nl)
y_03  = fltarr(nl)
nx_03 = fltarr(nl)

; get indicies of annuli, corresponding to banddef
for ib=0, nl-2 do begin
    if ib mod 10 eq 2 then print, 'ib = ', ib
    
    ; get the proxy for sample variance
    wh = where( (l_ks gt banddef[ib]) and (l_ks lt banddef[ib+1]) , nwh)
    x_04[ib]  = 1./sqrt( total( TF_04[wh] * KWEIGHT_04[wh] ) )
    nx_04[ib] = nwh
    y_04[ib]  = 1./sqrt( total( TF_04[wh] ) )
    
    ; get the proxy for sample variance
    x_03[ib]  = 1./sqrt( total( TF_03[wh] * KWEIGHT_03[wh]  ) )
    nx_03[ib] = nwh
    y_03[ib]  = 1./sqrt( total( TF_03[wh] ) )
    
endfor

window, 0
plot, banddef[v], (x_03/x_04)[v], yr=[0.9, 1.1], title='TF*kweight'

window, 1
plot, banddef[v], (y_03/y_04)[v], yr=[0.9, 1.1], title='TF only'

stop
END


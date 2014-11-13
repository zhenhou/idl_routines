;;;
; NAME: script_0531
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) comp_err, compare error bars from end2end run 03 and 04.
; 2) ra23h30dec50, compare error bars between 2008 and 2010.  2008
; wins.
; 3) make azrms_cut runlists.
;
; MODIFICATION HISTORY:
;  05/31/2012: (KTS) Created
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

nbig = 4320L
;l_ks = make_fft_grid(reso/60.*!dtor,nbig,nbig)*2.*!pi

; get stuff from f03
print,idx, ',  '+f03
restore,f03
cov_03 = reform(cov)
l = banddef-25.
nl = n_elements(l)
area_03 = total(mask_padded)*(reso/60.)^2.
tf03 = transfer

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
tf04 = transfer

cov_cond_04 = condition_cov_matrix(cov_04, nbins_off_diag, $
                                corr=this_corr, uncondcorr=this_uncorr, $
                                banddef=banddef, $
                                ellmin=ellmin, ellmax=ellmax, $
                                noaverage=noaverage,$
                                knowncorr=knowncorr)
these_diags2_04 = dblarr(nl)
for j=0,nl-1 do these_diags2_04[j] = sqrt(cov_cond_04[j,j])
diags2_04 = these_diags2_04


err_03 = diags2_03 * 0.
err_04 = diags2_04 * 0.

for ib=0, nl-2 do begin
    if ib mod 10 eq 2 then print, 'ib = ', ib
    
    ; get the proxy for sample variance
    wh = where( (ellkern gt banddef[ib]) and (ellkern lt banddef[ib+1]) , nwh)
    err_04[ib+1]  += total(1./tf04[wh])
    err_03[ib+1]  += total(1./tf03[wh])

endfor

exp = sqrt( (1/area_03) / (1/area_04) ) - 1.
print, 'expected: ', exp
exp_err = (err_03/err_04) - 1. + exp

; plots
whplot=indgen(46)+10
plot, l[whplot], (diags2_03/diags2_04)[whplot]-1, yr=[-0.15, 0.15], title='Error bars, field '+field, $
  xtitle='ell', ytitle='(err_run03 / err_run04) - 1'
oplot, l[whplot], (exp_err)[whplot], color=!red

fdir='/home/kstory/public_html/notebook/spt_lps12/'
stop
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ra23h30dec-55: 2008 or 2010?
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; compare errors between 2008 and 2010 data for ra23h30dec-55
PRO ra23h30dec55

; setup
f = lps12_fieldstruct()
end_dir = '/home/kstory/lps12/end2end/'

;field_01 = f[1].name
;field_02 = f[2].name
f01 = end_dir+'end_'+f[1].name+'_04_kweight.sav'
f02 = end_dir+'end_'+f[2].name+'_04_kweight.sav'

nbins_off_diag = 5
reso=1.0

; get stuff from f01
print, f01
restore,f01
cov_01 = reform(cov)
l = banddef-25.
nl = n_elements(l)
area_01 = total(mask_padded)*(reso/60.)^2.
tf01 = transfer

cov_cond_01 = condition_cov_matrix(cov_01, nbins_off_diag, $
                                corr=this_corr, uncondcorr=this_uncorr, $
                                banddef=banddef, $
                                ellmin=ellmin, ellmax=ellmax, $
                                noaverage=noaverage,$
                                knowncorr=knowncorr)
these_diags2_01 = dblarr(nl)
for j=0,nl-1 do these_diags2_01[j] = sqrt(cov_cond_01[j,j])
diags2_01 = these_diags2_01

; get stuff from f02
print, f02
restore,f02
cov_02 = reform(cov)
l = banddef-25.
nl = n_elements(l)
area_02 = total(mask_padded)*(reso/60.)^2.
tf02 = transfer

cov_cond_02 = condition_cov_matrix(cov_02, nbins_off_diag, $
                                corr=this_corr, uncondcorr=this_uncorr, $
                                banddef=banddef, $
                                ellmin=ellmin, ellmax=ellmax, $
                                noaverage=noaverage,$
                                knowncorr=knowncorr)
these_diags2_02 = dblarr(nl)
for j=0,nl-1 do these_diags2_02[j] = sqrt(cov_cond_02[j,j])
diags2_02 = these_diags2_02

; plots
whplot=indgen(46)+10
plot, l[whplot], (diags2_01/diags2_02)[whplot]-1, yr=[-0.15, 0.15], title='Error bars, field '+f[1].dir_name, $
  xtitle='ell', ytitle='(err_2008 / err_2010) - 1'

fdir='/home/kstory/public_html/notebook/spt_lps12/'
stop
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make azrms runlists for end2end 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO azrms_runlists, skeep
;skeep = '90'

; setup
f = lps12_fieldstruct()
goodfile_dir ='/home/kstory/lps12/jacks/jackgoodfiles/'
runlist_dir  = '/home/kstory/lps12/runlists/'

; loop over fields
for idx=0, 19 do begin
    
    new_runlist = runlist_dir+'runlist_azrms_'+skeep+'_'+f[idx].name+'.txt'
    f1 = goodfile_dir + 'goodfiles_'+f[idx].name+'_az_highrms_'+skeep+'.txt'
    f2 = goodfile_dir + 'goodfiles_'+f[idx].name+'_az_lowrms_'+skeep+'.txt'
    
    
    ; read jack files
    print, 'read f1: ' + f1
    readcol,/silent,f1,date_1,format='a'
    
    print, 'read f2: ' + f2
    readcol,/silent,f2,date_2,format='a'
    
    dates = [date_1, date_2]
    dates = dates[sort(dates)]
    nmaps = n_elements(dates)

    ; write new runlist
    print, 'write out new runlist: ' + new_runlist
    get_lun, lun1
    openw, lun1, new_runlist
    for ii=0, nmaps-1 do begin
        printf, lun1, dates[ii]
    endfor
    close, lun1
    free_lun,lun1
    
endfor


stop
END


; Make runlists with random 10% data cuts.
PRO rand_cut_runlists, seed=seed

; setup
keep_pct = 90
skeep = strtrim(string(keep_pct), 2)
if (n_elements(seed) eq 0) then seed = 1
sseed = strtrim(string(seed), 2)

; setup
f = lps12_fieldstruct()
runlist_dir  = '/home/kstory/lps12/runlists/'

; loop over fields
for idx=0, 19 do begin
    
    new_runlist = runlist_dir+'runlist_randcut_'+skeep+'_seed'+sseed+'_'+f[idx].name+'.txt'
    
    dates = get_lps12_runlist(idx, /xspec_dates)
    n_dates = n_elements(dates)

    ; Cut a random 10%
    rand_idx = randperm(n_dates, SEED=seed)
    rand_dates = dates[rand_idx]

    nkeep = floor((keep_pct/100.) * n_dates)
    new_dates = rand_dates[0:nkeep]
    new_dates = new_dates[sort(new_dates)]
    
    ; write new runlist
    print, 'write out new runlist: ' + new_runlist
    get_lun, lun1
    openw, lun1, new_runlist
    for ii=0, nkeep-1 do begin
        printf, lun1, new_dates[ii]
    endfor
    close, lun1
    free_lun,lun1
endfor

stop
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Run end2end on azrms_90 cut, in 2 batches
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; run fields 0 to 9
PRO end2end_0to9
tbegin = systime(0, /seconds)
for idx=0, 9 do begin
    print, 'run end2end_04, cut_name="90", field ', idx
    t0 = systime(0, /seconds)

    try_end2end_04, idx, run='04', /resume, /use_kweight, cut_name='azrms_90'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tbegin - tend)/60., ' minutes.'
END

; run fields 10 to 19
PRO end2end_10to19
; re-make clobbered end2end files first:
t0 = systime(0, /seconds)
try_end2end_04, 0, run='04', /resume, /use_kweight
print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'

t0 = systime(0, /seconds)
try_end2end_04, 10, run='04', /resume, /use_kweight
print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'

t0 = systime(0, /seconds)
try_end2end_04, 11, run='04', /resume, /use_kweight
print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'

; now run end2end on azrms_90
tbegin = systime(0, /seconds)
for idx=10, 19 do begin

    print, 'run end2end_04, cut_name="90", field ', idx
    t0 = systime(0, /seconds)

    try_end2end_04, idx, run='04', /resume, /use_kweight, cut_name='azrms_90'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tbegin - tend)/60., ' minutes.'
END




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Run end2end on randcut_90_seed1 and 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; seed1
PRO end2end_randcut_90_seed1
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    print, 'run end2end_04, cut_name="randcut_90_seed1", field ', idx
    t0 = systime(0, /seconds)

    try_end2end_04, idx, run='04', /resume, /use_kweight, cut_name='randcut_90_seed1'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tbegin - tend)/60., ' minutes.'
END

; seed2
PRO end2end_randcut_90_seed2
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    print, 'run end2end_04, cut_name="randcut_90_seed2", field ', idx
    t0 = systime(0, /seconds)

    try_end2end_04, idx, run='04', /resume, /use_kweight, cut_name='randcut_90_seed2'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tbegin - tend)/60., ' minutes.'
END

; seed3
PRO end2end_randcut_90_seed3
tbegin = long(systime(0, /seconds))
for idx=0, 19 do begin
    print, 'run end2end_04, cut_name="randcut_90_seed3", field ', idx
    t0 = systime(0, /seconds)

    try_end2end_04, idx, run='04', /resume, /use_kweight, cut_name='randcut_90_seed3'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = long(systime(0, /seconds))
print, 'The whole thing took ', (tbegin - tend)/60., ' minutes.'
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Run end2end on data in order to get all_mc_spectra
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; seed1
PRO end2end_all
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    print, 'run end2end_04, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_04, idx, run='04', /resume, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tbegin - tend)/60., ' minutes.'
END


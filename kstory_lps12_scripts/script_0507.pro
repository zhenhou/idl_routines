;;;
; NAME: script_0507
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) coadd remaining fields (17,18,19)
;
; MODIFICATION HISTORY:
;  05/07/2012: (KTS) Created
;;;

; Process a field
PRO process_field, idx
print, 'coadd maps, idx = ', idx
coadd_maps_0420, idx

print, 'make mask, idx = ', idx
make_mask_lps12, field_idx=idx

;print, 'make noise psd, idx = ', idx
;make_noise_psd_lps12, idx
;make_coupling_kernel, idx
;try_end2end, idx, run='03', /use_kweight, /resume
;make_twod_tfs, idx ??
;make_twod_kweights_lps12, idx ??
END

;;;;;;;;;;;;;;;;;;;;;;;;
; Processing for today
;;;;;;;;;;;;;;;;;;;;;;;;

; coadd all remaining maps
PRO coadd_all ; DONE
for ii=0, 19 do begin
    print, 'coadd maps, idx = ', ii
    coadd_maps_0420, ii
endfor
END
    
; make all apod masks
PRO masks_all ; DONE
print, 'masks, idx = ', 17 & make_mask_lps12, field_idx=17
print, 'masks, idx = ', 18 & make_mask_lps12, field_idx=18
print, 'masks, idx = ', 19 & make_mask_lps12, field_idx=19
END


PRO coupling_6789 ; DONE
for ii=6, 9 do begin
    print, 'make coupling kernel, field ', ii
    make_coupling_kernel, ii
endfor
END

PRO coupling_11121314 ; DONE
for ii=11, 14 do begin
    print, 'make coupling kernel, field ', ii
    make_coupling_kernel, ii
endfor
END

PRO coupling_1516171819 ; DONE
for ii=15, 19 do begin
    print, 'make coupling kernel, field ', ii
    make_coupling_kernel, ii
endfor
END

PRO noise_7to19 ; running (9:30pm)
for ii=7, 19 do begin
    print, 'make coupling kernel, field ', ii
    make_noise_psd_lps12, ii
endfor
END

PRO kweight_3 ; DONE
make_twod_kweights_lps12, 3, /save_files
END

PRO end_1 ; DONE
try_end2end, 1, run='03', /use_kweight, /resume
END

PRO process_2 ; running (9:30pm)
idx = 2
time_0 = systime(0,/seconds)

;print, 'Run try_end2end, idx=', idx
;try_end2end, idx, run='03', /resume
; time_1 = systime(0,/seconds)

; print, 'Run save_mode_coupling_lps12, idx=', idx
; save_mode_coupling_lps12, '03', idx_list=[-1,idx]
; time_2 = systime(0,/seconds)

print, 'Run make_twod_tfs, idx=', idx
make_twod_tfs, '03', /dosave, idx_list=[-1,idx]
time_3 = systime(0,/seconds)

print, 'Run make_twod_kweights, idx=', idx
make_twod_kweights_lps12, idx, /save_files
time_4 = systime(0,/seconds)

print, 'Re-Run try_end2end, idx=', idx
try_end2end, idx, run='03', /use_kweight, /resume
time_5 = systime(0,/seconds)

print, 'Run-times [min]:'
print, 'try_end2end, 1st go: ', (time_1-time_0) / 60.
print, 'save_mode_coupling_lps12: ', (time_2-time_1) / 60.
print, 'make_twod_tfs: ', (time_3-time_2) / 60.
print, 'make_twod_kweights: ', (time_4-time_3) / 60.
print, 'try_end2end, 2nd go: ', (time_5-time_4) / 60.

END




;;;;;;;;;;;;;;;;;;;;;;;;
; band powers
;;;;;;;;;;;;;;;;;;;;;;;;

; plot band powers from end2end
; Copied from combine_fields_lps12.pro
PRO check_bandpower
idx = 1
run = '03'

beam_err_scalefactor = 1.0
cal_error_power = 0.031

;----------------------
; Setup
;----------------------
f = lps12_fieldstruct()
idx_list =[0]
nfields = n_elements(idx_list)

reso=1.0

; directories
end_dir = '/home/kstory/lps12/end2end/'

; plotting
if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8

;----------------------
; Loop over all fields in idx_list, get dls
;----------------------

i = idx
field = f[i].name
fname = end_dir+'end_'+field+'_'+run
fname = fname + '_kweight.sav'
print,fname
restore,fname

; initialize
l = banddef-25.
nl = n_elements(l)
dls = dblarr(nfields,nl)
covs = dblarr(nfields,nl,nl)
covs_cond = dblarr(nfields,nl,nl)
corrs = dblarr(nfields,nl,nl)
uncorrs = dblarr(nfields,nl,nl)
diags = dblarr(nfields,nl)
diags2 = dblarr(nfields,nl)
area = dblarr(nfields)
years = fltarr(nfields)

nwf = n_elements(windowfunc[*,0])
wfs = dblarr(nfields,nwf,nl)

years = [2008., 2009., 2010., 2011.]

area[i] = total(mask_padded)*(reso/60.)^2.
cov = reform(cov)

wfs[i,*,*] = windowfunc

; Condition the covariance matrix
nbins_off_diag = 5
cov_cond = condition_cov_matrix(cov, nbins_off_diag, $
                                corr=this_corr, uncondcorr=this_uncorr, $
                                banddef=banddef, $
                                ellmin=ellmin, ellmax=ellmax, $
                                noaverage=noaverage,$
                                knowncorr=knowncorr)
    

covs_cond[i,*,*] = cov_cond
corrs[i,*,*] = this_corr
uncorrs[i,*,*] = this_uncorr

; Get the spectrum
dls[i,*]=spectrum
covs[i,*,*]=cov
these_diags = dblarr(nl)
for j=0,nl-1 do these_diags[j] = sqrt(cov[j,j])
diags[i,*] = these_diags

these_diags2 = dblarr(nl)
for j=0,nl-1 do these_diags2[j] = sqrt(cov_cond[j,j])
diags2[i,*] = these_diags2

if n_elements(hack_factor) eq 0 then hack_factor=1.
dls *= (hack_factor)

;----------------------
; create ell-dependent weights
;----------------------

w0 = dblarr(nfields,nl) ; weight array
for i=0,nfields-1 do w0[i,*] = area[i]
for i=0,nl-1 do w0[*,i] /= total(w0[*,i])

w1 = diags^(-2.)
for i=0,nl-1 do w1[*,i] /= total(w1[*,i])

w2 = diags2^(-2.)
for i=0,nl-1 do w2[*,i] /= total(w2[*,i])

w4 = dblarr(nfields,nl)
poly_order = 2
for i=0,nfields-1 do begin
    this_w = reform(w2[i,*]) ; let's smooth the weights based on the diagonals of the conditioned, beam-error-excluded cov matrices.
    whfit = where(finite(this_w),nfit)
    fit = poly_fit(banddef[whfit],this_w[whfit],poly_order)
    reco = dblarr(nl)
    for j=0,poly_order do reco += (fit[j]*banddef^(1d*j))
    w4[i,*] = reco
endfor
for i=0,nl-1 do w4[*,i] /= total(w4[*,i])

; let's use this ell-dependent weighting scheme:
;w2use = w4;ell-dependent, near-optimal
w2use = w0;area

;----------------------
; let's create the composite spectrum and window function
dl_all_init = dblarr(nl)
wf_all = dblarr(nwf,nl)
for i=0,nl-1 do begin
    dl_all_init[i] = total(dls[*,i]*w2use[*,i])
    tmp = dblarr(nwf)
    for j=0,nfields-1 do begin
        tmp += w2use[j,i]*wfs[j,*,i]
    endfor
    wf_all[*,i] = tmp
endfor
;l_wf = findgen(9000-50-1)+50
l_wf = findgen(4000-50-1)+50

;----------------------
; Composite Covariance Matricies
;----------------------

;----------------------
; 2008
year0 = 2008
whyear = where(years eq year0,nyear)
dl_2008 = dblarr(nl)
cov_2008 = dblarr(nl,nl)
for i=0,nl-1 do begin
    these_wi = reform(w2use[whyear,i])
    these_wi /= total(these_wi)
; let's create the composite spectrum
    dl_2008[i] = total(dls[whyear,i]*these_wi)
    for j=0,nl-1 do begin
        these_wj = reform(w2use[whyear,j])
        these_wj /= total(these_wj)
        cov_2008[i,j] = total(reform(covs_cond[whyear,i,j])*these_wi*these_wj)
    endfor
endfor
; now get the beam covariance for this year
get_beam_cov_matrix_lps12, l=l, year=year0, cov=beam_cov,scalefactor=beam_err_scalefactor
;tempp, let's add the calibration uncertainty to the beam cov
;beam_cov += (fltarr(nl,nl) + (0.072^2.))
;beam_cov += (fltarr(nl,nl) + (0.2^2.))
for j=0,nl-1 do begin
    for k=0,nl-1 do begin
        beam_cov[j,k] *= (dl_2008[j]*dl_2008[k])
    endfor
endfor
cov_2008_nobeam = cov_2008
cov_2008 += beam_cov

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; now combine the 2008 and 2009 DL and COV's.  We've already computed
; the composite spectrum, DL_ALL_INIT, so this will be a nice sanity check.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; first get the ell-dependent weight for each year.
nyears = 1
wyear = dblarr(nyears,nl)
wh2008 = where(years eq 2008)
for j=0,nl-1 do begin
    wyear[0,j] = total(w2use[wh2008,j])
endfor
for j=0,nl-1 do wyear[*,j] /= total(wyear[*,j])

dl_all = dblarr(nl)
cov_all = dblarr(nl,nl)
cov_all_nobeam = dblarr(nl,nl)
for i=0,nl-1 do begin
; let's create the composite spectrum
    dl_all[i] = dl_2008[i]*wyear[0,i]
    for j=0,nl-1 do begin
        cov_all_nobeam[i,j] = cov_2008_nobeam[i,j]*wyear[0,i]*wyear[0,j]

        cov_all[i,j] = cov_2008[i,j]*wyear[0,i]*wyear[0,j]
    endfor
endfor

; add calibration uncertainty to covariance matrix
cov_all_nocal = cov_all
for i=0,nl-1 do begin
    for j=0,nl-1 do begin
        cov_all[i,j] += ((cal_error_power^2.)*dl_all[i]*dl_all[j])
    endfor
endfor

; get diagonal errors
diag = dblarr(nl)
diag_nobeam = dblarr(nl)
for i=0,nl-1 do begin
    diag[i] = sqrt(cov_all[i,i])
    diag_nobeam[i] = sqrt(cov_all_nobeam[i,i])
endfor


dl2use = dl_all
diag2use = diag_nobeam
l2use = l

min_l_chisq = 700.
max_l_chisq = 3000.
;max_l_chisq = 1500.
;max_l_chisq = 2500.
;wh_chisq = where(l2use ge 500. and l2use le 3000., n_chisq)
wh_chisq = where(l2use ge min_l_chisq and l2use le max_l_chisq, n_chisq)


;whplot = where(l2use ge 500 and l2use le 3e3)
whplot = wh_chisq

;----------------------
; PLOTS
;----------------------

window,0
plot,l2use[whplot],dl2use[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=2, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars
errplot,l2use[whplot],(dl2use-diag2use)[whplot]*1d12,(dl2use+diag2use)[whplot]*1d12
clustered_and_sz = 10.
oplot,ellkern,SIMSPEC_INTERP+clustered_and_sz,color=!blue

oplot,ellkern,wmap7_th_cl(ellkern)*ellkern*(ellkern+1.)/2./!pi,color=!orange
legend,/top,/right,['Data','Sim Spectrum (WMAP7+Poisson+'+textoidl('10 \muK^2)'),'Sim Spectrum (WMAP7)'],textc=[!black,!blue,!orange],box=0

stop
END


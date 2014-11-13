;;;
; NAME: combine_fields_lps12.pro
; PURPOSE:
;   Combine spectra from fields into final bandpowers
;
; EXAMPLE CALL:
;   IDL> combine_fields_lps12, run='03', /use_kweight, idx_list=[0,1,3,4,5]
;   IDL> combine_fields_lps12, run='03', /use_kweight, /dosave
;
; INPUTS:
;   k11,                  Flag this to use only 2008, 2009 fields
;   
; NOTES:
;  1) output from try_end2end.pro
;
; TODO:
;
; MODIFICATION HISTORY:
;  05/04/2012: (KTS) Created from /home/rkeisler/ps09/combine_fields.pro
;  05/21/2012: (KTS) Fix bug.  Previously covariances for 2010 and 2011 were explicitly being set to the valuse for 2009.  Oops.
;  06/06/2012: (KTS) Add 'cut_name' option for azrms_90, randcut
;  06/08/2012: (KTS) Add by-year spectra to sav file
;  06/10/2012: (KTS) Add combined_calib
;  06/10/2012: (KTS) Use full cov_sv
;  06/11/2012: (KTS) Fix calibration, add k11 option
;  06/13/2012: (KTS) Fix calibration bug, use beam cov matrix from all years.
;  06/29/2012: (KTS) Change how covariance matricies are condidioned and combined.
;  07/16/2012: (KTS) Use new get_beam_correlation_matrix, add by-year
;  weights, remove all year calculations (now on get_beam_correlation_matrix).
;  07/17/2012: (KTS) move ellmax for conditioning upto 3100, so the last element does not get truncated to zero.
;  07/18/2012: (KTS) Fix K11 option to work with changes from 06/29
;  08/27/2012: (KTS) Fix bug - previously the cal error was not being included in the final cov matrix.
;  08/27/2012: (KTS) Change cal_error_power_in from 0.027 to 0.026
;  08/30/2012: (KTS) Add use1011 option for 2010, 2011 data only
;;;

PRO combine_fields_lps12, run=run, use_kweight=use_kweight, idx_list=idx_list, $
                          use_cov_sv=use_cov_sv, $
                          dosave=dosave, $  ; save combined spectrum
                          cut_name=cut_name, $
                          k11=k11, $
                          use1011=use1011, $
                          savename=savename_in, $
                          beam_err_scalefactor=beam_err_scalefactor_in, $
                          cal_error_power=cal_error_power_in, $
                          xr=xr, yr=yr, $
                          xtitle=xtitle, ytitle=ytitle


;----------------------
; Setup
;----------------------
f = lps12_fieldstruct()
if (n_elements(idx_list) eq 0) then idx_list = [0,1,indgen(17)+3]
if keyword_set(k11) then idx_list = [0,1,3,4,5] ; only 2008, 2009 data
if keyword_set(use1011) then idx_list = indgen(14)+6 ; only 2010, 2011 data
nfields = n_elements(idx_list)
reso=1.0

beam_err_scalefactor = n_elements(beam_err_scalefactor_in) gt 0 ? beam_err_scalefactor_in : 1.0
cal_error_power = n_elements(cal_error_power_in) gt 0 ? cal_error_power_in : 0.026 ;0.031

; Get the over-all calibration (Temperature units)
case run of
    '05' : combined_calib = 0.825 ; 1/1.47 = 0.680
    '08' : combined_calib = 0.825 ; 1/1.47 = 0.680
    '09' : combined_calib = 0.825 ; 1/1.47 = 0.680
    else : combined_calib = 1.0
endcase

; directories
end_dir = '/home/kstory/lps12/end2end/'
rundir  = '/home/kstory/lps12/end2end/run_'+run+'/'
if (run eq '03') then end_dir = '/home/kstory/lps12/end2end/run_03/'

; plotting
if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8

;----------------------
; Loop over all fields in idx_list, get dls
;----------------------

for i=0,nfields-1 do begin
    idx = idx_list[i]
    field = f[idx].name
    fname = end_dir+'end_'+field+'_'+run
    if keyword_set(use_kweight) then fname = fname + '_kweight'
    if keyword_set(cut_name) then fname = fname + '_' + cut_name
    fname = fname + '.sav'
    print,idx, ',  '+fname
    restore,fname
    if i eq 0 then begin
        l = banddef-25.
        nl = n_elements(l)
        dls = dblarr(nfields,nl)
        area = dblarr(nfields)
        years = fltarr(nfields)

        nwf = n_elements(windowfunc[*,0])
        wfs = dblarr(nfields,nwf,nl)
        wfs_sim = dblarr(nfields,nwf,nl)

        ; covariance matricies
        covs_mc        = dblarr(nfields,nl,nl)
        covs_data      = dblarr(nfields,nl,nl)
        covs_e2e       = dblarr(nfields,nl,nl) ; consistency checks only
    endif

    years[i] = f[idx].year
    area[i] = total(mask_padded)*(reso/60.)^2.
    wfs[i,*,*] = windowfunc
    wfs_sim[i,*,*] = sim_windowfunc

    ; Get the spectrum
    dls[i,*]=spectrum

    ; Get the covariance matricies
    covs_mc[i,*,*]=reform(sample_cov)
    covs_data[i,*,*]=reform(meas_cov)
    covs_e2e[i,*,*]=reform(cov)

endfor

;----------------------
; Sort out year stuff
;----------------------
yy=[2008, 2009, 2010, 2011]
nuniq_year = n_elements(yy)

wh2008 = where(years eq 2008, n08)
wh2009 = where(years eq 2009, n09)
wh2010 = where(years eq 2010, n10)
wh2011 = where(years eq 2011, n11)
wh0809 = where( (years eq 2008) or (years eq 2009), n0809)

;----------------------
; create area-based weights
;----------------------

w0 = dblarr(nfields,nl) ; weight array
for i=0,nfields-1 do w0[i,*] = area[i]
for i=0,nl-1 do w0[*,i] /= total(w0[*,i])
w2use = w0;area

; Get the weight for each year.
wyear = dblarr(nuniq_year)
wyear[0] = (n08 ne 0) ? total(w2use[wh2008]) : 0
wyear[1] = (n09 ne 0) ? total(w2use[wh2009]) : 0
wyear[2] = (n10 ne 0) ? total(w2use[wh2010]) : 0
wyear[3] = (n11 ne 0) ? total(w2use[wh2011]) : 0



;----------------------
; let's create the composite spectrum and window function
;----------------------
dl_all = dblarr(nl)
cov_all             = dblarr(nl,nl)
cov_all_nobeam      = dblarr(nl,nl)
cov_all_mc          = dblarr(nl,nl)
cov_all_data        = dblarr(nl,nl)
cov_all_e2e         = dblarr(nl,nl)
wf_all = dblarr(nwf,nl)
wf_all_sim = dblarr(nwf,nl)

; make the composite spectrum and window function
for i=0,nl-1 do begin
    dl_all[i] = total(dls[*,i]*w2use[*,i]) * (combined_calib^2.)

    for j=0, nwf-1 do wf_all[j,i] = total( wfs[*,j,i] * w2use[*,i] )
    for j=0, nwf-1 do wf_all_sim[j,i] = total( wfs_sim[*,j,i] * w2use[*,i] )
endfor
l_wf = findgen(4000-50-1)+50

; make the composite covariance matricies
for ifield=0, nfields-1 do begin
    for jbin=0,nl-1 do begin
        for kbin=0, nl-1 do begin
            cov_all_mc[jbin, kbin]   += covs_mc[ifield,jbin,kbin] * w2use[ifield,jbin]*w2use[ifield,kbin]
            cov_all_data[jbin, kbin] += covs_data[ifield,jbin,kbin] * w2use[ifield,jbin]*w2use[ifield,kbin] * (combined_calib^4.)
            cov_all_e2e[jbin, kbin]  += covs_e2e[ifield,jbin,kbin] * w2use[ifield,jbin]*w2use[ifield,kbin]
        endfor
    endfor
endfor


;----------------------
; add Sample Variance to covariance matrix
;----------------------
restore, end_dir+'cov_sv_'+run+'_kweight.sav' ; gets cov_sv
; Use simulated sample variance
if keyword_set(use_cov_sv) then begin
    cov_all        = cov_all_data + cov_sv

    ; 0809 only
    if keyword_set(k11) then begin
        cov_all        *= 0.    ; re-set variables
        for j=0, nl-1 do begin
            for k=0, nl-1 do begin
                cov_all[j,k] = cov_all_data[j,k] + total( cov_sv_y[*,j,k]*wyear[*]^2. )
            endfor
        endfor
    endif

    ; 1011 only
    if keyword_set(use1011) then begin
        cov_all        *= 0.    ; re-set variables
        for j=0, nl-1 do begin
            for k=0, nl-1 do begin
                cov_all[j,k] = cov_all_data[j,k] + total( cov_sv_y[*,j,k]*wyear[*]^2. )
            endfor
        endfor
    endif

; or just use the output from end2end
endif else begin
    cov_all        = cov_all_e2e
endelse


;----------------------
; Condition the combined cov
;----------------------
cov_all_uncond_nobeam = cov_all

nbins_off_diag = 5
ellmin=650
ellmax=3100

; condition
cov_all = condition_cov_matrix(cov_all, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)


;----------------------
; Get beam covariances
;----------------------
cov_all_nobeam = cov_all ; no beam, with conditioning
beamdir = '/data/sptdat/beams/rev3.1/'

; Note: this still works for 0809 only, because wyear=0 for 2010, 2011
beam_cor = get_beam_correlation_matrix(l=l, year=yy, weight_years=wyear, band1=150, band2=150, scalefactor=beam_err_scalefactor, beamdir=beamdir)
beamcal_cor = beam_cor + (cal_error_power)^2

; add calibration uncertainty and beam error
bccov = beam_cor*0.
for i=0, nbands-1 do begin
    for j=0, nbands-1 do begin
        bccov[i,j] = beamcal_cor[i,j]*dl_all[i]*dl_all[j]
    endfor
endfor

cov_all += bccov

;----------------------
; get diagonal errors
;----------------------
diag = dblarr(nl)
diag_nobeam = dblarr(nl)
for j=0,nl-1 do begin
    diag[j] = sqrt(cov_all[j,j])
    diag_nobeam[j] = sqrt(cov_all_nobeam[j,j])
endfor


;----------------------
; Set which variables to use
;----------------------
dl2use = dl_all
diag2use = diag_nobeam
l2use = l

min_l_chisq = 700.
max_l_chisq = 3000.
wh_chisq = where(l2use ge min_l_chisq and l2use le max_l_chisq, n_chisq)


;whplot = where(l2use ge 500 and l2use le 3e3)
whplot = wh_chisq

;----------------------
; PLOTS
;----------------------

; get the input theory spectrum
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)

window,0
plot,l2use[whplot],dl2use[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars
oplot,l_vec,dl_th,color=!red
errplot,l2use[whplot],(dl2use-diag2use)[whplot]*1d12,(dl2use+diag2use)[whplot]*1d12
;clustered_and_sz = 10.
;oplot,ellkern,SIMSPEC_INTERP+clustered_and_sz,color=!blue

;oplot,ellkern,wmap7_th_cl(ellkern)*ellkern*(ellkern+1.)/2./!pi,color=!orange
;legend,/top,/right,['Data','Sim Spectrum (WMAP7+Poisson+'+textoidl('10 \muK^2)'),'Sim Spectrum (WMAP7)'],textc=[!black,!blue,!orange],box=0

;----------------------
; Additional plots, skip for now
;----------------------

if 0 then begin ; skip this for now

; get ACBAR                             
color_acbar=!darkgreen
readcol,'acbar_reichardt.txt',band_num_acbar,dl_acbar,ddl_acbar,log_offset_acbar
l_acbar = [470., 608., 694., 763., 823., 884., 944., 1003., 1062., 1122., 1182., 1242., 1301., 1361., 1421., 1481., 1541., 1618., 1713., 1814., 1898., 2020., 2194., 2391., 2646.]
;oplot,l_acbar,dl_acbar,color=color_acbar,ps=2
;errplot,l_acbar,(dl_acbar-ddl_acbar),(dl_acbar+ddl_acbar),color=color_acbar

window,1

dl_th_wf = (wmap7_th_cl(l_wf) + (7.36d-6) + (5.21003d-6))*l_wf*(l_wf+1.)/2./!pi + clustered_and_sz
dl_th_wf_nolens = (wmap7_th_cl(l_wf,/nolensing) + (7.36d-6) + (5.21003d-6))*l_wf*(l_wf+1.)/2./!pi + clustered_and_sz

dl_th = dblarr(nl)
dl_th_nolens = dblarr(nl)
for j=0,nl-1 do begin
    dl_th[j] = total(dl_th_wf*wf_all[*,j])
    dl_th_nolens[j] = total(dl_th_wf_nolens*wf_all[*,j])
endfor

ratio = dl2use/dl_th*1d12
ratio_nolens = dl2use/dl_th_nolens*1d12
err_ratio = diag2use/dl_th*1d12

;ratio = dl2use/interpol(SIMSPEC_INTERP+clustered_and_sz,ellkern,l2use)*1d12
;err_ratio = diag2use/interpol(SIMSPEC_INTERP+clustered_and_sz,ellkern,l2use)*1d12

mean_ratio = mean(ratio[wh_chisq])
print,'MEAN_RATIO: ',mean_ratio
ratio /= mean_ratio
chisq = total((((ratio-1.)/err_ratio)^2.)[wh_chisq])
dof = n_chisq-2.
pte = mpchitest(chisq,dof)
nsigma = mpchitest(chisq,dof,/sigma)


icov = invert(cov_all[min(wh_chisq):max(wh_chisq), min(wh_chisq):max(wh_chisq)], status, /double)

;tempp
delta = (dl_all - dl_th/1d12)[wh_chisq]
;delta = (dl_all/mean((dl_all/dl_th*1d12)[wh_chisq]) - dl_th/1d12)[wh_chisq]
chisq_full = matrix_multiply(transpose(delta),matrix_multiply(icov,delta))
dof = n_chisq-2.
pte_full = mpchitest(chisq_full,dof)
nsigma_full = mpchitest(chisq_full,dof,/sigma)


delta = (dl_all/mean((dl_all/dl_th_nolens*1d12)[wh_chisq]) - dl_th_nolens/1d12)[wh_chisq]
chisq_full_nolens = matrix_multiply(transpose(delta),matrix_multiply(icov,delta))
dof = n_chisq-2.
pte_full_nolens = mpchitest(chisq_full_nolens,dof)
nsigma_full_nolens = mpchitest(chisq_full_nolens,dof,/sigma)

print,' '
print,'********************************'
print,'PTE_FULL: ',pte_full
print,'PTE_FULL_NOLENS: ',pte_full_nolens
lensing_sn = sqrt(chisq_full_nolens-chisq_full)
print,'Nominal Lensing S/N: ',lensing_sn
print,'********************************'
print,' '
plot,l2use[whplot],ratio[whplot],yr=[0.5,1.5],/yst,xr=xr,ps=2, /xst, $
  xtitle='ELL',ytitle='CL/CL_TH',title='PTE: '+sigfig(pte_full,3)+', '+sigfig(nsigma_full,3)+textoidl('\sigma'),chars=1.7
oplot,l2use[whplot],ratio[whplot]
errplot,l2use[whplot],(ratio-err_ratio)[whplot],(ratio+err_ratio)[whplot]
oplot,[1,1]*min_l_chisq,[-1,1]*9999.,lines=1
oplot,[1,1]*max_l_chisq,[-1,1]*9999.,lines=1
oplot,[-1,1]*9999999.,[1,1]

window,2
plot,[0],[0],xr=xr,yr=[0,.75], /yst, /xst, $
  xtitle='!12l!X!N', $
  ytitle='Weight'
colors = [!blue,!red,!darkgreen,!orange,!purple]
for i=0,4 do oplot,l2use,w0[i,*],color=colors[i],lines=2
for i=0,4 do oplot,l2use,w1[i,*],color=colors[i]
for i=0,4 do oplot,l2use,w2[i,*],color=colors[i]
;for i=0,4 do oplot,l2use,w3[i,*],color=colors[i]
for i=0,4 do oplot,l2use,w4[i,*],color=colors[i],thick=2

legend,/top,/left,strtrim(field_vec(),2),textc=colors,box=0,chars=1.7

endif 
;------------------------------

;stop
if keyword_set(dosave) then begin
    if n_elements(savename_in) eq 0 then begin
        tag = date_toolkit('now','file')
        savename = rundir+'combined_spectrum_'+tag
        if keyword_set(use_kweight) then savename = savename + '_kweight'
        if beam_err_scalefactor ne 1.0 then savename = savename + '_beamErr'
        if keyword_set(cut_name) then savename = savename + '_'+cut_name
        if keyword_set(k11) then savename = savename + '_0809'
        if keyword_set(use1011) then savename = savename + '_1011'
        savename = savename + '.sav'
    endif else savename=savename_in

    help_str = 'run_'+run+', combined_calib applied to dl_all, diag*, cov_all*.  dl in [K^2].'

    print,'Sure you want to save?'
    print,'We would be writing to: ',savename
    pause
    save,l,dl_all,diag,diag_nobeam, $
      cov_all,cov_all_nobeam, cov_all_uncond_nobeam, $
      cov_all_mc, $
      cov_all_data, $
      cov_all_e2e, $
      cov_sv, cov_sv_y, $
      wyear, w0, area, $
      cal_error_power, beam_err_scalefactor, $
      l_wf,wf_all, wf_all_sim, $
      calib, combined_calib, $
      help_str, $
      filename=savename
endif


stop
end


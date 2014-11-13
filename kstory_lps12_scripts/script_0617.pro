;;;
; NAME: script_0617
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Plots
;
; MODIFICATION HISTORY:
;  06/14/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;
; Plot error bars
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
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_205039_kweight_0809.sav' ; from comb_0809
dl0_0809     = dl_all
diag0_0809 = diag_nobeam

restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_205039_kweight_0809.sav' ; from comb2_0809
;dl2_0809     = dl_all_init
diag2_0809 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_rk = diag_nobeam

; cl_theory
readcol,'/home/kstory/lps12/theory_spectrum/Cls_100sims_ave_from_combined_map.txt',l_vec,cl_uK2
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
; ; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0617')

; ; compare against K11
window, 2
plot, l[whplot], ((err_rk )/diag0_0809)[whplot], xtitle=xtitle, ytitle='diag_k11/diag_0809', title='Error bars, comb_0809'
oplot, [0,10000], [1,1], linestyle=1
err=tvread(/png,/nodialog,filename=fdir+'err_k11v0809_0617')

;window, 3
;plot, l[whplot], ((err_rk )/diag2_0809)[whplot], xtitle=xtitle, ytitle='diag_k11/comb_0809', title='Error bars, comb_0809'

; window, 3
; plot, l[whplot], (dl_05/dl_07)[whplot], xtitle=xtitle, ytitle='dl_05 / dl_07'

; window, 4
; plot, l[whplot], (diag05/diag07)[whplot], xtitle=xtitle, ytitle='dl_05 / dl_07'


stop
END






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Idiot-check combining fields
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO comb_0809

; input arguments
use_kweight=1
run='07'
k11=1
dosave=1
idx_list = [0,1,3,4,5] ; only 2008, 2009 data

beam_err_scalefactor = n_elements(beam_err_scalefactor_in) gt 0 ? beam_err_scalefactor_in : 1.0
cal_error_power = n_elements(cal_error_power_in) gt 0 ? cal_error_power_in : 0.027 ;0.031

;----------------------
; Setup
;----------------------
f = lps12_fieldstruct()
nfields = n_elements(idx_list)

reso=1.0

; Get the over-all calibration (Temperature units)
case run of
    '05' : combined_calib = 0.825 ; 1/1.47 = 0.680
    else : combined_calib = 1.0
endcase

; DEBUGGING
;combined_calib = 1.0

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
;         covs = dblarr(nfields,nl,nl) ; covs_data
;         covs_cond = dblarr(nfields,nl,nl)
;         corrs = dblarr(nfields,nl,nl)
;         uncorrs = dblarr(nfields,nl,nl)
;         diags = dblarr(nfields,nl)
;         diags2 = dblarr(nfields,nl)
        area = dblarr(nfields)
        years = fltarr(nfields)

        nwf = n_elements(windowfunc[*,0])
        wfs = dblarr(nfields,nwf,nl)
        wfs_sim = dblarr(nfields,nwf,nl)

        ; covariance matricies
        covs_mc        = dblarr(nfields,nl,nl)
        covs_cond_mc   = dblarr(nfields,nl,nl)
        covs_data      = dblarr(nfields,nl,nl)
        covs_cond_data = dblarr(nfields,nl,nl)
        covs_e2e       = dblarr(nfields,nl,nl)
        covs_cond_e2e  = dblarr(nfields,nl,nl)
    endif

    years[i] = f[idx].year
    area[i] = total(mask_padded)*(reso/60.)^2.
    wfs[i,*,*] = windowfunc
    wfs_sim[i,*,*] = sim_windowfunc



    ; Condition the covariance matricies
    nbins_off_diag = 5
    sample_cov = reform(sample_cov)
    cov_cond_mc = condition_cov_matrix(sample_cov, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)
    covs_cond_mc[i,*,*] = cov_cond_mc
    ;corrs[i,*,*] = this_corr 
    ;uncorrs[i,*,*] = this_uncorr

    meas_cov = reform(meas_cov)
    cov_cond_data = condition_cov_matrix(meas_cov, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)
    covs_cond_data[i,*,*] = cov_cond_data

    cov_e2e = reform(cov)
    cov_cond_e2e = condition_cov_matrix(cov_e2e, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)
    covs_cond_e2e[i,*,*] = cov_cond_e2e
    ; DEBUGGING --------------------------

    ; Get the spectrum
    dls[i,*]=spectrum

    covs_mc[i,*,*]=sample_cov
    covs_data[i,*,*]=meas_cov
    covs_e2e[i,*,*]=cov_e2e

endfor

;----------------------
; Sort out year stuff
;----------------------
yy=[2008, 2009]
nuniq_year = n_elements(yy)

wh2008 = where(years eq 2008, n08)
wh2009 = where(years eq 2009, n09)

;----------------------
; create ell-dependent weights
;----------------------

w0 = dblarr(nfields,nl) ; weight array
for i=0,nfields-1 do w0[i,*] = area[i]
for i=0,nl-1 do w0[*,i] /= total(w0[*,i])
w2use = w0;area

;----------------------
; let's create the composite spectrum and window function
dl_all_init = dblarr(nl)
wf_all = dblarr(nwf,nl)
wf_all_sim = dblarr(nwf,nl)
for i=0,nl-1 do begin
    dl_all_init[i] = total(dls[*,i]*w2use[*,i]) * (combined_calib^2.)

    for j=0, nwf-1 do wf_all[j,i] = total( wfs[*,j,i] * w2use[*,i] )
    for j=0, nwf-1 do wf_all_sim[j,i] = total( wfs_sim[*,j,i] * w2use[*,i] )
endfor
l_wf = findgen(4000-50-1)+50


;----------------------
; Composite Covariance Matricies
;----------------------

; let's create the composite covariance matrices.  one subtlety here
; is that we'd like to apply the 2008 beam covariance matrix to the
; 5h+23h covariance matrix, and only apply it once.  same goes with
; the 2009 fields.

dl_y              = dblarr(nuniq_year, nl)
wf_y              = dblarr(nuniq_year, nwf, nl)
cov_y_mc          = dblarr(nuniq_year, nl, nl)
cov_y_nobeam_mc   = dblarr(nuniq_year, nl, nl)
cov_y_data        = dblarr(nuniq_year, nl, nl)
cov_y_nobeam_data = dblarr(nuniq_year, nl, nl)
cov_y_e2e         = dblarr(nuniq_year, nl, nl)
cov_y_nobeam_e2e  = dblarr(nuniq_year, nl, nl)

cov_y_final         = dblarr(nuniq_year, nl, nl)
cov_y_nobeam_final  = dblarr(nuniq_year, nl, nl)

for iy=0, nuniq_year-1 do begin
    this_y = yy[iy]
    whyear = where(years eq this_y, nyear)

    ; if no fields from this year were included in idx_list, then skip
    if (nyear eq 0) then continue

    ; loop over ell bins
    for jbin=0, nl-1 do begin
        these_wj = reform(w2use[whyear, jbin])
        these_wj /= total(these_wj)
        dl_y[iy, jbin] = total(dls[whyear, jbin]*these_wj)
        for knwf=0, nwf-1 do wf_y[iy, knwf, jbin] = total(wfs[whyear, knwf, jbin]*these_wj)


        ; loop over bins again for cov matrix
        for kbin=0, nl-1 do begin
            these_wk = reform(w2use[whyear, kbin])
            these_wk /= total(these_wk)
            cov_y_mc[iy, jbin, kbin]   = total(reform(covs_cond_mc[whyear,jbin,kbin])*these_wj*these_wk)
            cov_y_data[iy, jbin, kbin] = total(reform(covs_cond_data[whyear,jbin,kbin])*these_wj*these_wk)
            cov_y_e2e[iy, jbin, kbin]  = total(reform(covs_cond_e2e[whyear,jbin,kbin])*these_wj*these_wk)
        endfor
    endfor

    ; now get the beam covariance for this year
    undefine, beam_cov          ; paranoia (0522)
    get_beam_cov_matrix_lps12, l=l, year=this_y, cov=beam_cov,scalefactor=beam_err_scalefactor
    for jbin=0,nl-1 do begin
        for kbin=0,nl-1 do begin
            beam_cov[jbin,kbin] *= (dl_y[iy, jbin]*dl_y[iy, kbin])
        endfor
    endfor

    ;----------------------
    ; add calibration factor to dl and cov_data
    ;----------------------
    dl_y[iy,*]                *= (combined_calib^2.)
    cov_y_data[iy,*,*]        *= (combined_calib^4.)
    cov_y_nobeam_data[iy,*,*] *= (combined_calib^4.)

    ; store the composite covariance matricies in the arrays.
    cov_y_nobeam_mc[iy,*,*]  = cov_y_mc[iy,*,*]
    cov_y_mc[iy,*,*] += beam_cov
    cov_y_nobeam_data[iy,*,*]  = cov_y_data[iy,*,*]
    cov_y_data[iy,*,*] += beam_cov
    cov_y_nobeam_e2e[iy,*,*]  = cov_y_e2e[iy,*,*]
    cov_y_e2e[iy,*,*] += beam_cov
    
endfor ; i nuniq_year

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; now combine the DL and COV's for all years.  We've already computed
; the composite spectrum, DL_ALL_INIT, so this will be a nice sanity check.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; first get the ell-dependent weight for each year.
wyear = dblarr(nuniq_year,nl)
for j=0,nl-1 do begin
    wyear[0,j] = (n08 ne 0) ? total(w2use[wh2008,j]) : 0
    wyear[1,j] = (n09 ne 0) ? total(w2use[wh2009,j]) : 0
endfor

; normalize to be 1 in each ell bin
for j=0,nl-1 do wyear[*,j] /= total(wyear[*,j])

;----------------------
; make the composite spectrum
;----------------------
dl_all = dblarr(nl)
wf_all_check = dblarr(nwf, nl) ; DEBUGGING
cov_all             = dblarr(nl,nl)
cov_all_nobeam      = dblarr(nl,nl)

cov_all_mc          = dblarr(nl,nl)
cov_all_nobeam_mc   = dblarr(nl,nl)
cov_all_data        = dblarr(nl,nl)
cov_all_nobeam_data = dblarr(nl,nl)
cov_all_e2e         = dblarr(nl,nl)
cov_all_nobeam_e2e  = dblarr(nl,nl)

for iy=0, nuniq_year-1 do begin

    this_y = yy[iy]
    whyear = where(years eq this_y, nyear)

    ; if no fields from this year were included in idx_list, then skip
    if (nyear eq 0) then continue

    print, yy[iy], ', nyear = ', nyear
    for jbin=0,nl-1 do begin
        ; let's create the composite spectrum
        dl_all[jbin] += dl_y[iy, jbin]*wyear[iy, jbin]
        for knwf=0, nwf-1 do wf_all_check[knwf, jbin] += wf_y[iy,knwf,jbin]*wyear[iy,jbin] ; DEBUGGING

        ; Create composite covariance matricies
        for kbin=0, nl-1 do begin
            cov_all_nobeam_mc[jbin, kbin] += cov_y_nobeam_mc[iy,jbin,kbin] * wyear[iy,jbin]*wyear[iy,kbin]
            cov_all_mc[jbin, kbin]        += cov_y_mc[iy,jbin,kbin] * wyear[iy,jbin]*wyear[iy,kbin]
            cov_all_nobeam_data[jbin, kbin] += cov_y_nobeam_data[iy,jbin,kbin] * wyear[iy,jbin]*wyear[iy,kbin]
            cov_all_data[jbin, kbin]        += cov_y_data[iy,jbin,kbin] * wyear[iy,jbin]*wyear[iy,kbin]
            cov_all_nobeam_e2e[jbin, kbin] += cov_y_nobeam_e2e[iy,jbin,kbin] * wyear[iy,jbin]*wyear[iy,kbin]
            cov_all_e2e[jbin, kbin]        += cov_y_e2e[iy,jbin,kbin] * wyear[iy,jbin]*wyear[iy,kbin]
        endfor
    endfor

endfor

;----------------------
; add Sample Variance to covariance matrix
;----------------------
cov_all        = cov_all_e2e
cov_all_nobeam = cov_all_nobeam_e2e

; Add sample variance to cov_y_final
cov_y_nobeam_final = cov_y_nobeam_e2e
cov_y_final = cov_y_e2e

;----------------------
; add calibration uncertainty to covariance matrix
;----------------------
cov_all_nocal = cov_all
for j=0,nl-1 do begin
    for jj=0,nl-1 do begin
        cov_all[j,jj] += ((cal_error_power^2.)*dl_all[j]*dl_all[jj])
    endfor
endfor

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


stop
if keyword_set(dosave) then begin
    if n_elements(savename_in) eq 0 then begin
        tag = date_toolkit('now','file')
        savename = rundir+'combined_spectrum_'+tag
        if keyword_set(use_kweight) then savename = savename + '_kweight'
        if keyword_set(cut_name) then savename = savename + '_'+cut_name
        if keyword_set(k11) then savename = savename + '_0809'
        savename = savename + '.sav'
    endif else savename=savename_in

    help_str = 'run_'+run+', combined_calib applied to dl_all, diag*, cov_all*.  dl in [K^2].'

    print,'Sure you want to save?'
    print,'We would be writing to: ',savename
    pause
    save,l,dl_all,diag,diag_nobeam, $
      cov_all,cov_all_nobeam, cov_all_nocal, $
      cov_all_mc,cov_all_nobeam_mc, $
      cov_all_data,cov_all_nobeam_data, $
      cov_all_e2e,cov_all_nobeam_e2e, $
      cov_sv, cov_sv_y, $
      yy, dl_y, wf_y, cov_y_mc, cov_y_nobeam_mc, cov_y_data, cov_y_nobeam_data, cov_y_e2e, cov_y_nobeam_e2e, cov_y_final, cov_y_nobeam_final, $
      wyear, w0, area, $
      cal_error_power, beam_err_scalefactor, $
      l_wf,wf_all, wf_all_sim, $
      calib, combined_calib, $
      help_str, $
      filename=savename
endif
stop
END




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Second Idiot-check combining fields
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO comb2_0809

; input arguments
use_kweight=1
run='07'
k11=1
dosave=1
idx_list = [0,1,3,4,5] ; only 2008, 2009 data

beam_err_scalefactor = n_elements(beam_err_scalefactor_in) gt 0 ? beam_err_scalefactor_in : 1.0
cal_error_power = n_elements(cal_error_power_in) gt 0 ? cal_error_power_in : 0.027 ;0.031

;----------------------
; Setup
;----------------------
f = lps12_fieldstruct()
nfields = n_elements(idx_list)

reso=1.0

; Get the over-all calibration (Temperature units)
case run of
    '05' : combined_calib = 0.825 ; 1/1.47 = 0.680
    else : combined_calib = 1.0
endcase

; DEBUGGING
;combined_calib = 1.0

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
;         covs = dblarr(nfields,nl,nl) ; covs_data
;         covs_cond = dblarr(nfields,nl,nl)
;         corrs = dblarr(nfields,nl,nl)
;         uncorrs = dblarr(nfields,nl,nl)
;         diags = dblarr(nfields,nl)
;         diags2 = dblarr(nfields,nl)
        area = dblarr(nfields)
        years = fltarr(nfields)

        nwf = n_elements(windowfunc[*,0])
        wfs = dblarr(nfields,nwf,nl)
        wfs_sim = dblarr(nfields,nwf,nl)

        ; covariance matricies
        covs_mc        = dblarr(nfields,nl,nl)
        covs_cond_mc   = dblarr(nfields,nl,nl)
        covs_data      = dblarr(nfields,nl,nl)
        covs_cond_data = dblarr(nfields,nl,nl)
        covs_e2e       = dblarr(nfields,nl,nl)
        covs_cond_e2e  = dblarr(nfields,nl,nl)
    endif

    years[i] = f[idx].year
    area[i] = total(mask_padded)*(reso/60.)^2.
    wfs[i,*,*] = windowfunc
    wfs_sim[i,*,*] = sim_windowfunc



    ; Condition the covariance matricies
    nbins_off_diag = 5
    sample_cov = reform(sample_cov)
    cov_cond_mc = condition_cov_matrix(sample_cov, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)
    covs_cond_mc[i,*,*] = cov_cond_mc
    ;corrs[i,*,*] = this_corr 
    ;uncorrs[i,*,*] = this_uncorr

    meas_cov = reform(meas_cov)
    cov_cond_data = condition_cov_matrix(meas_cov, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)
    covs_cond_data[i,*,*] = cov_cond_data

    cov_e2e = reform(cov)
    cov_cond_e2e = condition_cov_matrix(cov_e2e, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)
    covs_cond_e2e[i,*,*] = cov_cond_e2e
    ; DEBUGGING --------------------------

    ; Get the spectrum
    dls[i,*]=spectrum

    covs_mc[i,*,*]=sample_cov
    covs_data[i,*,*]=meas_cov
    covs_e2e[i,*,*]=cov_e2e

endfor

;----------------------
; create ell-dependent weights
;----------------------

w0 = dblarr(nfields,nl) ; weight array
for i=0,nfields-1 do w0[i,*] = area[i]
for i=0,nl-1 do w0[*,i] /= total(w0[*,i])
w2use = w0;area

;----------------------
; let's create the composite spectrum and window function
dl_all_init = dblarr(nl)
wf_all = dblarr(nwf,nl)
wf_all_sim = dblarr(nwf,nl)

cov_all_nobeam      = dblarr(nl,nl)
cov_all_nobeam_mc   = dblarr(nl,nl)
cov_all_nobeam_data = dblarr(nl,nl)
cov_all_nobeam_e2e  = dblarr(nl,nl)

for i=0,nl-1 do begin

    ; bandpowers
    dl_all_init[i] = total(dls[*,i]*w2use[*,i]) * (combined_calib^2.)

    ; window functions
    for j=0, nwf-1 do wf_all[j,i] = total( wfs[*,j,i] * w2use[*,i] )
    for j=0, nwf-1 do wf_all_sim[j,i] = total( wfs_sim[*,j,i] * w2use[*,i] )

    ; covariance matricies
    for j=0,nl-1 do begin
        for field_itr=0, nfields-1 do begin
            cov_all_nobeam_mc[i,j] += reform(covs_cond_mc[field_itr,i,j]) * w0[field_itr,i] * w0[field_itr,j]
            cov_all_nobeam_data[i,j] += reform(covs_cond_data[field_itr,i,j]) * w0[field_itr,i] * w0[field_itr,j]
            cov_all_nobeam_e2e[i,j] += reform(covs_cond_e2e[field_itr,i,j]) * w0[field_itr,i] * w0[field_itr,j]
        endfor
    endfor

endfor
l_wf = findgen(4000-50-1)+50


;----------------------
; add Sample Variance to covariance matrix
;----------------------
cov_all_nobeam = cov_all_nobeam_e2e

;----------------------
; add calibration uncertainty to covariance matrix
;----------------------
; cov_all_nocal = cov_all
; for j=0,nl-1 do begin
;     for jj=0,nl-1 do begin
;         cov_all[j,jj] += ((cal_error_power^2.)*dl_all[j]*dl_all[jj])
;     endfor
; endfor

;----------------------
; get diagonal errors
;----------------------
diag_nobeam = dblarr(nl)
for j=0,nl-1 do begin
    diag_nobeam[j] = sqrt(cov_all_nobeam[j,j])
endfor

;----------------------
; Set which variables to use
;----------------------
dl2use = dl_all_init
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

if keyword_set(dosave) then begin
    if n_elements(savename_in) eq 0 then begin
        tag = date_toolkit('now','file')
        savename = rundir+'combined_spectrum_'+tag
        if keyword_set(use_kweight) then savename = savename + '_kweight'
        if keyword_set(cut_name) then savename = savename + '_'+cut_name
        if keyword_set(k11) then savename = savename + '_0809'
        savename = savename + '.sav'
    endif else savename=savename_in

    help_str = 'run_'+run+', combined_calib applied to dl_all, diag*, cov_all*.  dl in [K^2].'

    print,'Sure you want to save?'
    print,'We would be writing to: ',savename
    pause
    save,l,dl_all_init,diag,diag_nobeam, $
      cov_all,cov_all_nobeam, cov_all_nocal, $
      cov_all_mc,cov_all_nobeam_mc, $
      cov_all_data,cov_all_nobeam_data, $
      cov_all_e2e,cov_all_nobeam_e2e, $
      cov_sv, cov_sv_y, $
      yy, dl_y, wf_y, cov_y_mc, cov_y_nobeam_mc, cov_y_data, cov_y_nobeam_data, cov_y_e2e, cov_y_nobeam_e2e, cov_y_final, cov_y_nobeam_final, $
      wyear, w0, area, $
      cal_error_power, beam_err_scalefactor, $
      l_wf,wf_all, wf_all_sim, $
      calib, combined_calib, $
      help_str, $
      filename=savename
endif
stop
END





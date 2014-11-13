;;;
; NAME: pipeline_test3.pro
; PURPOSE:
;   Analyze sims with different input spectrum.
;
; NOTES:
; 1) Pipeline test3
;;;

; Write the theory spectrum to a .txt file
PRO write_theory_0
; get foregrounds
readcol,'/home/cr/cmb_models/foreground_sim09_150.txt', f_ell, f_dl
readcol,'~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt', ksz_ell, ksz_dl
readcol, '/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat', sim_ell, sim_dl

vec = indgen(15000)
wh_f = (where(f_ell eq sim_ell[0]))[0]
wh_ksz = (where(ksz_ell eq sim_ell[0]))[0]

l_tot = sim_ell
dl_tot = sim_dl[vec] + f_dl[vec+wh_f] + ksz_dl[vec+wh_ksz]

openw,lun,/get_lun,'/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_0.txt'
for ii=0, n_elements(dl_tot)-1 do begin
    printf,lun,l_tot[ii], dl_tot[ii]
endfor
free_lun,lun
END

;;;;;;;;;;;;;
; Write the theory spectrum to a .txt file
; CMB + slope
;;;;;;;;;;;;;
PRO write_theory_1
readcol,'/home/cr/cmb_models/foreground_sim09_150.txt', f_ell, f_dl
readcol,'~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt', ksz_ell, ksz_dl
readcol, '/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat', sim_ell, sim_dl

; Make each vector the same length
vec = indgen(15000)
wh_f = (where(f_ell eq sim_ell[0]))[0]
wh_ksz = (where(ksz_ell eq sim_ell[0]))[0]

l_tot = sim_ell[vec]
sim_dl = sim_dl[vec]
f_dl   = f_dl[vec+wh_f]
ksz_dl = ksz_dl[vec+wh_ksz]

; add a slope
sim_dl *= (l_tot/1500.)^0.02

; combine into final spectrum
dl_tot = sim_dl + f_dl + ksz_dl

; Write out file
openw,lun,/get_lun,'/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_1.txt'
for ii=0, n_elements(dl_tot)-1 do begin
    printf,lun,l_tot[ii], dl_tot[ii]
endfor
free_lun,lun
END



;;;;;;;;;;;;;
; Write the theory spectrum to a .txt file
; CMB + slope
;;;;;;;;;;;;;
PRO write_theory_2
readcol,'/home/cr/cmb_models/foreground_sim09_150.txt', f_ell, f_dl
readcol,'~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt', ksz_ell, ksz_dl
readcol, '/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat', sim_ell, sim_dl

; Make each vector the same length
vec = indgen(15000)
wh_f = (where(f_ell eq sim_ell[0]))[0]
wh_ksz = (where(ksz_ell eq sim_ell[0]))[0]

l_tot = sim_ell[vec]
sim_dl = sim_dl[vec]
f_dl   = f_dl[vec+wh_f]
ksz_dl = ksz_dl[vec+wh_ksz]

; add a extra point-source power
dl_extra = 9*(l_tot/3000.)^2.

; combine into final spectrum
dl_tot = sim_dl + f_dl + ksz_dl + dl_extra

stop
; Write out file
openw,lun,/get_lun,'/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_2.txt'
for ii=0, n_elements(dl_tot)-1 do begin
    printf,lun,l_tot[ii], dl_tot[ii]
endfor
free_lun,lun
END



;;;;;;;;;;;;;;;;;;;;;;;
;
; 1) Pipeline test3, spec0
;    specnum is 0, 1 or 2
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe3_dataTF, specnum, save_fig=save_fig
ss = strtrim(string(specnum),2)
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_09p3_spec'+ss+'_kweight.sav'
nsims = 48.

; get pipe2
restore, file
calib = 1.
dl_p = spectrum*(1d12) * (calib)
;cov_p = reform(sample_cov)
cov_p = reform(meas_cov) * 1d24
wf2use = sim_windowfunc
l_wf = findgen(4000-50-1)+50
l = banddef-25.
vec=indgen(47) + 9

; get errors
nbin = n_elements(dl_p)
diags = dblarr(nbin)
;fac = ( 1/sqrt(nsims) + 1/sqrt(100.))
fac = ( 1/nsims + 1/100.)
ccov = cov_p * fac
for ii=0, nbin-1 do begin
    diags[ii]        = sqrt(ccov[ii,ii])
endfor

; condition
nbins_off_diag = 5
ellmin=650
ellmax=3100
ccov = condition_cov_matrix(ccov, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)

ccov = ccov[vec[0]:vec[46], vec[0]:vec[46]]
icov = invert(ccov,/double)

;----------------------------
;Get the theory spectrum ss
readcol, '/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_'+ss+'.txt', l_th, dl_th
istart = where_closest(l_th,l_wf[0])
istop = where_closest(l_th, max(l_wf))
dl_thb = dl_th[istart:istop]
dl_thb = reform(dl_thb # wf2use)

; Get pixel window function
wpix = HEALPIXWINDOW(8192)
wpixb = wpix[l,0]

; divide out the pixel window function from the sims
dl_p /= wpixb^(2.)

; Re-scale
ratio1=(dl_p/dl_thb)
xx = mean(ratio1[vec])
dl_p /= xx

;----------------------------
; Plotting
;----------------------------
;colors
cc_spec = !red
cc_th = !purple


; Calculate chisq
;chisq = total( ( ((dl_p-dl_thb)[vec])/diags[vec] )^2 )
delta = (dl_p-dl_thb)[vec]
chisq = (delta ## (icov ## delta))[0]
print, 'chisq = ', chisq

!p.multi=[0,2,1]
window, 1, xsize=900,ysize=400
plot, l[vec], dl_p[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl,$
  title="Pipeline test3, spec"+ss
oplot, l[vec], dl_p[vec], color=cc_spec, linestyle=0
oplot, l[vec], dl_thb[vec], color=cc_th, linestyle=2
legend,['Calculated PS', 'Input spectrum'],colors=[cc_spec,cc_th],linestyle=[0,2],pos=[1500,2000]

;window, 2
yr=[0.93,1.07]
ratio=(dl_p/dl_thb)
plot, l[vec], ratio[vec], yr=yr,/yst, xtitle='ell',ytitle='Calculate / Input'
errplot, l[vec], ((dl_p-diags)/dl_thb)[vec], ((dl_p+diags)/dl_thb)[vec]
!p.multi=0

if keyword_set(save_fig) then begin
    ff = '/home/kstory/public_html/notebook/spt_lps12/pipe3_spec'+ss+'_0924'
    print, 'Save figure to: ', ff
    err=tvread(/png,/nodialog,filename=ff)
endif
stop
END




;;;;;;;;;;;;;;;;;;;;;;;
;
; 2) Pipeline test3, plot all spectra
;    
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe3_dataTF_2, save_fig=save_fig
svec = [0,1,2]
ratios = dblarr(3,58)

for sindex=0,2 do begin
    specnum = svec[sindex]

ss = strtrim(string(specnum),2)
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_09p3_spec'+ss+'_kweight.sav'
nsims = 48.

; get pipe2
restore, file
calib = 1.
dl_p = spectrum*(1d12) * (calib)
;cov_p = reform(sample_cov)
cov_p = reform(meas_cov) * 1d24
wf2use = sim_windowfunc
l_wf = findgen(4000-50-1)+50
l = banddef-25.
vec=indgen(47) + 9

; get errors
nbin = n_elements(dl_p)
diags = dblarr(nbin)
fac = ( 1/nsims + 1/100.)
ccov = cov_p * fac
for ii=0, nbin-1 do begin
    diags[ii]        = sqrt(ccov[ii,ii])
endfor

; condition
nbins_off_diag = 5
ellmin=650
ellmax=3100
ccov = condition_cov_matrix(ccov, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)

ccov = ccov[vec[0]:vec[46], vec[0]:vec[46]]
icov = invert(ccov,/double)

;----------------------------
;Get the theory spectrum ss
readcol, '/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_'+ss+'.txt', l_th, dl_th
istart = where_closest(l_th,l_wf[0])
istop = where_closest(l_th, max(l_wf))
dl_thb = dl_th[istart:istop]
dl_thb = reform(dl_thb # wf2use)

; Get pixel window function
wpix = HEALPIXWINDOW(8192)
wpixb = wpix[l,0]

; divide out the pixel window function from the sims
dl_p /= wpixb^(2.)

; Re-scale
ratio1=(dl_p/dl_thb)
xx = mean(ratio1[vec])
dl_p /= xx

print, 'Here, k = ', sindex
ratios[sindex,*] = (dl_p/dl_thb)
endfor


;----------------------------
; Plotting
;----------------------------
;colors
cc_spec = !red
cc_th = !purple


;window, 2
yr=[0.93,1.07]
;ratio=(dl_p/dl_thb)
plot, l[vec], ratios[0,vec], yr=yr,/yst, xtitle='ell',ytitle='Calculate / Input'
;errplot, l[vec], ((dl_p-diags)/dl_thb)[vec], ((dl_p+diags)/dl_thb)[vec]
oplot, l[vec], ratios[0,vec], thick=2,color=!red
oplot, l[vec], ratios[1,vec], thick=2,color=!blue
oplot, l[vec], ratios[2,vec], thick=2,color=!darkgreen
legend,['spec0','spec1','spec2'],color=[!red,!blue,!green],linestyle=[0,0,0],thick=[2,2,2]

if keyword_set(save_fig) then begin
    ff = '/home/kstory/public_html/notebook/spt_lps12/pipe3_allspec_0924'
    print, 'Save figure to: ', ff
    err=tvread(/png,/nodialog,filename=ff)
endif
stop
END






;;;;;;;;;;;;;;;;;;;;;;;
;
; 3) Pipeline test3, mcspec0
;    specnum is 1 or 2
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe3_mcspec0, specnum, save_fig=save_fig
ss = strtrim(string(specnum),2)
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_09p3_spec'+ss+'_mcspec0_kweight.sav'
nsims = 48.

; get pipe2
restore, file
calib = 1.
dl_p = spectrum*(1d12) * (calib)
;cov_p = reform(sample_cov)
cov_p = reform(meas_cov) * 1d24
wf2use = sim_windowfunc
l_wf = findgen(4000-50-1)+50
l = banddef-25.
vec=indgen(47) + 9

; get errors
nbin = n_elements(dl_p)
diags = dblarr(nbin)
fac = ( 1/nsims + 1/nsims)
ccov = cov_p * fac
for ii=0, nbin-1 do begin
    diags[ii]        = sqrt(ccov[ii,ii])
endfor

; condition
nbins_off_diag = 5
ellmin=650
ellmax=3100
ccov = condition_cov_matrix(ccov, nbins_off_diag, corr=this_corr, uncondcorr=this_uncorr, banddef=banddef, ellmin=ellmin, ellmax=ellmax, noaverage=noaverage, knowncorr=knowncorr)

ccov = ccov[vec[0]:vec[46], vec[0]:vec[46]]
icov = invert(ccov,/double)

;----------------------------
;Get the theory spectrum ss
readcol, '/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_'+ss+'.txt', l_th, dl_th
istart = where_closest(l_th,l_wf[0])
istop = where_closest(l_th, max(l_wf))
dl_thb = dl_th[istart:istop]
dl_thb = reform(dl_thb # wf2use)

; Get pixel window function
wpix = HEALPIXWINDOW(8192)
wpixb = wpix[l,0]

; divide out the pixel window function from the sims
dl_p /= wpixb^(2.)

; Re-scale
ratio1=(dl_p/dl_thb)
xx = mean(ratio1[vec])
dl_p /= xx

;----------------------------
; Plotting
;----------------------------
;colors
cc_spec = !red
cc_th = !purple


; Calculate chisq
;chisq = total( ( ((dl_p-dl_thb)[vec])/diags[vec] )^2 )
delta = (dl_p-dl_thb)[vec]
chisq = (delta ## (icov ## delta))[0]
print, 'chisq = ', chisq

!p.multi=[0,2,1]
window, 1, xsize=900,ysize=400
plot, l[vec], dl_p[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl,$
  title="Pipeline test3, spec"+ss+"_mcspec0"
oplot, l[vec], dl_p[vec], color=cc_spec, linestyle=0
oplot, l[vec], dl_thb[vec], color=cc_th, linestyle=2
legend,['Calculated PS', 'Input spectrum'],colors=[cc_spec,cc_th],linestyle=[0,2],pos=[1500,2000]

;window, 2
yr=[0.93,1.07]
ratio=(dl_p/dl_thb)
plot, l[vec], ratio[vec], yr=yr,/yst, xtitle='ell',ytitle='Calculate / Input'
errplot, l[vec], ((dl_p-diags)/dl_thb)[vec], ((dl_p+diags)/dl_thb)[vec]
!p.multi=0

if keyword_set(save_fig) then begin
    ff = '/home/kstory/public_html/notebook/spt_lps12/pipe3_spec'+ss+'_mcspec0_0924'
    print, 'Save figure to: ', ff
    err=tvread(/png,/nodialog,filename=ff)
endif
stop
END





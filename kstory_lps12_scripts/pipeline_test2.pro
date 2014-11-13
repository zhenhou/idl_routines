;;;
; NAME: pipeline_test2.pro
; PURPOSE:
;   Analyze sims with different input spectrum.
;
; NOTES:
; 1) Pipeline test2
;;;

FUNCTION get_theory, l_tot=l_tot, dl_tot=dl_tot, shift=shift
theoryspectrumfiles=[ ['/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat',$
                       '/home/cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt',$
                       '/home/cr/cmb_models/foreground_sim09_150.txt']]

; get foregrounds
readcol,'/home/cr/cmb_models/foreground_sim09_150.txt', f_ell, f_dl
readcol,'~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt', ksz_ell, ksz_dl
readcol, '/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat', sim_ell, sim_dl

vec = indgen(15000)
wh_f = (where(f_ell eq sim_ell[0]))[0]
wh_ksz = (where(ksz_ell eq sim_ell[0]))[0]

l_tot = sim_ell
dl_tot = sim_dl[vec] + f_dl[vec+wh_f] + ksz_dl[vec+wh_ksz]
RETURN, 1
END



; Write the theory spectrum to a .txt file
PRO write_theory
theoryspectrumfiles=[ ['/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat',$
                       '/home/cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt',$
                       '/home/cr/cmb_models/foreground_sim09_150.txt']]

; get foregrounds
readcol,'/home/cr/cmb_models/foreground_sim09_150.txt', f_ell, f_dl
readcol,'~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt', ksz_ell, ksz_dl
readcol, '/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat', sim_ell, sim_dl

vec = indgen(15000)
wh_f = (where(f_ell eq sim_ell[0]))[0]
wh_ksz = (where(ksz_ell eq sim_ell[0]))[0]

l_tot = sim_ell
dl_tot = sim_dl[vec] + f_dl[vec+wh_f] + ksz_dl[vec+wh_ksz]

openw,lun,/get_lun,'/data23/kstory/lps12/sims/pipetest2/dl_th_fromCR.txt'
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
openw,lun,/get_lun,'/data23/kstory/lps12/sims/pipetest2/dl_th_fromCR_1.txt'
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
openw,lun,/get_lun,'/data23/kstory/lps12/sims/pipetest2/dl_th_fromCR_2.txt'
for ii=0, n_elements(dl_tot)-1 do begin
    printf,lun,l_tot[ii], dl_tot[ii]
endfor
free_lun,lun
END



;;;;;;;;;;;;;;;;;;;;;;;
;;; Calculate errors on sim spec
;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION get_diags, invkernmat, all_data_spectra
nsims = 100
corspec = reform( invkernmat ## transpose(all_data_spectra) )
diagcov = dblarr(58)
for i=0, 57 do diagcov[i] = variance(corspec[*,i]) / nsims
diags = sqrt(diagcov)*1d24

RETURN, diags
END

;;;;;;;;;;;;;;;;;;;;;;;
;;; Calculate the covariance on sim spec
;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION get_cov, invkernmat, all_data_spectra
nrealizations = 100
;corspec = reform( invkernmat ## transpose(all_data_spectra) )
spectrum_2d=rebin(spectrum, nbands*nspectra, nrealizations)
cov1=(transpose(all_data_spectra-spectrum_2d)##(all_data_spectra-spectrum_2d))/nrealizations/(nrealizations-1)

RETURN, cov1
END



;;;;;;;;;;;;;;;;;;;;;;;
;;; 2) Pipeline test2, spec1
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_spec1
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_09p2_spec1_kweight_calib.sav'

; get pipe2
restore, file
calib = 1.
dl_p1 = spectrum*(1d12)^2. * (calib)
cov_p1 = cov
; get errors
diags = get_diags(invkernmat, all_data_spectra)

; get final combined
restore, edir+'run_09/combined_spectrum_20120828_170101_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

;----------------------------
;Get the theory spectrum 1
tmp = get_theory(l_tot=l_th, dl_tot=dl_th1)
istart = where_closest(l_th,l_wf[0])
istop = where_closest(l_th, max(l_wf))
dl_th1b = dl_th1[istart:istop]
dl_th1b = reform(dl_th1b # wf2use)


vec=indgen(47) + 9

; Calculate chisq
chisq = total( ( ((dl_p1-dl_th1b)[vec])/diags[vec] )^2 )
stop
;----------------------------
; Plotting
;----------------------------
;colors
cc_spec1 = !red
cc_th1 = !purple


!p.multi=[0,2,1]
window, 1, xsize=900,ysize=400
plot, l[vec], dl_p1[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl,$
  title="Pipeline test1"
oplot, l[vec], dl_p1[vec], color=cc_spec1, linestyle=0
oplot, l[vec], dl_th1b[vec], color=cc_th1, linestyle=2
legend,['Calculated PS', 'Input spectrum'],colors=[cc_spec1,cc_th1],linestyle=[0,2],pos=[1500,2000]

;window, 2
yr=[0.93,1.07]
ratio=(dl_p1/dl_th1b)
plot, l[vec], ratio[vec], yr=yr,/yst, xtitle='ell',ytitle='Calculate / Input'
errplot, l[vec], ((dl_p1-diags)/dl_th1b)[vec], ((dl_p1+diags)/dl_th1b)[vec]
!p.multi=0

ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_spec1_0830'
;err=tvread(/png,/nodialog,filename=ff)
stop
END





;;;;;;;;;;;;;;;;;;;;;;;
;
; 2) Pipeline test2, spec1
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_spec2
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_09p2_spec2_kweight_calib.sav'

; get pipe2
restore, file
calib = 1.
dl_p2 = spectrum*(1d12)^2. * (calib)
cov_p2 = cov
; get errors
diags = get_diags(invkernmat, all_data_spectra)


; get final combined
restore, edir+'run_09/combined_spectrum_20120828_170101_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

;----------------------------
;Get the theory spectrum 2
dl_th2 = shift(dl_th1, -10) ; shift by 10:
dl_th2b = dl_th2[istart:istop]
dl_th2b = dl_th2b # wf2use


;----------------------------
; Plotting
;----------------------------
;colors
cc_spec2 = !darkgreen
cc_th2 = !blue

vec=indgen(47) + 9

!p.multi=[0,2,1]
window, 1, xsize=900,ysize=400
plot, l[vec], dl_p2[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl,$
  title="Pipeline test2"
oplot, l[vec], dl_p2[vec], color=cc_spec2,linestyle=0
oplot, l[vec], dl_th2b[vec], color=cc_th2, linestyle=2
legend,['Calculated PS', 'Input spectrum'],colors=[cc_spec2,cc_th2],linestyle=[0,2],pos=[1500,2000]



;window, 2
yr=[0.93,1.07]
plot, l[vec], (dl_p2/dl_th2b)[vec], yr=yr,/yst,$
  xtitle='ell',ytitle='Calculate / Input'
errplot, l[vec], ((dl_p2-diags)/dl_th2b)[vec], ((dl_p2+diags)/dl_th2b)[vec]
!p.multi=0

ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_spec2_0830'
err=tvread(/png,/nodialog,filename=ff)
stop
END






;;;;;;;;;;;;;;;;;;;;;;;
;
; 3) Pipeline test2, all
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_all
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'

; get spec1
restore, edir+'end_ra3h30dec-60_09p2_spec1_kweight_calib.sav'
calib = 1.
dl_p1 = spectrum*(1d12)^2. * (calib)
cov_p1 = cov
diags1 = get_diags(invkernmat, all_data_spectra)

; get spec2
restore, edir+'end_ra3h30dec-60_09p2_spec2_kweight_calib.sav'
calib = 1.
dl_p2 = spectrum*(1d12)^2. * (calib)
cov_p2 = cov
diags2 = get_diags(invkernmat, all_data_spectra)

; get final combined
restore, edir+'run_09/combined_spectrum_20120828_170101_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

;----------------------------
;Get the theory spectrum 1
tmp = get_theory(l_tot=l_th, dl_tot=dl_th1)
istart = where_closest(l_th,l_wf[0])
istop = where_closest(l_th, max(l_wf))
dl_th1b = dl_th1[istart:istop]
dl_th1b = reform(dl_th1b # wf2use)

;----------------------------
;Get the theory spectrum 2
dl_th2 = shift(dl_th1, -10) ; shift by 10:
dl_th2b = dl_th2[istart:istop]
dl_th2b = dl_th2b # wf2use


vec=indgen(47) + 9

;----------------------------
; Calculate chisq
chisq1 = total( ( ((dl_p1-dl_th1b)[vec])/diags1[vec] )^2. )
chisq2 = total( ( ((dl_p2-dl_th2b)[vec])/diags2[vec] )^2. )

stop
;----------------------------
; Plotting
;----------------------------
;colors
cc_spec1 = !red
cc_spec2 = !darkgreen
cc_th1 = !purple
cc_th2 = !blue


!p.multi=[0,2,1]
window, 1, xsize=900,ysize=400
plot, l[vec], dl_p1[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl, title="Pipeline test, all"

; pipe1
oplot, l[vec], dl_p1[vec], color=cc_spec1, linestyle=0
oplot, l[vec], dl_th1b[vec], color=cc_th1, linestyle=2

; pipe2
oplot, l[vec], dl_p2[vec], color=cc_spec2,linestyle=0
oplot, l[vec], dl_th2b[vec], color=cc_th2, linestyle=2

legend,['Calculated 1', 'Input 1', 'Calculated 2', 'Input 2'],colors=[cc_spec1,cc_th1,cc_spec2,cc_th2],$
  linestyle=[0,2,0,2],pos=[1500,2000]


; Pannel 2
plot, l[vec], (dl_th1b/dl_th2b)[vec], $
  yr=[0.93,1.07],/yst,xtitle='ell',ytitle='Calculate / Input',title='input1 - input2'
oplot, l[vec], (dl_th1b/dl_th2b)[vec], thick=2

; oplot, l[vec], (dl_p1/dl_th1b)[vec],thick=2, color=!darkgreen
; oplot, l[vec], (dl_p2/dl_th2b)[vec],thick=2, color=!blue
!p.multi=0

ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_all_0830'
;err=tvread(/png,/nodialog,filename=ff)

;------------------
; Ratios
window, 2
plot, l[vec], (dl_th1b/dl_th2b)[vec], $
  yr=[0.93,1.07],/yst,xtitle='ell',ytitle='Calculate / Input'
oplot, l[vec], (dl_th1b/dl_th2b)[vec], thick=2

oplot, l[vec], (dl_p1/dl_th1b)[vec],thick=2, color=!darkgreen
oplot, l[vec], (dl_p2/dl_th2b)[vec],thick=2, color=!blue
legend,['th1/th2', 'spec1/th1', 'spec2/th2'],colors=[!black,!darkgreen,!blue], $
  linestyle=[0,0,0]
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_ratios_0830'
err=tvread(/png,/nodialog,filename=ff)

;------------------
; cross:
window, 3
plot, l[vec], (dl_th1b/dl_th2b)[vec], $
  yr=[0.93,1.07],/yst,xtitle='ell',ytitle='Calculate / Input'
oplot, l[vec], (dl_th1b/dl_th2b)[vec], thick=2, color=!red
oplot, l[vec], (dl_p1/dl_th2b)[vec],thick=2, color=!orange

oplot, l[vec], (dl_th2b/dl_th1b)[vec], thick=2, color=!blue, linestyle=2
oplot, l[vec], (dl_p2/dl_th1b)[vec],thick=2, color=!skyblue, linestyle=2

legend,['th1/th2', 'spec1/th2', 'th2/th1','spec2/th1'],colors=[!red,!orange,!blue,!skyblue], $
  linestyle=[0,0,2,2]
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_cross_ratios_0830'
err=tvread(/png,/nodialog,filename=ff)


stop
END






;;;;;;;;;;;;;;;;;;;;;;;
;
; 4) spec2, using sims=spec1
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_spec2_mcspec1
f = lps12_fieldstruct()
file = '/home/kstory/lps12/end2end/end_ra3h30dec-60_09p2_spec2_mcspec1_kweight.sav'
nsims = 100.

; get pipe2
restore, file
calib = 1.
dl_p2 = spectrum*(1d12) * (calib)
cov_p2 = reform(sample_cov)
wf2use = sim_windowfunc

; get errors
nbin = n_elements(dl_p2)
diags = dblarr(nbin)
for ii=0, nbin-1 do begin
    diags[ii]        = sqrt(cov_p2[ii,ii])
endfor
diags = diags*(1d12) / sqrt(nsims)

;----------------------------
;Get the theory spectrum 2
readcol, spectrumfiles[0], format='d,d', l_th1, dl_th1
dl_th2 = shift(dl_th1, -10) ; shift by 10:

l_wf = findgen(4000-50-1)+50
istart = where_closest(l_th1,l_wf[0])
istop = where_closest(l_th1, max(l_wf))

l_th2   = l_th1[istart:istop]
dl_th2  = dl_th2[istart:istop]

dl_th2b = dl_th2 # wf2use

dl_th2b_5 = ((shift(dl_th1, -5))[istart:istop])# wf2use
dl_th2b_7 = ((shift(dl_th1, -7))[istart:istop])# wf2use
dl_th2b_9 = ((shift(dl_th1, -9))[istart:istop])# wf2use
dl_th2b_11 = ((shift(dl_th1, -11))[istart:istop])# wf2use
dl_th2b_13 = ((shift(dl_th1, -13))[istart:istop])# wf2use
dl_th2b_15 = ((shift(dl_th1, -15))[istart:istop])# wf2use


;----------------------------
; Plotting
;----------------------------
;colors
cc_spec2 = !darkgreen
cc_th2 = !blue

l = banddef-25.
vec=indgen(47) + 9

; Calculate the chisq
chisq2 = total( ( ((dl_p2-dl_th2b)[vec])/diags[vec] )^2. )
print, 'chisq = ', chisq2

!p.multi=[0,2,1]
window, 1, xsize=900,ysize=400
plot, l[vec], dl_p2[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl,$
  title="Pipeline test2"
oplot, l[vec], dl_p2[vec], color=cc_spec2,linestyle=0
oplot, l[vec], dl_th2b[vec], color=cc_th2, linestyle=2
errplot, l[vec], (dl_p2-diags)[vec], (dl_p2+diags)[vec]
legend,['Calculated PS', 'Input spectrum'],colors=[cc_spec2,cc_th2],linestyle=[0,2],pos=[1500,2000]

;window, 2
yr=[0.98,1.02]
plot, l[vec], (dl_p2/dl_th2b)[vec], yr=yr,/yst,$
  xtitle='ell',ytitle='Calculate / Input'
errplot, l[vec], ((dl_p2-diags)/dl_th2b)[vec], ((dl_p2+diags)/dl_th2b)[vec]
!p.multi=0

ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_spec2_mcspec1_0914'
err=tvread(/png,/nodialog,filename=ff)


window, 2
yr=[0.96,1.04]
plot, l[vec], (dl_p2/dl_th2b)[vec], yr=yr,/yst,$
  xtitle='ell',ytitle='Calculate / Input'

oplot, l[vec], (dl_p2/dl_th2b_5)[vec], color=!red
oplot, l[vec], (dl_p2/dl_th2b_7)[vec], color=!blue
oplot, l[vec], (dl_p2/dl_th2b_9)[vec], color=!green
oplot, l[vec], (dl_p2/dl_th2b_11)[vec], color=!skyblue
oplot, l[vec], (dl_p2/dl_th2b_13)[vec], color=!purple
oplot, l[vec], (dl_p2/dl_th2b_15)[vec], color=!orange

stop
END





;;;;;;;;;;;;;;;;;;;;;;;
;
; 5) What ell can we resolve?
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO chisq_spec2_mcspec1
f = lps12_fieldstruct()
file = '/home/kstory/lps12/end2end/end_ra3h30dec-60_09p2_spec2_mcspec1_kweight.sav'
nsims = 100.

; get pipe2
restore, file
calib = 1.
dl_p2 = spectrum*(1d12) * (calib)
cov_p2 = reform(sample_cov)
wf2use = sim_windowfunc

; get errors
nbin = n_elements(dl_p2)
diags = dblarr(nbin)
for ii=0, nbin-1 do begin
    diags[ii]        = sqrt(cov_p2[ii,ii])
endfor
diags = diags*(1d12) / sqrt(nsims)

;----------------------------
;Get the theory spectrum 2
readcol, spectrumfiles[0], format='d,d', l_th1, dl_th1
dl_th2 = shift(dl_th1, -10) ; shift by 10:

l_wf = findgen(4000-50-1)+50
istart = where_closest(l_th1,l_wf[0])
istop = where_closest(l_th1, max(l_wf))

l_th2   = l_th1[istart:istop]
dl_th2  = dl_th2[istart:istop]

dl_th2b = dl_th2 # wf2use

; dl_th2b_5 = ((shift(dl_th1, -5))[istart:istop])# wf2use
; dl_th2b_7 = ((shift(dl_th1, -7))[istart:istop])# wf2use
; dl_th2b_9 = ((shift(dl_th1, -9))[istart:istop])# wf2use
; dl_th2b_10 = ((shift(dl_th1, -9))[istart:istop])# wf2use
; dl_th2b_11 = ((shift(dl_th1, -11))[istart:istop])# wf2use
; dl_th2b_13 = ((shift(dl_th1, -13))[istart:istop])# wf2use
; dl_th2b_15 = ((shift(dl_th1, -15))[istart:istop])# wf2use

; Calculate chisq
l = banddef-25.
vec=indgen(47) + 9

chisq2 = total( ( ((dl_p2-dl_th2b)[vec])/diags[vec] )^2. )
print, 'chisq = ', chisq2


sh=indgen(11)+5 & nsh = n_elements(sh)
chisq = fltarr(nsh)
nsigma = fltarr(nsh)
dof = n_elements(vec)

cc = kscolor_array()

window, 0
dl_th_tmp = ((shift(dl_th1, -1*sh[0]))[istart:istop])# wf2use
plot, l[vec], (dl_th_tmp - dl_p2)[vec], xtitle='ell',ytitle='dl_shift - dl_spec1, [uK^2]',title='Absolute difference in spectra'

for ii=0, nsh-1 do begin
    dl_th_tmp = ((shift(dl_th1, -1*sh[ii]))[istart:istop])# wf2use
    chisq[ii] = total( ( ((dl_p2-dl_th_tmp)[vec])/diags[vec] )^2. )
    nsigma[ii] = mpchitest(chisq[ii],dof,/sigma)

    oplot, l[vec], (dl_th_tmp - dl_p2)[vec], color=cc[ii]
endfor
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_diff_0918'
err=tvread(/png,/nodialog,filename=ff)

window, 1
plot, sh-10, nsigma, xtitle='shift in ell', ytitle='sigma', title='Change in chisq due to shift in ell'
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_nsigma_0918'
err=tvread(/png,/nodialog,filename=ff)

window, 2
plot, sh-10, chisq, xtitle='shift in ell', ytitle='chisq (47 dof)', title='Change in chisq due to shift in ell'
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_0918'
err=tvread(/png,/nodialog,filename=ff)

stop
END





;;;;;;;;;;;;;;;;;;;;;;;
;
; 6) Convolve the theory spectrum with a gaussian
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO chisq_spec2_mcspec1_gauss
f = lps12_fieldstruct()
file = '/home/kstory/lps12/end2end/end_ra3h30dec-60_09p2_spec2_mcspec1_kweight.sav'
nsims = 100.

; get pipe2
restore, file
calib = 1.
dl_p2 = spectrum*(1d12) * (calib)
cov_p2 = reform(sample_cov)
wf2use = sim_windowfunc

; get errors
nbin = n_elements(dl_p2)
diags = dblarr(nbin)
for ii=0, nbin-1 do begin
    diags[ii]        = sqrt(cov_p2[ii,ii])
endfor
diags = diags*(1d12) / sqrt(nsims)

;----------------------------
;Get the theory spectrum 2
readcol, spectrumfiles[0], format='d,d', l_th1, dl_th1
dl_th2 = shift(dl_th1, -10) ; shift by 10:

l_wf = findgen(4000-50-1)+50
istart = where_closest(l_th1,l_wf[0])
istop = where_closest(l_th1, max(l_wf))

l_th2   = l_th1[istart:istop]
dl_th2  = dl_th2[istart:istop]

dl_th2b = dl_th2 # wf2use


; Plotting stuff
cc = kscolor_array()
l = banddef-25.
vec=indgen(47) + 9

; Gaussian smoothing kernels
;gauss_smooth(cl_th,l_fwhm,lbinsize)
l_fwhm = indgen(10)*5 + 1
nsh = n_elements(l_fwhm)
chisq = fltarr(nsh)
nsigma = fltarr(nsh)
dof = n_elements(vec)

window, 0
plot, l_th2, dl_th2, xtitle='ell',ytitle='dl[uK^2]',/yl

window, 1
dl_th_tmp = gauss_smooth(dl_th2,l_fwhm[nsh-1],1)
dl_th_tmp = dl_th_tmp# wf2use
plot, l[vec], (dl_th_tmp - dl_p2)[vec], xtitle='ell',ytitle='dl_smooth - dl_spec1, [uK^2]',title='Absolute difference in spectra'

for ii=0, nsh-1 do begin
    dl0 = gauss_smooth(dl_th2,l_fwhm[ii],1)
    dl1 = dl0# wf2use
    chisq[ii] = total( ( ((dl_p2-dl1)[vec])/diags[vec] )^2. )
    nsigma[ii] = mpchitest(chisq[ii],dof,/sigma)

    wset, 0
    oplot, l_th2, dl0, color=cc[ii]

    wset, 1
    oplot, l[vec], (dl1 - dl_p2)[vec],color=cc[ii]

;    stop
endfor
;ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_diff_gauss_0919'
;err=tvread(/png,/nodialog,filename=ff)

window, 2
plot, l_fwhm, nsigma, xtitle='shift in ell', ytitle='sigma', title='Change in chisq due to Gaussin smoothing'
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_nsigma_gauss_0919'
err=tvread(/png,/nodialog,filename=ff)

window, 3
plot, l_fwhm, chisq, xtitle='shift in ell', ytitle='chisq (47 dof)', title='Change in chisq due to Gaussin smoothing'
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_gauss_0919'
err=tvread(/png,/nodialog,filename=ff)

stop
END





;;;;;;;;;;;;;;;;;;;;;;;
;
; 7) Stretching the theory spectrum
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO chisq_spec2_mcspec1_stretch
f = lps12_fieldstruct()
file = '/home/kstory/lps12/end2end/end_ra3h30dec-60_09p2_spec2_mcspec1_kweight.sav'
nsims = 100.

; get pipe2
restore, file
calib = 1.
dl_p2 = spectrum*(1d12) * (calib)
cov_p2 = reform(sample_cov)
wf2use = sim_windowfunc

; get errors
nbin = n_elements(dl_p2)
diags = dblarr(nbin)
for ii=0, nbin-1 do begin
    diags[ii]        = sqrt(cov_p2[ii,ii])
endfor
diags = diags*(1d12) / sqrt(nsims)

;----------------------------
;Get the theory spectrum 2
readcol, spectrumfiles[0], format='d,d', l_th1, dl_th1
dl_th2 = shift(dl_th1, -10) ; shift by 10:

l_wf = findgen(4000-50-1)+50
istart = where_closest(l_th1,l_wf[0])
istop = where_closest(l_th1, max(l_wf))

l_th2   = l_th1[istart:istop]
dl_th2  = dl_th2[istart:istop]

dl_th2b = dl_th2 # wf2use


; Plotting stuff
cc = kscolor_array()
l = banddef-25.
vec=indgen(47) + 9

; Gaussian smoothing kernels
;gauss_smooth(cl_th,l_fwhm,lbinsize)
ss = findgen(15)*0.0002 + 1.0002
nss = n_elements(ss)
chisq = fltarr(nss)
nsigma = fltarr(nss)
dof = n_elements(vec)

; l-vec for interpolating 
l_interp = l_th2-500

window, 0
plot, l_th2, dl_th2, xtitle='ell',ytitle='dl[uK^2]',/yl

window, 1
dl0 = interpol(dl_th2, l_interp, l_interp/ss[nss-1])
dl1 = dl0# wf2use
plot, l[vec], (dl1 - dl_p2)[vec], xtitle='ell',ytitle='dl_smooth - dl_spec1, [uK^2]',title='Absolute difference in spectra'

for ii=0, nss-1 do begin
    dl0 = interpol(dl_th2, l_interp, l_interp/ss[ii])
    dl1 = dl0# wf2use
    chisq[ii] = total( ( ((dl_p2-dl1)[vec])/diags[vec] )^2. )
    nsigma[ii] = mpchitest(chisq[ii],dof,/sigma)

;     wset, 0
;     oplot, l_th2, dl0, color=cc[ii]

    wset, 1
    oplot, l[vec], (dl1 - dl_p2)[vec],color=cc[ii]

;    stop
endfor
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_diff_stretch_0919'
;ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_spec_stretch_0919'
err=tvread(/png,/nodialog,filename=ff)

window, 2
plot, ss, nsigma, psym=-2,xtitle='stretch-factor', ytitle='sigma', title='Change in chisq due to Stretching; dl/(factor)'
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_nsigma_stretch_0919'
err=tvread(/png,/nodialog,filename=ff)

window, 3
plot, ss, chisq, xtitle='stretch-factor', ytitle='chisq (47 dof)', title='Change in chisq due to Stretching; dl/(factor)'
ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_chisq_stretch_0919'
err=tvread(/png,/nodialog,filename=ff)

stop
END






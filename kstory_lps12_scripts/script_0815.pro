;;;
; NAME: script_0815
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Pipeline test2
; 3) testing dl_all plots
;;;

FUNCTION get_theory
;;; get dl_th
; readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_uK2
; dl_uK2_2 = dl_uK2[50:3998]
; l_vec_2  = l_vec[50:3998]
; dl_th_2 = dl_uK2_2 # wf2use
theoryspectrumfiles=[ ['/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat',$
                       '/home/cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt',$
                       '/home/cr/cmb_models/foreground_sim09_150.txt']]

; get foregrounds
readcol,'/home/cr/cmb_models/foreground_sim09_150.txt', f_ell, f_dl
readcol,'~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt', ksz_ell, ksz_dl
readcol, '/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat', sim_ell, sim_dl

wh = indgen(15000)
dl_tot = sim_dl[wh] + f_dl[wh] + ksz_dl[wh]

RETURN, dl_tot
END


;;;;;;;;;;;;;;;;;;;;;;;
;;; 2) Pipeline test2
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_08p2_kweight.sav'

; get pipe2
restore, file
calib = 1.
;calib = 0.825
;calib = 0.76
;calib = 0.85
dl_p2 = spectrum*(1d12)^2. * (calib)
cov_p2 = cov

; get final combined
restore, edir+'run_08/combined_spectrum_20120717_174249_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

; readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_uK2
; dl_uK2_2 = dl_uK2[50:3998]
; l_vec_2  = l_vec[50:3998]
; dl_th_2 = dl_uK2_2 # wf2use

;----------------------------
;Get the theory spectrum
dl_th = get_theory()
x_th = indgen(3001) + 500
dl_th_plot = dl_th[x_th]

dl_th_2 = dl_th[50:3998]
dl_th_2 = reform(dl_th_2 # wf2use)

; Plotting
vec=indgen(47) + 9

window, 1
plot, l[vec], dl_p2[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl
oplot, l[vec], dl_p2[vec], color=!red, linestyle=3
;oplot, l[vec], dl_all[vec]

oplot, l[vec], dl_th_2[vec], color=!green, linestyle=0
oplot, x_th, dl_th_plot, color=!blue, linestyle=2

window, 2
plot, l[vec], (dl_th_2/dl_p2)[vec]
stop
END




; freq=150
; pi=acos(-1.d0)
; n=4096 ;x2 = 13824 factors into 2s and 3s
; resolution=(0.5/60d0)/180d0*pi
; waka=pi/(2d0*n*resolution)
; df=float(waka)/pi
; dtheta=1./df/2./n
; spectrum05=complexarr(n*2,n*2) ; will be combined for all years
; gar=intarr(2)
; calc=dblarr(25000)
; calccmb=dblarr(25000)
; ll=dindgen(25000)+2
; dlfac = ll*(ll+1.)/(2*!pi)
; fwhm2sig = 1./sqrt(8.*alog(2.))

; scmb=11 
; sps=23
; sps2=31
; calc = dlfac / ((3000.*3001.)/(2.*!pi))
;                                 ; I've normalized with respect to
;                                 ; D3000=1
; calcell = calc /( ll/3000.) ; likewise D3000=1

; facs = estimate_pointsources_sim09()
; matrix = [[facs.ps[0],facs.ps[3],facs.ps[4]],$
;           [facs.ps[3],facs.ps[1],facs.ps[5]],$
;           [facs.ps[4],facs.ps[5],facs.ps[2]]]
; ellmatrix = [[facs.ellclus[0],facs.ellclus[3],facs.ellclus[4]],$
;              [facs.ellclus[3],facs.ellclus[1],facs.ellclus[5]],$
;              [facs.ellclus[4],facs.ellclus[5],facs.ellclus[2]]]
; evecs=matrix
; trired,evecs,evals,e
; triql,evals,e,evecs

; ;at freq=150
; ii0=1

; amps = sqrt(evals)
; effs=fltarr(3)
; effs[0] = amps[0]*evecs[ii0,0]
; effs[1] = amps[1]*evecs[ii0,1]
; effs[2] = amps[2]*evecs[ii0,2]

; ;check if planned N makes sense
; jj=abs(evals)
; ind = where(jj gt 1e-4*max(jj),nbig,complement=cnd)
; if nbig ne 2 then stop
; if nbig lt 3 then effs[cnd]=0

; evecs=ellmatrix
; trired,evecs,evals,e
; triql,evals,e,evecs

; jj=abs(evals)
; ind = where(jj gt 1e-4*max(jj),nbig,complement=cnd)
; if nbig ne 1 then stop
; ampsclus = sqrt(evals)

; effsclus=fltarr(3)
; effsclus[0] = ampsclus[0]*evecs[ii0,0]
; effsclus[1] = ampsclus[1]*evecs[ii0,1]
; effsclus[2] = ampsclus[2]*evecs[ii0,2]
; if nbig lt 3 then effsclus[cnd]=0
; ;stop

; calccmb= spt_th_cl_run3_2009(ll+10,freq,/onlycmb)

; ;also read in kSZ
; a = read_ascii('~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt')
; ksz=ll*0.0
; ksz[0:16300]=a.field1[1,2:16302]
; ;smoothly go to zero
; ksz[16301:24999] = ksz[16300] - (ksz[16300]/8699.)*(1+findgen(8699))

; ; get foregrounds
; readcol,'/home/cr/cmb_models/foreground_sim09_150.txt', f_ell, foreground
; readcol,'~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt', ksz_ell, ksz_vec

; calccmb += ksz
; ;calccmb += foreground

; RETURN, calccmb
; END

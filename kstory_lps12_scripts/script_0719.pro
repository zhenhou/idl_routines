;;;
; NAME: script_0719
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Test SV averaging
;
; MODIFICATION HISTORY:
;  07/19/2012: (KTS) Created
;;;

;;; Check the transfer function levels for the paper
pro tf

edir = '/data23/kstory/lps12/end2end/'
f = lps12_fieldstruct()

low  = 1400
high = 2200

for i=7,19 do begin
    restore, edir+'end_'+f[i].name+'_08_kweight.sav'
    tf = transfer

    wh = where(ellkern ge 650) & ss = wh[0]
    wh = where(ellkern ge low and ellkern le high)


    plot, ellkern, transfer, xr=[0,3500], yr=[0,1]
    oplot, [low,low],[0,1], color=!red
    oplot, [high,high],[0,1], color=!red
    oplot, ellkern[wh], transfer[wh], color=!green

    print, i, ', field '+f[i].name
    print, ellkern[ss], transfer[ss]
    print, mean(transfer[wh])

endfor

stop
END


;;; Find where cov_sv becomes smaller than cov_data
pro error
restore, '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
plot, l, d_sv*1d12, /ylog, yr=[0.1, 1000]
oplot, l, d_data*1d12, color=!red

window, 2
plot, l[40:52], (d_sv[40:52])*1d12
oplot, l[40:52], (d_data[40:52])*1d12, color=!red
end


;;; Make PS plot
pro ps_plot1
restore, '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
; get the input theory spectrum
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)

xr = [500,3.5e3]
yr = [40,4e3]
xtitle= '!12l!X!N'
ell = '!12l'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
ytitle= ell+'('+ell+'+1)C'+'!D!12l!X!N'+textoidl('/2\pi [\muK^2]')
chars=1.8
whplot = indgen(45) + 10

dl2use = dl_all
diag2use=diag_nobeam
l2use = l

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SPT data points only
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
window,0
plot,l2use[whplot],dl2use[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars
errplot,l2use[whplot],(dl2use-diag2use)[whplot]*1d12,(dl2use+diag2use)[whplot]*1d12
xyouts,2000.,1100,'SPT,!C Full 2500 deg'+textoidl('^2'),color=color_data,chars=3.,font=1

; SPT plus other data sets
; get ACBAR                             
color_acbar=!darkgreen
readcol,'acbar_reichardt.txt',band_num_acbar,dl_acbar,ddl_acbar,log_offset_acbar
l_acbar = [470., 608., 694., 763., 823., 884., 944., 1003., 1062., 1122., 1182., 1242., 1301., 1361., 1421., 1481., 1541., 1618., 1713., 1814., 1898., 2020., 2194., 2391., 2646.]

; get WMAP7 binned power spectrum                                                                                   
readcol,'wmap_binned_tt_spectrum_7yr_v4p1.txt',l_wmap,llo,lhi,dl_wmap,ddl_wmap,ddl_noise,ddl_sv
;wh_wmap = where(l_wmap le 800.)
wh_wmap = where(l_wmap le 8000. and l_wmap gt 2)


window, 1
plot,l2use[whplot],dl2use[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars
;ACBAR
errplot,l2use[whplot],(dl2use-diag2use)[whplot]*1d12,(dl2use+diag2use)[whplot]*1d12
errplot,l_acbar,(dl_acbar-ddl_acbar),(dl_acbar+ddl_acbar), color=color_acbar

;WMAP
errplot,l_wmap[wh_wmap],(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap+ddl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color
errplot,l_wmap[wh_wmap],(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color,width=0.008;centers
oplot,[0,0],[1d-20,1d20],color=!black,thick=1
oplot,[-1e9,1e9],[1,1]*40.,color=!black,thick=1

stop
end


;;;;;;;;;;;;;;;;;;
; Procedure to make over-plot
PRO reg_plot
pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
xmargin=[8.5,2]
ymargin=[4,2]

readcol,pdir+'wmap7cmb_bestfit_dl_k11.txt',lth_cmb,dlth_cmb


;size=3.9
size=6.
;chars=0.9
chars=1.5
lchars=chars*1.6
lfont=1
filename = 'figs/dl_all.eps'
openplotps,filen=filename,xs=size*1.8,ys=0.8*size,/eps
!p.multi = [0,2,1]
;window, 1, xsize=1000, ysize=400
setcolors2,/sys
thick=4

xtitle='!12l!X!N/1000'
xtitle='!12l!X!N'
ytitle='!12l(l+1)C!D!12l!3!N/2!7p!X'+textoidl(' [\muK^2]')
sf=1.
xr=[200,3400]/sf
xxx = 200. & xr = [650-xxx,3000+xxx]
yr=[30,4e3]
lines_neg = 2
nsmooth=5

; Make the plot
plot,[0],[0],/nodata,xr=xr,yr=yr,/xst,/yst,xtitle=xtitle,ytitle=ytitle, $
  chars=chars,/yl,xtickint=1000/sf,xmargin=xmargin,ymargin=ymargin

; Get the data
restore,file

dl_all *= (1d12)
diag_nobeam *= (1d12)


istart=9
istop=55
dl_all = dl_all[istart:istop]
diag_nobeam = diag_nobeam[istart:istop]
l = l[istart:istop]

color_this_work = !black
;color_this_work = !blue
oplot,l/sf,dl_all,ps=3,color=color_this_work
errplot,l/sf,(dl_all-diag_nobeam),(dl_all+diag_nobeam),thick=thick,color=color_this_work
;oplot,lth_cmb/sf,dlth_cmb,lines=2

xtxt=2060
ytxt=2000*0.9
dytxt = 0.63
bigspt = 1.05
xyouts,xtxt,ytxt,'SPT,',chars=lchars*bigspt,font=lfont,color=color_this_work
xyouts,xtxt,ytxt*dytxt,'Full Survey',chars=lchars*bigspt,font=lfont,color=color_this_work

;;;;;;;;;;;;;;;
; RIGHT PANEL
;;;;;;;;;;;;;;;
color_quad = !darkgreen
color_acbar = !blue
color_act = !darkorange
color_shiro = !orange


; load quad
readcol,pdir+'quad_lowell.txt',l_quad1,dl_quad1,delta_dl_quad1,dl_lcdm_quad1
readcol,pdir+'quad_highell.txt',l_quad2,dl_quad2,delta_dl_quad2,dl_lcdm_quad2
l_quad = [l_quad1, l_quad2]
dl_quad = [dl_quad1, dl_quad2]
ddl_quad = [delta_dl_quad1, delta_dl_quad2]

; load acbar
readcol,pdir+'cl_dc1_v3.dat',iband,dl_acbar,delta_dl_acbar,log_normal_offset
ledge=[[350,550],[550,650],[650,730],[730,790],[790,850],[850,910],[910,970],[970,1030],[1030,1090],[1090,1150],[1150,1210],[1210,1270],[1270,1330],[1330,1390],[1390,1450],[1450,1510],[1510,1570],[1570,1650],[1650,1750],[1750,1850],[1850,1950],[1950,2100],[2100,2300],[2300,2500],[2500,3000]]
llo=reform(ledge[0,*])
lhi=reform(ledge[1,*])
l_acbar=0.5*(llo+lhi)
ddl_acbar = delta_dl_acbar

; load act (das et al)
readcol,pdir+'act_150_use.txt',l_act,dl_act,ddl_act

; load spt shirokoff et al
readcol,pdir+'shiro.txt',llo_shiro,lhi_shiro,l_shiro,dl150,ddl150,dl150220,ddl150220,dl220,ddl220
dl_shiro = dl150
ddl_shiro = ddl150

xtitle='!12l!X!N/1000'
xtitle='!12l!X!N'
ytitle='!12l(l+1)C!D!12l!3!N/2!7p!X'+textoidl(' [\muK^2]')
sf=1.
xr=[200,3400]/sf
xxx = 200. & xr = [650-xxx,3000+xxx]
yr=[30,4e3]
lines_neg = 2
nsmooth=5

; Plot #2
plot,[0],[0],/nodata,xr=xr,yr=yr,/xst,/yst,xtitle=xtitle,ytitle=ytitle, $
  chars=chars,/yl,xtickint=1000/sf,xmargin=xmargin,ymargin=ymargin

restore,file

dl_all *= (1d12)
diag_nobeam *= (1d12)

istart=9
istop=55
dl_all = dl_all[istart:istop]
diag_nobeam = diag_nobeam[istart:istop]
l = l[istart:istop]

errplot,l_acbar,dl_acbar-ddl_acbar,dl_acbar+ddl_acbar,color=color_acbar,thick=thick
errplot,l_quad,dl_quad-ddl_quad,dl_quad+ddl_quad,color=color_quad,thick=thick
errplot,l_act,dl_act-ddl_act,dl_act+ddl_act,color=color_act,thick=thick
errplot,l_shiro,dl_shiro-ddl_shiro,dl_shiro+ddl_shiro,color=color_shiro,thick=thick

xyouts,xtxt,ytxt*dytxt^0.,'ACBAR',chars=lchars,font=lfont,color=color_acbar
xyouts,xtxt,ytxt*dytxt^1.,'QUaD',chars=lchars,font=lfont,color=color_quad
xyouts,xtxt,ytxt*dytxt^2.,'ACT',chars=lchars,font=lfont,color=color_act
xyouts,xtxt,ytxt*dytxt^3.,'SPT, S10',chars=lchars,font=lfont,color=color_shiro

closeps
spawn,'convert '+filename+' figs/dl_all.pdf'

stop
END


;;;;;;;;;;;;;;;;;;
; Procedure to make over-plot
PRO plot_wmap

pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'

readcol,pdir+'wmap7cmb_bestfit_dl_k11.txt',lth_cmb,dlth_cmb

stop
END



;;;;;;;;;;;;;;;;;;
; Make best-fit spectra for SPT+WMAP
;;;;;;;;;;;;;;;;;;
pro make_bestfit_s12_spectra, filename, data=data, $
                              dl_total=dl_total, dl_cmb=dl_cmb, l_cmb=l_cmb, $
                              noplot=noplot, stopit=stopit, ratio=ratio, err_ratio=err_ratio, $
                              l_ratio=l, doeps=doeps, alens=alens, $
                              dosave=dosave, $
                              in_l_max_scalar=in_l_max_scalar, $
                              in_k_eta_max_scalar=in_k_eta_max_scalar

pdir = '/home/kstory/lps12/scripts/plotting/'
if n_elements(filename) eq 0 then filename=pdir+'cons_chain_baseline_k11.sav'
if n_elements(data) eq 0 then $
  data = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
  ;data='/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'

;a = read_ascii(filename)
;y = a.field01
restore,filename
y = echain

;no_lensing = STRMATCH(filename, '*unlensed*')
;print,no_lensing

lnl = reform(y[1,*])
min_lnl = min(lnl)
wh_min = where(lnl eq min(lnl), nmin)
if nmin gt 1 then begin
    print,'yo: there were more than one points in the chain with the minimum log likelihood'
endif
wh_min = wh_min[0]

ombh2 = y[2,wh_min]
omch2 = y[3,wh_min]
yhe = y[10,wh_min]
neff = y[9,wh_min]
tau = y[5,wh_min]


omch2_orig = omch2


;k = keyword_set(alens) ? 1 : 0
k=1
hubble = y[42+k,wh_min]
as = exp(y[14+k,wh_min])/1d10
ns = y[11+k,wh_min]
nrun = y[13+k,wh_min]
nt = y[12+k,wh_min]
r = y[15+k,wh_min]

in_l_max_scalar = 13000.

;goto, yoyo


r = 0. & alens=1. & omnuh2=0. & omch2 = omch2_orig - omnuh2
params_to_camb, no_lensing=no_lensing, $
  ombh2=ombh2, omch2=omch2, $
  hubble=hubble, yhe=yhe, neff=neff, $
  as=as, ns=ns, nrun=nrun, nt=nt, r=r, $
  tau=tau, $
  in_l_max_scalar=in_l_max_scalar, $
  in_k_eta_max_scalar=in_k_eta_max_scalar, $
  l=l, dl=tt,phi=phi, $
  ee=ee, bb=bb, te=te, $
  alens=1., accuracy_level=2., omnuh2=omnuh2
tab = [[l],[tt],[te],[ee],[bb]]
writetab,transpose(tab),pdir+'wmap_s12_bestfit_cmbonly_dl_lcdm.txt'

r = 0. & alens=0. & omnuh2=0. & omch2 = omch2_orig - omnuh2
params_to_camb, no_lensing=no_lensing, $
  ombh2=ombh2, omch2=omch2, $
  hubble=hubble, yhe=yhe, neff=neff, $
  as=as, ns=ns, nrun=nrun, nt=nt, r=r, $
  tau=tau, $
  in_l_max_scalar=in_l_max_scalar, $
  in_k_eta_max_scalar=in_k_eta_max_scalar, $
  l=l, dl=tt,phi=phi, $
  ee=ee, bb=bb, te=te, $
  alens=0., accuracy_level=2., omnuh2=omnuh2
tab = [[l],[tt],[te],[ee],[bb]]
writetab,transpose(tab),pdir+'wmap_s12_bestfit_cmbonly_dl_lcdm_unlensed.txt'


stop

r = 0. & alens=1. & omnuh2=0. & omch2 = omch2_orig - omnuh2
params_to_camb, no_lensing=no_lensing, $
  ombh2=ombh2, omch2=omch2, $
  hubble=hubble, yhe=yhe, neff=neff, $
  as=as, ns=ns, nrun=nrun, nt=nt, r=r, $
  tau=tau, $
  in_l_max_scalar=in_l_max_scalar, $
  in_k_eta_max_scalar=in_k_eta_max_scalar, $
  l=l, dl=tt,phi=phi, $
  ee=ee, bb=bb, te=te, $
  alens=1., accuracy_level=2., omnuh2=omnuh2
tab = [[l],[tt],[te],[ee],[bb]]
writetab,transpose(tab),pdir+'wmap_s12_bestfit_r0p00_summnu0p0eV.txt'


r = 0.02 & alens=1. & omnuh2=0. & omch2 = omch2_orig - omnuh2
params_to_camb, no_lensing=no_lensing, $
  ombh2=ombh2, omch2=omch2, $
  hubble=hubble, yhe=yhe, neff=neff, $
  as=as, ns=ns, nrun=nrun, nt=nt, r=r, $
  tau=tau, $
  in_l_max_scalar=in_l_max_scalar, $
  in_k_eta_max_scalar=in_k_eta_max_scalar, $
  l=l, dl=tt,phi=phi, $
  ee=ee, bb=bb, te=te, $
  alens=1., accuracy_level=2., omnuh2=omnuh2
tab = [[l],[tt],[te],[ee],[bb]]
writetab,transpose(tab),pdir+'wmap_s12_bestfit_r0p02_summnu0p0eV.txt'


r = 0. & alens=1. & omnuh2=0.1/94. & omch2 = omch2_orig - omnuh2
params_to_camb, no_lensing=no_lensing, $
  ombh2=ombh2, omch2=omch2, $
  hubble=hubble, yhe=yhe, neff=neff, $
  as=as, ns=ns, nrun=nrun, nt=nt, r=r, $
  tau=tau, $
  in_l_max_scalar=in_l_max_scalar, $
  in_k_eta_max_scalar=in_k_eta_max_scalar, $
  l=l, dl=tt,phi=phi, $
  ee=ee, bb=bb, te=te, $
  alens=1., accuracy_level=2., omnuh2=omnuh2
tab = [[l],[tt],[te],[ee],[bb]]
writetab,transpose(tab),pdir+'wmap_s12_bestfit_r0p00_summnu0p1eV.txt'

yoyo:
r = 0. & alens=1. & omnuh2=0.5/94. & omch2 = omch2_orig - omnuh2
params_to_camb, no_lensing=no_lensing, $
  ombh2=ombh2, omch2=omch2, $
  hubble=hubble, yhe=yhe, neff=neff, $
  as=as, ns=ns, nrun=nrun, nt=nt, r=r, $
  tau=tau, $
  in_l_max_scalar=in_l_max_scalar, $
  in_k_eta_max_scalar=in_k_eta_max_scalar, $
  l=l, dl=tt,phi=phi, $
  ee=ee, bb=bb, te=te, $
  alens=1., accuracy_level=2., omnuh2=omnuh2
tab = [[l],[tt],[te],[ee],[bb]]
writetab,transpose(tab),pdir+'wmap_s12_bestfit_r0p00_summnu0p5eV.txt'



stop


l_cmb = findgen(4e3)+2.
cl_cmb=wmap7_th_cl(l_cmb,dl=dl_cmb)

; get the residual Poisson power
poisson3000 = y[18+k,wh_min]
help,poisson3000
dl_poisson = poisson3000*((l_cmb/3000.)^2.)

; get the power from SZ
b=read_ascii('~rkeisler/paramfits/cosmomc.s10/data/dl_tsz_sehgal.txt')
yy = b.field1
; match the ell range with that of the cmb
istart = (where_closest(reform(yy[0,*]),2))[0]
istop = istart + n_elements(l_cmb) -1
yy = yy[*,istart:istop]
l_sz = reform(yy[0,*])
dl_sz = reform(yy[1,*])
wh_3000 = (where_closest(l_sz,3000))[0]
sz3000 = y[16+k,wh_min]
dl_sz = dl_sz/dl_sz[wh_3000]*sz3000

; get the power from spatially correlated galaxies and tSZ and kSZ
b=read_ascii('~rkeisler/paramfits/cosmomc.k10/ptsrc/dl_dogleg_flat_0p8.txt')
yy = b.field1
; match the ell range with that of the cmb
istart = (where_closest(reform(yy[0,*]),2))[0]
istop = istart + n_elements(l_cmb) -1
yy = yy[*,istart:istop]
l_flat = reform(yy[0,*])
dl_flat = reform(yy[1,*])
wh_3000 = (where_closest(l_flat,3000))[0]
flat3000 = y[19+k,wh_min]
dl_flat = dl_flat/dl_flat[wh_3000]*flat3000

; combine
dl_total = dl_cmb + dl_poisson + dl_sz + dl_flat



; copies for plotting down to very low ell
dl_totalo = dl_total
dl_cmbo = dl_cmb
l_cmbo = l_cmb
;a = transpose([[l_cmbo],[dl_totalo]]) & writetab,a,'k11_bestfit_total_dl_lcdm.txt'
;b = transpose([[l_cmbo],[dl_cmbo]]) & writetab,b,'k11_bestfit_cmbonly_dl_lcdm.txt'
;stop

; load the spt low-ell spectrum
restore,data
istart = (where_closest(l_cmb,min(l_wf)))[0]
istop = (where_closest(l_cmb,max(l_wf)))[0]
dl_cmb = dl_cmb[istart:istop]
dl_poisson = dl_poisson[istart:istop]
dl_flat = dl_flat[istart:istop]
dl_total = dl_total[istart:istop]
l_cmb = l_cmb[istart:istop]


; only use a subset of the bandpowers
istart=9
istop=55
;stop
cov_all = cov_all[istart:istop,istart:istop]
cov_all_nobeam = cov_all_nobeam[istart:istop,istart:istop]
;cov_all_nocal = cov_all_nocal[istart:istop,istart:istop]
cov_all_nocal = cov_all_nobeam
diag = diag[istart:istop]
diag_nobeam = diag_nobeam[istart:istop]
dl_all = dl_all[istart:istop]
l = l[istart:istop]
wf_all = wf_all[*,istart:istop]
;stop
;print,dl_all
nl = n_elements(dl_all)
dl_th = dblarr(nl)
for i=0,nl-1 do dl_th[i] = total(dl_total*wf_all[*,i])

; convert the data from K to uK2
dl_all *= (1d12)
diag *= (1d12)
diag_nobeam *= (1d12)
cov_all *= (1d24)
cov_all_nobeam *= (1d24)
cov_all_nocal *= (1d24)

icov = invert(cov_all,/double)
icov_nocal = invert(cov_all_nocal,/double)
icov_nobeam = invert(cov_all_nobeam,/double)

; find the calibration factor which is provides the best fit between
; the data and the theory spectrum.
dcal = 0.0001
mincal = 1./1.3
maxcal = 1.*1.3
ncal = ceil((maxcal-mincal)/dcal)
cal = findgen(ncal)*dcal + mincal
chisq = dblarr(ncal)
chisq2 = dblarr(ncal)
chisq3 = dblarr(ncal)
for i=0,ncal-1 do begin
    dl_tmp = dl_all*cal[i]
    delta = dl_tmp-dl_th
    chisq[i] = (delta ## (icov ## delta))[0]
    chisq2[i] = (delta ## (icov_nocal ## delta))[0]
    chisq3[i] = (delta ## (icov_nobeam ## delta))[0]
endfor


chisq_tmp = chisq
min_chisq = min(chisq_tmp)
wh_min=(where(chisq_tmp eq min_chisq))[0]
best_cal = cal[wh_min]
dcal = abs(best_cal - cal[where_closest(chisq_tmp,min_chisq+1.)])
print,'full cov, BEST_CAL: ',best_cal,' +/- ',dcal

chisq_tmp = chisq2
min_chisq = min(chisq_tmp)
wh_min=(where(chisq_tmp eq min_chisq))[0]
best_cal = cal[wh_min]
dcal = abs(best_cal - cal[where_closest(chisq_tmp,min_chisq+1.)])
print,'diag+beam cov, BEST_CAL: ',best_cal,' +/- ',dcal

chisq_tmp = chisq3
min_chisq = min(chisq_tmp)
wh_min=(where(chisq_tmp eq min_chisq))[0]
best_cal = cal[wh_min]
dcal = abs(best_cal - cal[where_closest(chisq_tmp,min_chisq+1.)])
print,'diag cov, BEST_CAL: ',best_cal,' +/- ',dcal


ratio = dl_all*cal[wh_min]/dl_th
err_ratio = diag_nobeam/dl_th



stop
if ~keyword_set(noplot) then begin
;window,0,xsize=2.*400.,ysize=370.
;!p.multi=[0,2,1]
window,0
color_cmb = !eggplant
color_poisson = !darkgreen
color_flat = !orange
color_total = !black
color_data = !blue
thick1 = 1
thick2 = 1
plot,[0],[0],/nodata,xtitle='!12l!X!N', $
  ytitle='D!D!12l!X!N '+textoidl('(\muK^2)'), chars=1.8, $
  xr=[0,3e3],yr=[40,8000],/yst,/xst,/yl
;oplot,l,dl_all*cal[wh_min],color=color_data
;oplot,l,dl_all*cal[wh_min],color=color_data,ps=1
errplot,l,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=2
;pause
;stop
oplot,l_cmb,dl_cmb,color=color_cmb,thick=thick1
oplot,l_cmb,dl_poisson,color=color_poisson,thick=thick1
oplot,l_cmb,dl_flat,color=color_flat,thick=thick1
oplot,l_cmb,dl_total,color=color_total,thick=thick2
errplot,l,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=2



window,1
plot,l,ratio,/yn,/nodata,xtitle='!12l!X!N',ytitle='Data/Theory',yr=[-1,1]*0.2+1,/yst
oplot,[-1,1]*9999999999.,[1,1],thick=2
oplot,l,ratio
errplot,l,(ratio-err_ratio),(ratio+err_ratio)
endif



nbins = n_elements(l)
print,'MIN_CHISQ:',min_chisq
print,'# of bins:',nbins
eff_nparam = 3. + 0.*6. ;3 for foreground and secondaries, plus however "free" the other 6 params are.
eff_nbins = nbins 
eff_dof = eff_nbins - eff_nparam
print,'eff # of bins:',eff_nbins
print,'eff # of params:',eff_nparam
print,'eff DOF: ',eff_dof
pte = mpchitest(min_chisq,eff_dof)
nsigma = mpchitest(min_chisq,eff_dof,/sigma)
print,'PTE: ',pte
print,'NSIGMA: ',nsigma


if keyword_set(doeps) then begin
    size=6.0
    outname = 'figs/bestfit2.eps'
    outname_pdf = 'figs/bestfit2.pdf'
    openplotps,filen=outname,xs=size,ys=0.8*size,/eps
    setcolors2,/sys
    color_data=!blue
;    color=!darkgreen
    wmap_color=!redorange
    wmap_color=!medorange
    chars=1.3
    thick=6.0
    wmap_thick=4.0
    thick1 = 3.0
    thick2 = 3.0
    lchars=0.95
    !p.multi = 0
;xtitle='!12l!X!N/1000'
    xtitle='!12l!X!N'
    ytitle='!12l(l+1)C!D!12l!3!N/2!7p!X'+textoidl(' [\muK^2]')
    sf=1
    plot,[0],[0],/nodata,xtitle=xtitle, $
      ytitle=ytitle, $
      xr=[0,3.1e3]/sf,yr=[40,7000],/yst,/xst,/yl,chars=chars
    
; get WMAP7 binned power spectrum                                                                                   
    readcol,'wmap_binned_tt_spectrum_7yr_v4p1.txt',l_wmap,llo,lhi,dl_wmap,ddl_wmap,ddl_noise,ddl_sv
;wh_wmap = where(l_wmap le 800.)
    wh_wmap = where(l_wmap le 8000. and l_wmap gt 2)
    errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap+ddl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color
    errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color,width=0.008 ;centers
    oplot,[0,0],[1d-20,1d20],color=!black,thick=1
    oplot,[-1e9,1e9],[1,1]*40.,color=!black,thick=1
    
;oplot,l,dl_all*cal[wh_min],color=color_data
;oplot,l,dl_all*cal[wh_min],color=color_data,ps=1
    errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
;pause
    oplot,l_cmbo/sf,dl_cmbo,color=color_cmb,thick=thick1,lines=2
    oplot,l_cmbo/sf,dl_cmbo,thick=thick1,lines=2
    
    
;oplot,l_cmb,dl_poisson,color=color_poisson,thick=thick1
;oplot,l_cmb,dl_flat,color=color_flat,thick=thick1
    oplot,l_cmbo/sf,dl_totalo,color=color_total,thick=thick2
    errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
    
    xyouts,500./sf,3700,'WMAP7',color=wmap_color,chars=chars*1.5,font=1
    xyouts,1600./sf,1100,'SPT,!CThis work',color=color_data,chars=chars*1.5,font=1
    
    closeps
    spawn,'convert '+outname+' '+outname_pdf
    
endif

if keyword_set(dosave) then begin
    save,l_cmb,dl_total,dl_cmb,dl_poisson,dl_sz,dl_flat,filename='bestfit_WMAP7_SPT_lowell.sav'
endif


if keyword_set(stopit) then stop
END


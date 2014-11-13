;;;
; NAME: script_0822
; PURPOSE:
;   General script
;
; NOTES:
; 1) Make ns + other models plot
; 2) plots for cr
;;;

;;;;;;;;;;;;;;;;;;;
; Use this one
PRO ns_plot

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data, CMB+BAO

; nmu, sum_nu, omegaK, w

; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c4_lcdm_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c4_lcdm_camb_w7s12.paramnames'

; neff
type = 'c12_lcdm_neff_camb_w7s12_BAO'
dir_neff = '/data23/hou/lps12/paramfits/chains_final/'+type+'/chains/'
files_neff = file_search(dir_neff+type+'_*.txt')
pname_neff = dir_neff+type+'.paramnames'

; mnu
type = 'c55_lcdm_mnu_pico_w7s12_BAO'
dir_mnu = '/data23/hou/lps12/paramfits/chains_final/'+type+'/post_chains/'
files_mnu = file_search(dir_mnu+type+'_*.txt')
pname_mnu = dir_mnu+type+'.paramnames'

; omk
type = 'c25_lcdm_omk_camb_w7s12_BAO'
dir_omk = '/data23/hou/lps12/paramfits/chains_final/'+type+'/post_chains/'
files_omk = file_search(dir_omk+type+'_*.txt')
pname_omk = dir_omk+type+'.paramnames'

; w
type = 'c30_lcdm_w_camb_w7s12_BAO'
dir_w = '/data23/hou/lps12/paramfits/chains_final/'+type+'/post_chains/'
files_w = file_search(dir_w+type+'_*.txt')
pname_w = dir_w+type+'.paramnames'

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.3
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.20
ddx=.75
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1700 ; xyouts label
ytxt=2300
dytxt = 0.63
bigspt = 1.
thickspt = 1
lfont = 0
ls = [0,1,2,3,4]

;*******************************
; colors
cc = [ pscolors[8], $ ; blue
       pscolors[2], $ ; blue
       pscolors[5], $ ; forest
       pscolors[10],$ ; magenta
       pscolors[3]]   ; orange


!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[0.92, 1.02]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, 0.902, 0.45, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 0.975, -0.16, charsize=1.5, charthick=3, font=0, '!3n!Ds!X!N'

oplot,[0,0],[0,0],color=3
;subsamp = 5

; LCDM
plot_like1dname,files_cmb,pname_cmb,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[0], color=cc[0], /oplot
xyouts,0.985,0.95,'LCDM',charsize=csize, color=cc[0],font=lfont

; neff
plot_like1dname,files_neff,pname_neff,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[1], color=cc[1], /oplot
xyouts,0.985,0.9,'LCDM+neff',charsize=csize, color=cc[1],font=lfont

; mnu
plot_like1dname,files_mnu,pname_mnu,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[2], color=cc[2], /oplot
xyouts,0.985,0.85,'LCDM+mnu',charsize=csize, color=cc[2],font=lfont

; omk
plot_like1dname,files_omk,pname_omk,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[3], color=cc[3], /oplot
xyouts,0.985,0.8,'LCDM+omk',charsize=csize, color=cc[3],font=lfont

; w
plot_like1dname,files_w,pname_w,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[4], color=cc[4], /oplot
xyouts,0.985,0.75,'LCDM+w',charsize=csize, color=cc[4],font=lfont


; vertical line
;oplot, [1.,1.], [0,1], linestyle=1, thick=3

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/ns_extra_0822.pdf'
stop
END







;;;;;;;;;;;;;;;;
; Plots for CR, 1

pro plot_for_cr_1, filename, data=data, $
                   dl_total=dl_total, dl_cmb=dl_cmb, l_cmb=l_cmb, $
                   noplot=noplot, stopit=stopit, ratio=ratio, err_ratio=err_ratio, $
                   l_ratio=l, doeps=doeps, alens=alens, $
                   dosave=dosave, $
                   in_l_max_scalar=in_l_max_scalar, $
                   in_k_eta_max_scalar=in_k_eta_max_scalar

doeps=1
pdir = '/home/kstory/lps12/scripts/plotting/'
if n_elements(filename) eq 0 then $
  filename='/data/kstory/projects/lps12/scripts/paramfits/cons_chain_c4_lcdm_camb_w7s12.sav'
  ;filename=pdir+'cons_chain_baseline_k11.sav'
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
k=0
hubble = y[33+k,wh_min]
as = exp(y[18+k,wh_min])/1d10
ns = y[15+k,wh_min]
nrun = y[17+k,wh_min]
nt = y[16+k,wh_min]
r = y[19+k,wh_min]

in_l_max_scalar = 13000.

;goto, yoyo


; r = 0. & alens=1. & omnuh2=0. & omch2 = omch2_orig - omnuh2
; params_to_camb, no_lensing=no_lensing, $
;   ombh2=ombh2, omch2=omch2, $
;   hubble=hubble, yhe=yhe, neff=neff, $
;   as=as, ns=ns, nrun=nrun, nt=nt, r=r, $
;   tau=tau, $
;   in_l_max_scalar=in_l_max_scalar, $
;   in_k_eta_max_scalar=in_k_eta_max_scalar, $
;   l=l, dl=tt,phi=phi, $
;   ee=ee, bb=bb, te=te, $
;   alens=1., accuracy_level=2., omnuh2=omnuh2
; tab = [[l],[tt],[te],[ee],[bb]]
; writetab,transpose(tab),pdir+'wmap_s12_bestfit_cmbonly_dl_lcdm.txt'


;stop


; OBSOLETE
;l_cmb = findgen(4e3)+2.
;cl_cmb=wmap7_th_cl(l_cmb,dl=dl_cmb)

; Get the best-fit line to this data, which was created above
readcol,pdir+'wmap_s12_bestfit_cmbonly_dl_lcdm.txt', l_cmb, dl_cmb
l_cmb = l_cmb[0:4e3]
dl_cmb = dl_cmb[0:4e3]

; get the residual Poisson power
poisson3000 = 19.8966 ;y[18+k,wh_min]
help,poisson3000
dl_poisson = poisson3000*((l_cmb/3000.)^2.)

; get the power from SZ
;b=read_ascii('~rkeisler/paramfits/cosmomc.s10/data/dl_tsz_sehgal.txt')
b=read_ascii('/home/cr/paramfits/cosmomc.r11/ptsrc/dl_shaw_tsz_s10_153ghz.txt')
yy = b.field1
; match the ell range with that of the cmb
istart = (where_closest(reform(yy[0,*]),2))[0]
istop = istart + n_elements(l_cmb) -1
yy = yy[*,istart:istop]
l_sz = reform(yy[0,*])
dl_sz = reform(yy[1,*])
wh_3000 = (where_closest(l_sz,3000))[0]
sz3000 = y[16+k,wh_min]
;sz3000 = y[21+k,wh_min]
dl_sz = dl_sz/dl_sz[wh_3000]*sz3000

;;; get the power from spatially correlated galaxies and tSZ and kSZ
clust_3000 = 5.41430
dl_flat = clust_3000*((l_cmb/3000.)^(0.8))


; b=read_ascii('~rkeisler/paramfits/cosmomc.k10/ptsrc/dl_dogleg_flat_0p8.txt')
; yy = b.field1
; ; match the ell range with that of the cmb
; istart = (where_closest(reform(yy[0,*]),2))[0]
; istop = istart + n_elements(l_cmb) -1
; yy = yy[*,istart:istop]
; l_flat = reform(yy[0,*])
; dl_flat = reform(yy[1,*])
; wh_3000 = (where_closest(l_flat,3000))[0]
; flat3000 = y[19+k,wh_min]
; dl_flat = dl_flat/dl_flat[wh_3000]*flat3000

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
    outname = 'figs/bestfit2_forcr1.eps'
    outname_pdf = 'figs/bestfit2_forcr1.pdf'
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
    ytitle='!12l(l+1)C!Dl!3!N/2!7p!X!N (!7l!6K!3!U2!X!N)'
    sf=1
    plot,[0],[0],/nodata,xtitle=xtitle, $
      ytitle=ytitle, $
      xr=[0,3.1e3]/sf,yr=[40,7000],/yst,/xst,/yl,chars=chars, charthick=3, color=!black,$
      xthick=5, ythick=5
    
; get WMAP7 binned power spectrum                                                                                   
    readcol,'wmap_binned_tt_spectrum_7yr_v4p1.txt',l_wmap,llo,lhi,dl_wmap,ddl_wmap,ddl_noise,ddl_sv
;wh_wmap = where(l_wmap le 800.)
    wh_wmap = where(l_wmap le 8000. and l_wmap gt 2)
    errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap+ddl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color
    errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color,width=0.008 ;centers
;     oplot,[0,0],[1d-20,1d20],color=!black,thick=1
;     oplot,[-1e9,1e9],[1,1]*40.,color=!black,thick=1
    
;     errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
; ;pause
;     oplot,l_cmbo/sf,dl_cmbo,color=color_cmb,thick=thick1,lines=2
;     oplot,l_cmbo/sf,dl_cmbo,thick=thick1,lines=2
    
    
;     oplot,l_cmbo/sf,dl_totalo,color=color_total,thick=thick2
;     errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
    
    xyouts,500./sf,3700,'WMAP7',color=wmap_color,chars=chars*1.5,font=1
;     xyouts,1600./sf,1100,'SPT, Full Survey',color=color_data,chars=chars*1.5,font=1
    
    closeps
    spawn,'convert '+outname+' '+outname_pdf
    
endif

; if keyword_set(dosave) then begin
;     save,l_cmb,dl_total,dl_cmb,dl_poisson,dl_sz,dl_flat,filename='/home/kstory/lps12/scripts/plotting/bestfit_WMAP7_lps12.sav'
; endif


if keyword_set(stopit) then stop
END






;;;;;;;;;;;;;;;;
; Plots for CR, 2

pro plot_for_cr_2, filename, data=data, $
                   dl_total=dl_total, dl_cmb=dl_cmb, l_cmb=l_cmb, $
                   noplot=noplot, stopit=stopit, ratio=ratio, err_ratio=err_ratio, $
                   l_ratio=l, doeps=doeps, alens=alens, $
                   dosave=dosave, $
                   in_l_max_scalar=in_l_max_scalar, $
                   in_k_eta_max_scalar=in_k_eta_max_scalar

doeps=1
pdir = '/home/kstory/lps12/scripts/plotting/'
if n_elements(filename) eq 0 then $
  filename='/data/kstory/projects/lps12/scripts/paramfits/cons_chain_c4_lcdm_camb_w7s12.sav'
  ;filename=pdir+'cons_chain_baseline_k11.sav'
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
k=0
hubble = y[33+k,wh_min]
as = exp(y[18+k,wh_min])/1d10
ns = y[15+k,wh_min]
nrun = y[17+k,wh_min]
nt = y[16+k,wh_min]
r = y[19+k,wh_min]

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


;stop


; OBSOLETE
;l_cmb = findgen(4e3)+2.
;cl_cmb=wmap7_th_cl(l_cmb,dl=dl_cmb)

; Get the best-fit line to this data, which was created above
readcol,pdir+'wmap_s12_bestfit_cmbonly_dl_lcdm.txt', l_cmb, dl_cmb
l_cmb = l_cmb[0:4e3]
dl_cmb = dl_cmb[0:4e3]

; get the residual Poisson power
poisson3000 = 19.8966 ;y[18+k,wh_min]
help,poisson3000
dl_poisson = poisson3000*((l_cmb/3000.)^2.)

; get the power from SZ
;b=read_ascii('~rkeisler/paramfits/cosmomc.s10/data/dl_tsz_sehgal.txt')
b=read_ascii('/home/cr/paramfits/cosmomc.r11/ptsrc/dl_shaw_tsz_s10_153ghz.txt')
yy = b.field1
; match the ell range with that of the cmb
istart = (where_closest(reform(yy[0,*]),2))[0]
istop = istart + n_elements(l_cmb) -1
yy = yy[*,istart:istop]
l_sz = reform(yy[0,*])
dl_sz = reform(yy[1,*])
wh_3000 = (where_closest(l_sz,3000))[0]
sz3000 = y[16+k,wh_min]
;sz3000 = y[21+k,wh_min]
dl_sz = dl_sz/dl_sz[wh_3000]*sz3000

;;; get the power from spatially correlated galaxies and tSZ and kSZ
clust_3000 = 5.41430
dl_flat = clust_3000*((l_cmb/3000.)^(0.8))


; b=read_ascii('~rkeisler/paramfits/cosmomc.k10/ptsrc/dl_dogleg_flat_0p8.txt')
; yy = b.field1
; ; match the ell range with that of the cmb
; istart = (where_closest(reform(yy[0,*]),2))[0]
; istop = istart + n_elements(l_cmb) -1
; yy = yy[*,istart:istop]
; l_flat = reform(yy[0,*])
; dl_flat = reform(yy[1,*])
; wh_3000 = (where_closest(l_flat,3000))[0]
; flat3000 = y[19+k,wh_min]
; dl_flat = dl_flat/dl_flat[wh_3000]*flat3000

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
    outname = 'figs/bestfit2_forcr2.eps'
    outname_pdf = 'figs/bestfit2_forcr2.pdf'
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
    ytitle='!12l(l+1)C!Dl!3!N/2!7p!X!N (!7l!6K!3!U2!X!N)'
    sf=1
    plot,[0],[0],/nodata,xtitle=xtitle, $
      ytitle=ytitle, $
      xr=[0,3.1e3]/sf,yr=[40,7000],/yst,/xst,/yl,chars=chars, charthick=3, color=!black,$
      xthick=5, ythick=5
    
; get WMAP7 binned power spectrum                                                                                   
    readcol,'wmap_binned_tt_spectrum_7yr_v4p1.txt',l_wmap,llo,lhi,dl_wmap,ddl_wmap,ddl_noise,ddl_sv
;wh_wmap = where(l_wmap le 800.)
    wh_wmap = where(l_wmap le 8000. and l_wmap gt 2)
    errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap+ddl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color
    errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color,width=0.008 ;centers
;     oplot,[0,0],[1d-20,1d20],color=!black,thick=1
;     oplot,[-1e9,1e9],[1,1]*40.,color=!black,thick=1
    
    errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
;pause
;     oplot,l_cmbo/sf,dl_cmbo,color=color_cmb,thick=thick1,lines=2
;     oplot,l_cmbo/sf,dl_cmbo,thick=thick1,lines=2
    
    
;     oplot,l_cmbo/sf,dl_totalo,color=color_total,thick=thick2
;     errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
    
    xyouts,500./sf,3700,'WMAP7',color=wmap_color,chars=chars*1.5,font=1
    xyouts,1600./sf,1100,'SPT, Full Survey',color=color_data,chars=chars*1.5,font=1
    
    closeps
    spawn,'convert '+outname+' '+outname_pdf
    
endif

; if keyword_set(dosave) then begin
;     save,l_cmb,dl_total,dl_cmb,dl_poisson,dl_sz,dl_flat,filename='/home/kstory/lps12/scripts/plotting/bestfit_WMAP7_lps12.sav'
; endif


if keyword_set(stopit) then stop
END




;;;;;;;;;;;;;;;;;;
; Procedure to make bandpower comparison plot for CR
; 
PRO plot_dl_all_2

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data
pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
uk2mk = 1d-6

; S12
restore,file
dl_all_lps12 = dl_all*(1d12)
diag_nobeam_lps12 = diag_nobeam*(1d12)
l_lps12 = l

istart=9
istop=55
dl_all_lps12 = dl_all_lps12[istart:istop]
diag_nobeam_lps12 = diag_nobeam_lps12[istart:istop]
l_lps12 = l_lps12[istart:istop]
cl4_lps12 = dl_all_lps12 *(l_lps12^4.) / (l_lps12*(l_lps12+1)) * uk2mk
dcl4_lps12 = diag_nobeam_lps12 *(l_lps12^4.) / (l_lps12*(l_lps12+1)) * uk2mk

; load acbar
readcol,pdir+'cl_dc1_v3.dat',iband,dl_acbar,delta_dl_acbar,log_normal_offset
ledge=[[350,550],[550,650],[650,730],[730,790],[790,850],[850,910],[910,970],[970,1030],[1030,1090],[1090,1150],[1150,1210],[1210,1270],[1270,1330],[1330,1390],[1390,1450],[1450,1510],[1510,1570],[1570,1650],[1650,1750],[1750,1850],[1850,1950],[1950,2100],[2100,2300],[2300,2500],[2500,3000]]
llo=reform(ledge[0,*])
lhi=reform(ledge[1,*])
l_acbar=0.5*(llo+lhi)
ddl_acbar = delta_dl_acbar

cl4_acbar = dl_acbar *(l_acbar^4.) / (l_acbar*(l_acbar+1)) * uk2mk
dcl4_acbar = ddl_acbar *(l_acbar^4.) / (l_acbar*(l_acbar+1)) * uk2mk
; drop last 2 bins:
cl4_acbar = cl4_acbar[0:n_elements(cl4_acbar)-2]
dcl4_acbar = dcl4_acbar[0:n_elements(cl4_acbar)-2]


; ; load  K11
; restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
; cal_rk = 1;0.76
; dl_all_k11      = dl_all*1d12*cal_rk^2.
; diag_nobeam_k11 = diag_nobeam*1d12*cal_rk^2.
; l_k11 = l

; dl_all_k11 = dl_all_k11[istart:istop]
; diag_nobeam_k11 = diag_nobeam_k11[istart:istop]
; l_k11 = l_k11[istart:istop]
; cl4_k11 = dl_all_k11 *(l_k11^4.) / (l_k11*(l_k11+1)) * uk2mk
; dcl4_k11 = diag_nobeam_k11 *(l_k11^4.) / (l_k11*(l_k11+1)) * uk2mk

sf=1.

; load act (das et al)
readcol,pdir+'act_150_use.txt',l_act,dl_act,ddl_act
l_act   = l_act[3:39]
dl_act  = dl_act[3:39]
ddl_act = ddl_act[3:39]
cl4_act = dl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk
dcl4_act = ddl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk

; load theory
; restore, '/home/kstory/lps12/scripts/plotting/bestfit_WMAP7_lps12.sav'
; cl4_fit = dl_total *(l_cmb^4.) / (l_cmb*(l_cmb+1)) * uk2mk

readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)
cl4_th = dl_th *(l_vec^4.) / (l_vec*(l_vec+1)) * uk2mk
cl4_fit = cl4_th
l_cmb = l_vec

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.3
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.20
ddx=.75
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3  ; points
pthick=5 ; points
thick=3 ; points
xtxt=1700 ; xyouts label
ytxt=2300
dytxt = 0.63
bigspt = 1.5
thickspt = 1.7
lfont = 0
ls = [0,2]

;*******************************
; colors; use setup_ps_plotting
color_th = pscolors[0]        ; black
color_act = pscolors[2]       ; red
color_acbar = pscolors[5]       ; green
color_this_work = pscolors[0] ; black

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[450, 3200]
!y.window=[y0,y0+ddy]
!y.crange=[0, 2000]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=5
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=5
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=5
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=5

; axes labels
xtitle='!12l!X!N'
ytitle='!12l!U4!N!12C!Dl!3!N/2!7p!X!N (!6mK!3!U2!X!N)'
xyouts, 50, 1000, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle
xyouts, 1800, -330, charsize=1.4, charthick=3, xtitle

oplot,[0,0],[0,0],color=3

wh = where(l_cmb ge 450 and l_cmb lt 3200)
oplot,l_cmb[wh],cl4_fit[wh],color=color_th, thick=thick

errplot,l_act,cl4_act-dcl4_act,cl4_act+dcl4_act,thick=pthick,color=color_act
errplot,l_acbar,(cl4_acbar-dcl4_acbar),(cl4_acbar+dcl4_acbar),thick=pthick,color=color_acbar
errplot,l_lps12,(cl4_lps12-dcl4_lps12),(cl4_lps12+dcl4_lps12),thick=7,color=color_this_work


xyouts,xtxt,1800,         'ACT',chars=bigspt,charthick=thickspt,font=lfont,color=color_act
xyouts,xtxt,1650,   'ACBAR',chars=bigspt, charthick=thickspt,font=lfont,color=color_acbar
xyouts,xtxt,1500,'SPT, Full Survey',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work

;*******************************
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/dl_all2.pdf'

stop
END









;;;;;;;;;;;;;;;;;;;;;;
; Get CDF limits on ns
; Power-law [LCDM]
;   Print constraints for all LCDM parameters using this procedure
;   type in ['s12', 'w7', 'cmb', 'extra']
;;;;;;;;;;;;;;;;;;;;;;
PRO nrun_cdf, type=type, use_r=use_r

if n_elements(type) eq 0 then type='cmb'

;;; LCDM+nrun
; CMB+BAO
dir = '/data23/hou/lps12/paramfits/chains_final/c13_lcdm_nrun_camb_w7s12_BAO/chains/'
files = file_search(dir+'c13_lcdm_nrun_camb_w7s12_BAO*.txt')
pname = dir+'c13_lcdm_nrun_camb_w7s12_BAO.paramnames'


subsamp=10000
nskip=5000
stale=1
plot_like1dname,files,pname,'nrun',subsamp=subsamp,nskip=nskip,scale=scale,/cdf,/stopit

;; Run the following on the command line when the program stops
;wh = where(bins ge 1.0) & print, ff[wh[0]] & print, gauss_cvf(ff[wh[0]])

stop
END



PRO ns_cdf, type=type, use_r=use_r
; neff
type = 'c12_lcdm_neff_camb_w7s12_BAO'
dir = '/data23/hou/lps12/paramfits/chains_final/'+type+'/chains/'
files = file_search(dir+type+'_*.txt')
pname = dir+type+'.paramnames'

subsamp=10000
nskip=5000
stale=1
plot_like1dname,files,pname,'ns',subsamp=subsamp,nskip=nskip,scale=scale,/cdf,/stopit

;; Run the following on the command line when the program stops
;wh = where(bins ge 1.0) & print, ff[wh[0]] & print, gauss_cvf(ff[wh[0]])
stop
END

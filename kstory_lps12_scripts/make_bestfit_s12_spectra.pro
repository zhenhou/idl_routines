;;;
; NAME: make_bestfit_s12_spectra
; PURPOSE:
;   Make the best-fit SPT+WMAP spectra, and plot
;
; NOTES:
;
; MODIFICATION HISTORY:
;  07/20/2012: (KTS) Created
;  08/16/2012: (KTS) Update foreground components
;  08/22/2012: (KTS) Change plot to use best-fit line, rather than wmap7 line.
;  08/23/2012: (KTS) Re-write function to use CR's ML valuse
;  09/12/2012: (KTS) Use 0828 bandpowers and chains
;;;



;-------------------------------
;;;;;;;;;;;;;;;;;;
; Make best-fit spectra for SPT+WMAP
;   output, bestfit2.pdf
;;;;;;;;;;;;;;;;;;
pro make_bestfit_s12_spectra, doeps=doeps, dosave=dosave, stopit=stopit, $
                              data=data

pdir = '/home/kstory/lps12/scripts/plotting/'

if n_elements(data) eq 0 then $
  data = '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120828_170101_kweight.sav'

; Get ML values, from params_c2_lcdm_pico_w7s12_ML.ini
poisson3000 = 19.8883
clust_3000 = 4.37709
sz3000 = 5.10884

; Get ML values, from params_c4_lcdm_camb_w7s12_ML.ini
;poisson3000 = 20.0976 & clust_3000 = 5.46899 & sz3000 = 4.18801 ; run_08

;---------------------
; Get the best-fit spectrum including foregrounds
;---------------------
; readcol, '/home/kstory/lps12/best_fit/suxp_cls', cl_total, te, ee
; ncl = n_elements(cl_total)
; l_total = indgen(ncl)+2
; dl_total = cl_total*l_total*(l_total+1)/(2*!pi)
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', l_total, dl_total

;---------------------
; Get the foregrounds
;---------------------

;---------------------
; get the residual Poisson power
help,poisson3000
dl_poisson = poisson3000*((l_total/3000.)^2.)

;---------------------
; get the power from SZ
b=read_ascii('/home/cr/paramfits/cosmomc.r11/ptsrc/dl_shaw_tsz_s10_153ghz.txt')
yy = b.field1
; match the ell range with that of the cmb
istart = (where_closest(reform(yy[0,*]),2))[0]
istop = istart + n_elements(l_total) -1
yy = yy[*,istart:istop]
l_sz = reform(yy[0,*])
dl_sz = reform(yy[1,*])
wh_3000 = (where_closest(l_sz,3000))[0]
help, sz3000
dl_sz = dl_sz/dl_sz[wh_3000]*sz3000

wh_lstart = where( l_sz eq l_total[0])
l_sz = l_sz[wh_lstart:*]
dl_sz = dl_sz[wh_lstart:*]

;---------------------
; get the power from spatially correlated galaxies and tSZ and kSZ
help, clust_3000
dl_flat = clust_3000*((l_total/3000.)^(0.8))

; subtract foregrounds
dl_cmb = dl_total - dl_poisson - dl_sz - dl_flat



; copies for plotting down to very low ell
dl_totalo = dl_total
dl_cmbo = dl_cmb
l_cmbo = l_total

; load the spt low-ell spectrum
restore,data
istart = (where_closest(l_total,min(l_wf)))[0]
istop = (where_closest(l_total,max(l_wf)))[0]
dl_cmb = dl_cmb[istart:istop]
dl_poisson = dl_poisson[istart:istop]
dl_flat = dl_flat[istart:istop]
dl_total = dl_total[istart:istop]
l_total = l_total[istart:istop]


; only use a subset of the bandpowers
istart=9
istop=55
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


;---------------------
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



;---------------------
; Plot
;---------------------

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
  color=!black, xthick=1, ythick=1, $
  xr=[0,3e3],yr=[40,8000],/yst,/xst,/yl
;oplot,l,dl_all*cal[wh_min],color=color_data
;oplot,l,dl_all*cal[wh_min],color=color_data,ps=1
errplot,l,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=2
;pause
;stop
oplot,l_total,dl_cmb,color=color_cmb,thick=thick1
oplot,l_total,dl_poisson,color=color_poisson,thick=thick1
oplot,l_total,dl_flat,color=color_flat,thick=thick1
oplot,l_total,dl_total,color=color_total,thick=thick2
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
    xtitle='!12l!X!N'
    ytitle='!12D!Dl!X!N (!7l!6K!3!U2!X!N)'
;     xtitle='!12l!X!N'
;     ytitle='!8D!D!12l!X!N (!9m!X!NK!U2!X!N)'
    ;ytitle='!8l(l+1)C!D!8l!3!N/2!9p!X!N [!9m!X!NK!U2!X!N]'
    ;font=0
    sf=1
    plot,[0],[0],/nodata,xtitle=xtitle, $
      ytitle=ytitle, $
      color=!black, xthick=5, ythick=5, $
      xr=[0,3.1e3]/sf,yr=[40,7000],/yst,/xst,/yl,chars=chars, charthick=3

    
; get WMAP7 binned power spectrum                                                                                   
    readcol,'wmap_binned_tt_spectrum_7yr_v4p1.txt',l_wmap,llo,lhi,dl_wmap,ddl_wmap,ddl_noise,ddl_sv
    wh_wmap = where(l_wmap le 8000. and l_wmap gt 2)
    errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap+ddl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color
    errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color,width=0.008 ;centers
    oplot,[0,0],[1d-20,1d20],color=!black,thick=1
    oplot,[-1e9,1e9],[1,1]*40.,color=!black,thick=1
    
    errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
;pause
    oplot,l_cmbo/sf,dl_cmbo,color=color_cmb,thick=thick1,lines=2
    oplot,l_cmbo/sf,dl_cmbo,thick=thick1,lines=2
    

    oplot,l_cmbo/sf,dl_totalo,color=color_total,thick=thick2
    errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
    
    xyouts,500./sf,3700,'WMAP7',color=wmap_color,chars=chars*1.5,font=0
    xyouts,1500./sf,1100,'SPT',color=color_data,chars=chars*1.5,font=0
    
    closeps
    spawn,'convert '+outname+' '+outname_pdf
    
endif

if keyword_set(dosave) then begin
    save,l_total,dl_total,dl_cmb,dl_poisson,dl_sz,dl_flat,best_cal,filename='/home/kstory/lps12/scripts/plotting/bestfit_WMAP7_lps12.sav'
endif


if keyword_set(stopit) then stop
END













;;;;;;;;;;;;;;;;;;
; Make best-fit spectra for SPT+WMAP
;    Code originally from RK
;;;;;;;;;;;;;;;;;;
pro make_bestfit_s12_spectra_rk, filename, data=data, $
                              dl_total=dl_total, dl_cmb=dl_cmb, l_cmb=l_cmb, $
                              noplot=noplot, stopit=stopit, ratio=ratio, err_ratio=err_ratio, $
                              l_ratio=l, doeps=doeps, alens=alens, $
                              dosave=dosave, $
                              in_l_max_scalar=in_l_max_scalar, $
                              in_k_eta_max_scalar=in_k_eta_max_scalar

pdir = '/home/kstory/lps12/scripts/plotting/'
if n_elements(filename) eq 0 then $
  filename='/data/kstory/projects/lps12/scripts/paramfits/cons_chain_c4_lcdm_camb_w7s12.sav'

if n_elements(data) eq 0 then $
  data = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'

restore,filename
y = echain

lnl = reform(y[1,*])
min_lnl = min(lnl)
wh_min = where(lnl eq min(lnl), nmin)
if nmin gt 1 then begin
    print,'yo: there were more than one points in the chain with the minimum log likelihood'
endif
wh_min = wh_min[0]

;---------------------
; From RK
; ombh2 = y[2,wh_min]
; omch2 = y[3,wh_min]
; yhe = y[10,wh_min]
; neff = y[9,wh_min]
; tau = y[5,wh_min]
; omch2_orig = omch2
; ;k = keyword_set(alens) ? 1 : 0
; k=1
; hubble = y[42+k,wh_min]
; as = exp(y[14+k,wh_min])/1d10
; ns = y[11+k,wh_min]
; nrun = y[13+k,wh_min]
; nt = y[12+k,wh_min]
; r = y[15+k,wh_min]
;---------------------

ombh2 = y[2,wh_min]
omch2 = y[3,wh_min]
yhe = y[10,wh_min]
neff = y[9,wh_min]
tau = y[5,wh_min]
omch2_orig = omch2

hubble = y[33,wh_min]
as = exp(y[18,wh_min])/1d10
ns = y[15,wh_min]
nrun = y[17,wh_min]
nt = y[16,wh_min]
r = y[19,wh_min]

in_l_max_scalar = 13000.

stop
;---------------------
; Calculate the best-fit CMB spectrum 
;---------------------

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
writetab,transpose(tab),pdir+'wmap_s12_bestfit_rk_cmbonly_dl_lcdm.txt'

;stop

;---------------------
; Make the best-fit spectrum including foregrounds
;---------------------

;---------------------
; Get the best-fit CMB spectrum (created above)
readcol,pdir+'wmap_s12_bestfit_rk_cmbonly_dl_lcdm.txt', l_cmb, dl_cmb
l_cmb = l_cmb[0:4e3]
dl_cmb = dl_cmb[0:4e3]

;---------------------
; get the residual Poisson power
poisson3000 =  y[22,wh_min]
help,poisson3000
dl_poisson = poisson3000*((l_cmb/3000.)^2.)

;---------------------
; get the power from SZ
b=read_ascii('/home/cr/paramfits/cosmomc.r11/ptsrc/dl_shaw_tsz_s10_153ghz.txt')
yy = b.field1
; match the ell range with that of the cmb
istart = (where_closest(reform(yy[0,*]),2))[0]
istop = istart + n_elements(l_cmb) -1
yy = yy[*,istart:istop]
l_sz = reform(yy[0,*])
dl_sz = reform(yy[1,*])
wh_3000 = (where_closest(l_sz,3000))[0]
sz3000 = y[20,wh_min]
help, sz3000
dl_sz = dl_sz/dl_sz[wh_3000]*sz3000

wh_lstart = where( l_sz eq l_cmb[0])
l_sz = l_sz[wh_lstart:*]
dl_sz = dl_sz[wh_lstart:*]

;---------------------
; get the power from spatially correlated galaxies and tSZ and kSZ
clust_3000 = y[23,wh_min]
help, clust_3000
dl_flat = clust_3000*((l_cmb/3000.)^(0.8))

; combine
dl_total = dl_cmb + dl_poisson + dl_sz + dl_flat



; copies for plotting down to very low ell
dl_totalo = dl_total
dl_cmbo = dl_cmb
l_cmbo = l_cmb

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


;---------------------
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



;---------------------
; Plot
;---------------------

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
  color=!black, xthick=1, ythick=1, $
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
    outname = 'figs/bestfit2_rk.eps'
    outname_pdf = 'figs/bestfit2_rk.pdf'
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
      color=!black, xthick=5, ythick=5, $
      xr=[0,3.1e3]/sf,yr=[40,7000],/yst,/xst,/yl,chars=chars, charthick=3
    
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
    xyouts,1600./sf,1100,'SPT, Full Survey',color=color_data,chars=chars*1.5,font=1
    
    closeps
    spawn,'convert '+outname+' '+outname_pdf
    
endif

if keyword_set(dosave) then begin
    save,l_cmb,dl_total,dl_cmb,dl_poisson,dl_sz,dl_flat,filename='/home/kstory/lps12/scripts/plotting/bestfit_WMAP7_lps12_rk.sav'
endif


if keyword_set(stopit) then stop
END








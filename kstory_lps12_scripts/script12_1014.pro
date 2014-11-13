;;;
; Make best-fit spectrum with linear y-axis
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

doeps = 1

; Get ML values, from params_c2_lcdm_pico_w7s12_ML.ini
poisson3000 = 19.8883
clust_3000 = 4.37709
sz3000 = 5.10884

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
    outname = 'figs/ps_linear.eps'
    outname_pdf = 'figs/ps_linear.pdf'
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
    ;ytitle='!12l(l+1)C!Dl!3!N/2!7p!X!N (!7l!6K!3!U2!X!N)'
    ytitle='!12D!Dl!X!N (!7l!6K!3!U2!X!N)'
    sf=1
    plot,[0],[0],/nodata,xtitle=xtitle, $
      ytitle=ytitle, $
      color=!black, xthick=5, ythick=5, $
      xr=[0,3.1e3]/sf,yr=[40,6000],/yst,/xst,chars=chars, charthick=3
    
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
    
    xyouts,500./sf,3700,'WMAP7',color=wmap_color,chars=chars*1.5,font=1
    xyouts,1600./sf,1300,'SPT-SZ Survey',color=color_data,chars=chars*1.5,font=1
    xyouts,1600./sf,900,'Story et al., in prep',color=color_data,chars=chars*1.5,font=1
    
    closeps
    spawn,'convert '+outname+' '+outname_pdf
    
endif

if keyword_set(dosave) then begin
    save,l_total,dl_total,dl_cmb,dl_poisson,dl_sz,dl_flat,best_cal,filename='/home/kstory/lps12/scripts/plotting/bestfit_WMAP7_lps12.sav'
endif


if keyword_set(stopit) then stop
END



;;;
; NAME: script_0823
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Pipeline test2
; 2) make best-fit lines
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
;;; 2) Pipeline test2, spec1
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_spec1
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_08p2_spec1_kweight_calib.sav'

; get pipe2
restore, file
calib = 1.
;calib = 0.825
;calib = 0.76
calib = 0.85
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
oplot, l[vec], dl_p2[vec], color=!red, linestyle=0
;oplot, l[vec], dl_all[vec]

oplot, l[vec], dl_th_2[vec], color=!green, linestyle=0
oplot, x_th, dl_th_plot, color=!blue, linestyle=2

window, 2
plot, l[vec], (dl_th_2/dl_p2)[vec]
stop
END



;;;;;;;;;;;;;;;;;;;;;;;
;;; 2) Pipeline test2, spec2
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_spec2
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_08p2_spec2_kweight_calib.sav'

; get pipe2
restore, file
calib = 1.
;calib = 0.825
;calib = 0.76
;calib = 0.85
dl_p2 = spectrum*(1d12)^2. * (calib)
cov_p2 = cov

; get spec1
restore, edir+'end_ra3h30dec-60_08p2_spec1_kweight_calib.sav'
dl_p2_spec1 = spectrum*(1d12)^2. * (calib)

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
;x_th = x_th - 10 ; shift back 10
dl_th_plot = dl_th[x_th]

dl_th_2 = dl_th[50:3998]
dl_th_2 = reform(dl_th_2 # wf2use)

; Plotting
vec=indgen(47) + 9

window, 1
plot, l[vec], dl_p2[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl
oplot, l[vec], dl_p2[vec], color=!red, linestyle=0
;oplot, l[vec], dl_all[vec]

;plot, l[vec], dl_p2_spec1[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl
oplot, l[vec], dl_p2_spec1[vec], color=!blue, linestyle=1

; oplot, l[vec], dl_th_2[vec], color=!green, linestyle=0
; oplot, x_th, dl_th_plot, color=!purple, linestyle=2

window, 2
plot, l[vec], (dl_th_2/dl_p2)[vec]
stop
END




;;;;;;;;;;;;;;;;;;;
; Make fields overplot
;;;;;;;;;;;;;;;;;;;
PRO mk_fields_plot
;;;
; Make a plot of the spt fields for observation
; Copied from TC: spt:/home/tcrawfor/spt_analysis/temp/mkplot_spt_obsreg_2010version_batch
; 09/23/2010 (KTS): Created
;;;

; Dust plot
; dusttemp=readfits('/data/tcrawfor/Idl/forecast/SFD_i100_healpix_256.fits')
; orthview,/online,dusttemp,coord=['g','q'],grat=[15,10],colt=3,min=-1.,max=20.,eul_mat=emtemp,rot=[195,225,0],/half,title=' ',/nobar

; WMAP plot
read_fits_map,'/data/tcrawfor/WMAP_maps/wmap_ilc_7yr_v4.fits',dusttemp
dusttemp = dusttemp[*,0]
orthview,/online,/nest,dusttemp*1e3,coord=['g','q'],grat=[15,10],colt=13,min=-200.,max=200.,eul_mat=emtemp,rot=[195,225,0],/half,title=' ',/nobar


;legstr=['4000 sq. deg.','Published fields (Lueker 09','Finished fields','Future fields (2010, 2011)']
legstr=['2008 fields','2009 fields','2010 fields','2011 fields']
loadct,39

;DEVICE, SET_FONT = 'Helvetica*Bold' 
;DEVICE, SET_FONT = 'Palatino*Bold*Italic' 

; color names, 190=yellow, 140=green, 250=red, 40=purple, 80=blue
;color_full=190
color_2008=140
color_2009=250
color_2010=80
color_2011=190
legend,legstr,line=0,psym=0,thick=2,color=[color_2008,color_2009,color_2010,color_2011],position=[-1.01,1.25],charsize=1.9,charthick=1.5,font=-1

; Plot the full 4000 field
; ra_full=[20.,7.] & dec_full=[-65.,-30.]
; oplot_rectangle,ra_full,dec_full,emtemp,proj='ORTH',/half,color=color_full,thick=2

uthick=6
lthick=1

; Plot 2008 fields
oplot_rectangle,[5.,6.],[-60.,-50.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[5.,6.],[-60.,-50.],emtemp,proj='ORTH',/half,color=color_2008,thick=lthick,/crop
oplot_rectangle,[23.,0.],[-60.,-50.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[23.,0.],[-60.,-50.],emtemp,proj='ORTH',/half,color=color_2008,thick=lthick,/crop

; Plot 2009 fields
oplot_rectangle,[20.,22.],[-65.,-55.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[2.,5.],[-65.,-55.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[20.,22.],[-55.,-45.],emtemp,proj='ORTH',/half,thick=uthick,/crop

oplot_rectangle,[20.,22.],[-65.,-55.],emtemp,proj='ORTH',/half,color=color_2009,thick=lthick,/crop
oplot_rectangle,[2.,5.],[-65.,-55.],emtemp,proj='ORTH',/half,color=color_2009,thick=lthick,/crop
oplot_rectangle,[20.,22.],[-55.,-45.],emtemp,proj='ORTH',/half,color=color_2009,thick=lthick,/crop

; Plot 2010 and 2011 fields
oplot_rectangle,[0.,2.],[-65.,-55.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[0.,25.]/15.,[-55.,-45.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[25.,50.]/15.,[-55.,-45.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[50.,75.]/15.,[-55.,-45.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[5.,6.],[-50.,-40.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[6.,7.],[-60.,-50.],emtemp,proj='ORTH',/half,thick=uthick,/crop

oplot_rectangle,[0.,2.],[-65.,-55.],emtemp,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
oplot_rectangle,[0.,25.]/15.,[-55.,-45.],emtemp,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
oplot_rectangle,[25.,50.]/15.,[-55.,-45.],emtemp,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
oplot_rectangle,[50.,75.]/15.,[-55.,-45.],emtemp,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
oplot_rectangle,[5.,6.],[-50.,-40.],emtemp,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
oplot_rectangle,[6.,7.],[-60.,-50.],emtemp,proj='ORTH',/half,color=color_2010,thick=lthick,/crop

; Plot 2011 fields
oplot_rectangle,[20.,22.],[-45.,-40.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[22.,0.],[-65.,-60.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[22.,23.],[-60.,-50.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[22.,0.],[-50.,-40.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[0.,2.],[-45.,-40.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[2.,5.],[-45.,-40.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[5.,7.],[-65.,-60.],emtemp,proj='ORTH',/half,thick=uthick,/crop
oplot_rectangle,[6.,7.],[-50.,-40.],emtemp,proj='ORTH',/half,thick=uthick,/crop

oplot_rectangle,[20.,22.],[-45.,-40.],emtemp,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
oplot_rectangle,[22.,0.],[-65.,-60.],emtemp,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
oplot_rectangle,[22.,23.],[-60.,-50.],emtemp,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
oplot_rectangle,[22.,0.],[-50.,-40.],emtemp,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
oplot_rectangle,[0.,2.],[-45.,-40.],emtemp,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
oplot_rectangle,[2.,5.],[-45.,-40.],emtemp,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
oplot_rectangle,[5.,7.],[-65.,-60.],emtemp,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
oplot_rectangle,[6.,7.],[-50.,-40.],emtemp,proj='ORTH',/half,color=color_2011,thick=lthick,/crop

;annotate,load_file='/home/kstory//projects/field_obs/preview_2010_08/annotate_preview_0923.dat'
;annotate,load_file='/home/tcrawfor/SPTData/annotate_obsreg.dat'
annotate,load_file='/home/kstory/lps12/scripts/plotting/annotate_fields_lps12.dat'


END




;-------------------------------
;;;;;;;;;;;;;;;;;;
; Make best-fit spectra for SPT+WMAP
;;;;;;;;;;;;;;;;;;
pro make_bestfit_s12_spectra, doeps=doeps, dosave=dosave, stopit=stopit, $
                              data=data

pdir = '/home/kstory/lps12/scripts/plotting/'

if n_elements(data) eq 0 then $
  data = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'

; Get ML values, from params_c4_lcdm_camb_w7s12_ML.ini
poisson3000 = 20.0976
clust_3000 = 5.46899
sz3000 = 4.18801

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
    xyouts,1600./sf,1100,'SPT, Full Survey',color=color_data,chars=chars*1.5,font=1
    
    closeps
    spawn,'convert '+outname+' '+outname_pdf
    
endif

if keyword_set(dosave) then begin
    save,l_total,dl_total,dl_cmb,dl_poisson,dl_sz,dl_flat,filename='/home/kstory/lps12/scripts/plotting/bestfit_WMAP7_lps12.sav'
endif


if keyword_set(stopit) then stop
END








;;;;
; Write best-fit model to txt file
;;;;
PRO write_bestfit
readcol, '/home/kstory/lps12/best_fit/suxp_cls', cl_total, te, ee
ncl = n_elements(cl_total)
l_total = indgen(ncl)+2
dl_total = cl_total*l_total*(l_total+1)/(2*!pi)

get_lun, lun1
openw, lun1, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt'
for i=0, ncl-1 do begin
    printf, lun1, l_total[i], dl_total[i]
endfor
close, lun1
free_lun, lun1
END






;-------------------------------
;;;;;;;;;;;;;;;;;;
; Plot the residuals relative to the WMAP best-fit
;;;;;;;;;;;;;;;;;;
pro plot_residuals_wmapfit, doeps=doeps, dosave=dosave, stopit=stopit, $
                              data=data

pdir = '/home/kstory/lps12/scripts/plotting/'

if n_elements(data) eq 0 then $
  data = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'


if n_elements(filename) eq 0 then filename='/home/rkeisler/ps09/cons_chain_wmaponly.sav'
;if n_elements(filename) eq 0 then filename='/home/rkeisler/ps09/cons_chain_baseline.sav'
if n_elements(data) eq 0 then $
  data = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
;  data='/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
; data='/data/rkeisler/ps09/1.24/combined_spectrum_20110207_035848_kweight.sav'

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

;k = keyword_set(alens) ? 1 : 0
k=1
hubble = y[42+k,wh_min]
as = exp(y[14+k,wh_min])/1d10
ns = y[11+k,wh_min]
nrun = y[13+k,wh_min]
nt = y[12+k,wh_min]
r = y[15+k,wh_min]

; get the CMB spectrum
params_to_camb, no_lensing=no_lensing, $
  ombh2=ombh2, omch2=omch2, $
  hubble=hubble, yhe=yhe, neff=neff, $
  as=as, ns=ns, nrun=nrun, nt=nt, r=r, $
  tau=tau, $
  in_l_max_scalar=in_l_max_scalar, $
  in_k_eta_max_scalar=in_k_eta_max_scalar, $
  l=l_cmb, dl=dl_cmb,phi=phi
;stop



;params_to_camb,ombh2=0.02258, omch2=0.1109, hubble=71., yhe=0.2486,neff=3.04,as=2.43e-9,ns=(1.-0.037),nrun=0.,nt=0.,r=0.,tau=0.088,l=l_cmb,dl=dl_wmap
;stop

;params_to_camb, no_lensing=1, $
;  ombh2=ombh2, omch2=omch2, $
;  hubble=hubble, yhe=yhe, neff=neff, $
;  as=as, ns=ns, nrun=nrun, nt=nt, r=r, $
;  tau=tau, $
;  l=l_cmb, dl=dl_cmb_nolens
;stop

;tempp
;l_cmb = findgen(4e3)+2.
;cl_cmb=wmap7_th_cl(l_cmb,dl=dl_cmb)

; ; Get ML values, from params_c4_lcdm_camb_w7s12_ML.ini
poisson3000 = 20.0976
clust_3000  = 5.46899
sz3000      = 4.18801

; get the residual Poisson power
help,poisson3000
dl_poisson = poisson3000*((l_cmb/3000.)^2.)

l_total = l_cmb

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

; combine
dl_total = dl_cmb + dl_poisson + dl_sz + dl_flat


get_lun, lun1
openw, lun1, '/home/kstory/lps12/best_fit/lcdm_WMAP_ML_dl_fromRK.txt'
for i=0, n_elements(dl_total)-1 do begin
    printf, lun1, l_cmb[i], dl_cmb[i], dl_total[i]
endfor
close, lun1
free_lun, lun1



; ; Get ML values, from params_c4_lcdm_camb_w7s12_ML.ini
; poisson3000 = 59.6585 ;20.0976
; clust_3000  = 43.7112 ;5.46899
; sz3000      = 55.1913 ;4.18801

;---------------------
; Get the best-fit spectrum including foregrounds
;---------------------
; readcol, '/home/kstory/lps12/best_fit/suxp_cls', cl_total, te, ee
; ncl = n_elements(cl_total)
; l_total = indgen(ncl)+2
; dl_total = cl_total*l_total*(l_total+1)/(2*!pi)
; readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', l_total, dl_total
; dl_foreground = dl_total - dl_cmb

readcol, '/home/kstory/lps12/best_fit/lcdm_WMAP_ML_dl_fromRK.txt', l_cmb, dl_cmb, dl_total
l_total = l_cmb

; copies for plotting down to very low ell
dl_totalo = dl_total
dl_cmbo = dl_cmb
l_cmbo = l_total

; load the spt low-ell spectrum
restore,data
istart = (where_closest(l_total,min(l_wf)))[0]
istop = (where_closest(l_total,max(l_wf)))[0]
dl_cmb = dl_cmb[istart:istop]
;dl_poisson = dl_poisson[istart:istop]
;dl_flat = dl_flat[istart:istop]
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
; oplot,l_total,dl_poisson,color=color_poisson,thick=thick1
; oplot,l_total,dl_flat,color=color_flat,thick=thick1
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


; if keyword_set(doeps) then begin
;     size=6.0
;     outname = 'figs/bestfit2.eps'
;     outname_pdf = 'figs/bestfit2.pdf'
;     openplotps,filen=outname,xs=size,ys=0.8*size,/eps
;     setcolors2,/sys
;     color_data=!blue
; ;    color=!darkgreen
;     wmap_color=!redorange
;     wmap_color=!medorange
;     chars=1.3
;     thick=6.0
;     wmap_thick=4.0
;     thick1 = 3.0
;     thick2 = 3.0
;     lchars=0.95
;     !p.multi = 0
; ;xtitle='!12l!X!N/1000'
;     xtitle='!12l!X!N'
;     ytitle='!12l(l+1)C!Dl!3!N/2!7p!X!N (!7l!6K!3!U2!X!N)'
;     sf=1
;     plot,[0],[0],/nodata,xtitle=xtitle, $
;       ytitle=ytitle, $
;       color=!black, xthick=5, ythick=5, $
;       xr=[0,3.1e3]/sf,yr=[40,7000],/yst,/xst,/yl,chars=chars, charthick=3
    
; ; get WMAP7 binned power spectrum                                                                                   
;     readcol,'wmap_binned_tt_spectrum_7yr_v4p1.txt',l_wmap,llo,lhi,dl_wmap,ddl_wmap,ddl_noise,ddl_sv
;     wh_wmap = where(l_wmap le 8000. and l_wmap gt 2)
;     errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap+ddl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color
;     errplot,l_wmap[wh_wmap]/sf,(dl_wmap-ddl_wmap)[wh_wmap],(dl_wmap)[wh_wmap],thick=wmap_thick,color=wmap_color,width=0.008 ;centers
;     oplot,[0,0],[1d-20,1d20],color=!black,thick=1
;     oplot,[-1e9,1e9],[1,1]*40.,color=!black,thick=1
    
;     errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
; ;pause
;     oplot,l_cmbo/sf,dl_cmbo,color=color_cmb,thick=thick1,lines=2
;     oplot,l_cmbo/sf,dl_cmbo,thick=thick1,lines=2
    

;     oplot,l_cmbo/sf,dl_totalo,color=color_total,thick=thick2
;     errplot,l/sf,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=thick
    
;     xyouts,500./sf,3700,'WMAP7',color=wmap_color,chars=chars*1.5,font=1
;     xyouts,1600./sf,1100,'SPT, Full Survey',color=color_data,chars=chars*1.5,font=1
    
;     closeps
;     spawn,'convert '+outname+' '+outname_pdf
    
; endif

; if keyword_set(dosave) then begin
;     save,l_total,dl_total,dl_cmb,dl_poisson,dl_sz,dl_flat,filename='/home/kstory/lps12/scripts/plotting/bestfit_WMAP7_lps12.sav'
; endif


if keyword_set(stopit) then stop
END




;;;
; plot wmap-only best-fit minus CMB best-fit
PRO plot_bestfit_subtract

readcol, '/home/kstory/lps12/best_fit/lcdm_WMAP_ML_dl_fromRK.txt', l_wmap, dl_wmap, dl_wmap_total
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', l_cmb, dl_cmb_total


; plot, l_cmb, (dl_cmb_total - dl_wmap_total) / dl_cmb_total, $
;   ytitle = '(dl_cmb_total - dl_wmap_total) / dl_cmb_total', xtitle='ell'
; err=tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/bestfit_subtract_0823')

plot, l_cmb, (dl_cmb_total - dl_wmap_total), $
  ytitle = '(dl_cmb_total - dl_wmap_total)', xtitle='ell'
err=tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/bestfit_subtract_0824')

stop
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Plot best-fit wmap spectrum and wmap+spt spectrum on bandpowers
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO plot_bestfit_2, doeps=doeps, dosave=dosave, stopit=stopit, $
                    data=data
pdir = '/home/kstory/lps12/scripts/plotting/'

if n_elements(data) eq 0 then $
  data = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'

readcol, '/home/kstory/lps12/best_fit/lcdm_WMAP_ML_dl_fromRK.txt', l_wmap, dl_wmap, dl_wmap_total

; Get ML values, from params_c4_lcdm_camb_w7s12_ML.ini
poisson3000 = 20.0976
clust_3000 = 5.46899
sz3000 = 4.18801

;---------------------
; Get the best-fit spectrum including foregrounds
;---------------------
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
color_wmap = !red
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
oplot,l_total,dl_flat,color=color_flat,thick=thick1,linestyle=2
oplot,l_total,dl_total,color=color_total,thick=thick2
oplot,l_wmap,dl_wmap_total,color=color_wmap,thick=thick2,linestyle=2
errplot,l,dl_all*cal[wh_min]-diag_nobeam,dl_all*cal[wh_min]+diag_nobeam,color=color_data,thick=2
legend,['ML wmap','ML wmap+spt'],color=[color_wmap, color_total],pos=[1500,4000],linestyle=[2,0]

; rr = 10.^(ratio-1) * 1000.
; err_rr = 10.^(err_ratio-1.) * 1000.
; oplot,l,rr
; errplot,l,(rr-err_rr),(rr+err_rr)
;stop

;window, 2
;yr=[-1,1]*0.2
rr = 10^(ratio-1.)
err_minus = 10^(ratio-1-err_ratio)
err_plus = 10^(ratio-1+err_ratio)
;plot,l,rr*1000,/yl,yr=10^(3+yr),/yst
oplot,l,rr*1000
errplot,l,err_minus*1000, err_plus*1000
;stop

window,1
; plot residuals
plot,l,ratio,/yn,/nodata,xtitle='!12l!X!N',ytitle='Data/Theory',yr=[-1,1]*0.2+1,/yst,xr=[0,3000]
oplot,[-1,1]*9999999999.,[1,1],thick=2
oplot,l,ratio
errplot,l,(ratio-err_ratio),(ratio+err_ratio)

; plot theory line
dd = alog(dl_total) 
dd1 = dd - min(dd)
mxmn = maxmin(dd1)
dvec = (dd1 / mxmn[0]) * (1.2-0.8) + 0.8

dd_total = alog(dl_total) & dd_total = dd_total - min(dd)
dvec_total = (dd_total / mxmn[0]) * (1.2-0.8) + 0.8

dd_cmb = alog(dl_cmb) & dd_cmb = dd_cmb - min(dd)
dvec_cmb = (dd_cmb / mxmn[0]) * (1.2-0.8) + 0.8

dd_wmap = alog(dl_wmap_total) & dd_wmap = dd_wmap - min(dd)
dvec_wmap = (dd_wmap / mxmn[0]) * (1.2-0.8) + 0.8

oplot, l_total, dvec_total
;oplot, l_total, dvec_cmb, color=!green
oplot, l_wmap, dvec_wmap, color=!red, linestyle=2


; plot subtraction
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', l_ml, dl_ml_total
; oplot, l_ml, (dl_ml_total - dl_wmap_total) / dl_ml_total + 1, color=!green
xx = (dl_ml_total - dl_wmap_total)


oplot, l_ml, (dl_ml_total - dl_wmap_total) / dl_ml_total, color=!green

legend,['ML wmap','ML wmap+spt','ml_cmb-ml_wmap'],color=[color_wmap, color_total,!green],pos=[1200,1.15],linestyle=[2,0,0]

err=tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/w7s12_ps_resid_0824')
err=tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/w7s12_ps_resid2_0824')

stop
endif

if keyword_set(stopit) then stop
END

;;;
; NAME: script_0616_f2f
; PURPOSE:
;   Make plots for the face-to-face
;
; NOTES:
;
; MODIFICATION HISTORY:
;  06/16/2012: (KTS) Created
;;;


; Single map plot
PRO single_map
d = krf('/home/kstory/lps12/maps/20120420/ra22h30dec-55/map_ra22h30dec-55_150_20110411_131255.fits')
tv_spt_map, d.map.map, scale=1, /forcesize, psfile='/home/kstory/lps12/figs/map_0616.ps'

END


; Power spectrum plots
PRO ps_plots

; run_07
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_144542_kweight.sav'
dl_07 = dl_all
diag07 = diag_nobeam

; One Field
restore, '/home/kstory/lps12/end2end/end_ra22h30dec-55_07_kweight.sav'
dl_1 = reform(spectrum)
cov = reform(cov)
diag1 = dblarr(58)
for i=0, 57 do diag1[i] = sqrt(cov[i,i])

; 0809
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_172413_kweight_0809.sav'
dl_0809 = dl_all
diag0809 = diag_nobeam

; cl_theory
readcol,'/home/kstory/lps12/theory_spectrum/Cls_100sims_ave_from_combined_map.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)

; Plots
;!p.color = !black
;!p.background = !white
xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=3;1.8
whplot = indgen(45) + 10
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

; Full 2500 spectrum
window, 0, xsize=1000, ysize=700
plot,l[whplot],dl_07[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  chars=chars, charthick=2,$
  xtitle=xtitle,ytitle=ytitle,title='Full 2500 deg^2, run_07'
oplot,l_vec,dl_th,color=!red, thick=2
errplot,l[whplot],(dl_07-diag07)[whplot]*1d12,(dl_07+diag07)[whplot]*1d12, thick=2
;err=tvread(/png,/nodialog,filename=fdir+'ps_2500_0617')

; Full 2500 spectrum, no theory
window, 10, xsize=1000, ysize=700
plot,l[whplot],dl_07[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  chars=chars, charthick=2,$
  xtitle=xtitle,ytitle=ytitle,title='Full 2500 deg^2, run_07'
;oplot,l_vec,dl_th,color=!red, thick=2
errplot,l[whplot],(dl_07-diag07)[whplot]*1d12,(dl_07+diag07)[whplot]*1d12, thick=2
err=tvread(/png,/nodialog,filename=fdir+'ps_2500_noth_0617')

; Single field 
window, 1, xsize=1000, ysize=700
plot,l[whplot],dl_1[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  chars=chars, charthick=2,$
  xtitle=xtitle,ytitle=ytitle,title='PS, ra22h30dec-55'
oplot,l_vec,dl_th,color=!red, thick=2
errplot,l[whplot],(dl_1-diag1)[whplot]*1d12,(dl_1+diag1)[whplot]*1d12, thick=2
;err=tvread(/png,/nodialog,filename=fdir+'ps_ra22h30dec-55_0617')

; 0809
window, 2, xsize=1000, ysize=700
plot,l[whplot],dl_0809[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  chars=chars, charthick=2,$
  xtitle=xtitle,ytitle=ytitle,title='PS, 2008-09 fields'
oplot,l_vec,dl_th,color=!red, thick=2
errplot,l[whplot],(dl_0809-diag0809)[whplot]*1d12,(dl_0809+diag0809)[whplot]*1d12, thick=2
;err=tvread(/png,/nodialog,filename=fdir+'ps_0809_0617')

END





; errors
PRO err_plots

; run_07
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_144542_kweight.sav'
dl_07 = dl_all
diag07 = diag_nobeam
diag07_wbeam_wcal = diag

; get diags
diag_data_07 = dblarr(58)
diag_data_wbeam_07 = dblarr(58)
diag_mc_07 = dblarr(58)
diag_sv_07 = dblarr(58)
diag07_wbeam = dblarr(58)

for i=0, 57 do begin
    diag_data_07[i] = sqrt(cov_all_nobeam_data[i,i])
    diag_data_wbeam_07[i] = sqrt(cov_all_data[i,i])
    diag_mc_07[i] = sqrt(cov_all_nobeam_mc[i,i])
    diag_sv_07[i] = sqrt(cov_sv[i,i])
    diag07_wbeam[i] = sqrt( cov_all_nocal[i,i])
endfor

; ; One Field
; restore, '/home/kstory/lps12/end2end/end_ra22h30dec-55_07_kweight.sav'
; dl_1 = reform(spectrum)
; cov = reform(cov)
; diag1 = dblarr(58)
; for i=0, 57 do diag1[i] = sqrt(cov[i,i])

; ; 0809
; restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_172413_kweight_0809.sav'
; dl_0809 = dl_all
; diag0809 = diag_nobeam


; Plots
!p.color = !black
!p.background = !white
xr = [500,3e3]
yr = [1,1e2]
xtitle= '!12l!X!N'
ytitle= textoidl('\delta')+'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=3;1.8
whplot = indgen(45) + 10
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

; Full 2500 spectrum
window, 0, xsize=1000, ysize=700
plot,l[whplot],diag07[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,/xst, $
  chars=chars, charthick=2, thick=3, $
  xtitle=xtitle,ytitle=ytitle,title='errors, Full 2500 deg^2'
oplot, l[whplot],diag_sv_07[whplot]*1d12, thick=3, color=!red
oplot, l[whplot],diag_data_07[whplot]*1d12, thick=3, color=!blue
legend_str = ['err_tot', 'err_sv', 'err_noise']
legend, legend_str, linestyle=0, thick=3, chars=chars, charthick=charthick, colors=[!black, !red, !blue], position=[1900, 80]
;err=tvread(/png,/nodialog,filename=fdir+'err_2500_0617')


; Full 2500 spectrum, with beams
window, 1, xsize=1000, ysize=700
plot,l[whplot],diag07[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,/xst, $
  chars=chars, charthick=2, thick=3, $
  xtitle=xtitle,ytitle=ytitle,title='errors, Full 2500 deg^2'
oplot, l[whplot],diag07_wbeam[whplot]*1d12, thick=3, linestyle=3
oplot, l[whplot],diag_sv_07[whplot]*1d12, thick=3, color=!red
oplot, l[whplot],diag_data_07[whplot]*1d12, thick=3, color=!blue
oplot, l[whplot],diag_data_wbeam_07[whplot]*1d12, thick=3, color=!purple
oplot, l[whplot],diag07_wbeam_wcal[whplot]*1d12, thick=3, color=!green, linestyle=3
legend_str = ['err_tot', 'err_sv', 'err_noise', 'err_noise_with_beam', 'err_tot_with_beam', 'err_tot_cal&beam']
legend, legend_str, linestyle=[0,0,0,0,3,3], thick=3, chars=chars, charthick=charthick, colors=[!black, !red, !blue, !purple, !black, !green], position=[1500, 80]
;err=tvread(/png,/nodialog,filename=fdir+'err_wbeam_2500_0617')

; Plot beam errors
window, 2, xsize=1000, ysize=700
ytitle2= textoidl('\delta')+'D!D!12l!X!N'+textoidl('_{, nobeam}') + '/ ' + $
  textoidl('\delta')+'D!D!12l!X!N'+textoidl('_{, total}')
plot,l[whplot],(diag07/diag07_wbeam)[whplot],xr=xr,yr=[0.9,1],/yst,/xst, $
  chars=chars, charthick=2, thick=3, $
  xtitle=xtitle,ytitle=ytitle2,title='Beam Errors, Full 2500 deg^2'
;err=tvread(/png,/nodialog,filename=fdir+'beam_frac_err_2500_0617')

stop
END


PRO plot_0809
; plotting setup
if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
whplot = indgen(45) + 10

; run_07
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_144542_kweight.sav'
dl_07 = dl_all
diag07 = diag_nobeam

; run_07; 0809 only
restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_205039_kweight_0809.sav' ; from comb_0809
dl0_0809     = dl_all
diag0_0809 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
dl_rk = dl_all
err_rk = diag_nobeam

; Plots
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

window, 2
plot, l[whplot], ((err_rk )/diag0_0809)[whplot], thick=3, chars=1.8,charthick=2,$
  xtitle=xtitle, ytitle='err_k11/err_0809', title='Error bars, comb_0809'
oplot, [0,10000], [1,1], linestyle=1, thick=3
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v0809_0617')

window, 3
plot, l[whplot], ((err_rk )/diag07)[whplot], thick=3, chars=1.8,charthick=2,$
  xtitle=xtitle, ytitle='err_k11/err_07', title='Error bars, K11 v.s. 2500'
oplot, [0,10000], [sqrt(2500./790.),sqrt(2500./790.)], linestyle=1, thick=3
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v07_0617')

; window, 4
; plot, l[whplot], (dl0_0809/dl_rk)[whplot], thick=3, chars=1.8,charthick=2,$
;   xtitle=xtitle, ytitle='err_k11/err_0809', title='Error bars, comb_0809'
; oplot, [0,10000], [1,1], linestyle=1, thick=3
; xx = dl0_0809/dl_rk
; errplot,l[whplot],(xx-diag0_0809/dl_rk)[whplot],(xx+diag0_0809/dl_rk)[whplot]
; ;err=tvread(/png,/nodialog,filename=fdir+'err_k11v0809_0617')

; Ratio of bandpowers from k11 to 0809
window, 5, xsize=800, ysize=500
ytitle5= 'D!D!12l!X!N'+textoidl('_{, 0809}') + '/ ' + 'D!D!12l!X!N'+textoidl('_{, k11}')
plot, l[whplot], (dl0_0809/dl_rk)[whplot], thick=3, chars=2,charthick=2,yr=[0.9, 1.1],/yst,$
  xtitle=xtitle, ytitle=ytitle5, title='Ratio of bandpowers'
oplot, [0,10000], [1,1], linestyle=1, thick=3
xx = dl0_0809/dl_rk
errplot,l[14],(xx-diag0_0809/dl0_0809)[14],(xx+diag0_0809/dl_rk)[14], thick=3
err=tvread(/png,/nodialog,filename=fdir+'ps_k11v0809_0617')

stop
END




;; Kweight input plot
PRO kweight_inputs
field_idx = 3
if n_elements(l_fwhm) eq 0 then l_fwhm=1000.
if n_elements(calib) eq 0 then calib = 0.76 ; use 2008 value

; directories
out_dir = '/home/kstory/lps12/twod_kweights/'
tf_dir = '/home/kstory/lps12/twod_tfs/'
noise_dir = '/home/kstory/lps12/noise_psds/'

;------------------------
; Setup
;------------------------
print, 'MAKE_TWOD_KWEIGHTS_LPS12: setup'

f = lps12_fieldstruct()
field_name = f[field_idx].name
info = get_lps12_fieldinfo(field_idx)
nbig = info.nbig

reso = 1. ; arcminutes per pixel in data maps

; get the l-grid
l = make_fft_grid(reso/60.*!dtor,nbig,nbig)*2.*!pi
lbinsize = l[1] - l[0]
nbin = round(5000./lbinsize)
lbin = findgen(nbin)*lbinsize + 0.

; precompute 2d indices for the l-bins
for j=0,nbin-2 do begin
    ex=execute('wh'+strtrim(floor(j),2)+'=where(l ge lbin[j] and l lt lbin[j+1],nwh)')
    ex=execute('nwh'+strtrim(floor(j),2)+'=nwh')
endfor

;------------------------
; Get Inputs
;------------------------
print, 'MAKE_TWOD_KWEIGHTS_LPS12: get inputs'

; Get Theoretical Cl's
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
cl_th = interpol(cl_uK2, l_vec, l)

; get noise psd
restore, noise_dir+'noise_psd_'+field_name+'.sav'
; psd is in units of K-rad
psd_uK = psd * 1d6
noise = ( psd_uK )^2


; get the TF for this field
restore, tf_dir+'tf_'+field_name+'.sav'
tf = TF_W_BEAM

;----------------------
; calculate the weight_2d
;----------------------
print, 'MAKE_TWOD_KWEIGHTS_LPS12: make weight_2d'

; find where tf is not NaN, and tf ne 0 (can't divide by 0)
whgood = where( (finite(tf) ne 0) and (tf gt 0), nwh) 
if (nwh ne 0) then begin
    ntmp = dblarr(nbig, nbig)
    ntmp[whgood] = noise[whgood] / tf[whgood]
endif

; plot the inputs
!p.color = !black
!p.background = !white
chars=3 & charthick=2
xtitle= '!12l!X!N'
vec = indgen(1000)

window, 5, xsize=900, ysize=600
plot, l[vec,vec], cl_th[vec,vec], /ylog, $
  chars=chars, charthick=charthick, thick=2, $
  xtitle=xtitle, ytitle='Cl [uK-rad]^2', title=field_name, xra=[0,4000]
oplot, l[vec,vec], noise[vec,vec], thick=2, color=!red
oplot, l[vec,vec], ntmp[vec,vec], thick=2, color=!orange
legend_str = ['Cl_th', 'noise_psd', 'noise_TF']
legend, legend_str, linestyle=0, thick=2, chars=chars, charthick=charthick, colors=[!black, !red, !orange], position=[2000, 0.7e4]

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
err=tvread(/png,/nodialog,filename=fdir+'kweight_inputs_0617')

stop
END




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Jackknife plots
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO chisq_table
arr = [ $
0.1393, 0.5334, 0.5523, 0.0822, 0.9552, $
0.9650, 0.5555, 0.6717, 0.9665, 0.0224, $
0.8381, 0.1755, 0.4908, 0.0466, 0.6786, $
0.2745, 0.8699, 0.6104, 0.2275, 0.3083, $
0.5485, 0.6187, 0.3400, 0.0835, 0.3568, 0.2040, $
0.5012, 0.4659, 0.7422, 0.1672, 0.5017, $
0.0044, 0.4418, 0.2814, 0.6283, 0.2047, 0.3848, $
0.7700, 0.4719, 0.3487, 0.3230, 0.9465, $
0.0170, 0.5379, 0.4829, 0.8302, 0.8606, $
0.7039, 0.0776, 0.0528, 0.6015, 0.4677, $
0.0334, 0.8418, 0.4710, 0.4999, 0.5956, $
0.1021, 0.0184, 0.1247, 0.9541, 0.5770, 0.7648, $
0.5828, 0.8485, 0.5769, 0.5716, 0.5398, $
0.8299, 0.7155, 0.6691, 0.5063, 0.4251, $
0.7017, 0.4197, 0.9999, 0.9422, 0.3485, $
0.4593, 0.4407, 0.0166, 0.4185, 0.5270, $
0.1974, 0.4574, 0.8304, 0.6733, 0.0742, $
0.4764, 0.2211, 0.0776, 0.0393, 0.4657, 0.6201, $
0.8674, 0.7632, 0.8993, 0.4114, 0.4907, 0.1756, $
0.2829, 0.8023, 0.6627, 0.6754, 0.7462, 0.9983 ]

hh = histogram(arr, nbins=22, locations=xx, min=-0.05, max=1.05)
plot, xx, hh, psym=10, thick=2, charthick=2, xr=[-0.05,1.05],/xst,$
  xtitle='PTE', ytitle='num. jacks', title='PTE distribution for all individual-field jacks'

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
err=tvread(/png,/nodialog,filename=fdir+'pte_run05_0617')
stop
END




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot azrms 10% cuts
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_azrms


files = ['/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav', $
         '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_145116_kweight_azrms_90.sav', $
         '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_145300_kweight_randcut_90_seed1.sav', $
         '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_145423_kweight_randcut_90_seed2.sav', $
         '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_150600_kweight_randcut_90_seed3.sav']
nfiles = 5
nbins = 58

; make arrays
dls   = dblarr(nfiles, nbins)
diags = dblarr(nfiles, nbins)

for ii=0, nfiles-1 do begin
    restore, files[ii]
    dls[ii,*] = dl_all
    diags[ii,*] = diag_nobeam
endfor


v1 = indgen(50)+8
fdir = '/home/kstory/public_html/notebook/spt_lps12/'
color_arr = [!white, !red, !blue, !lavender, !green]

window, 1, xsize=700, ysize=500
plot, l[v1], (dls[0,*])[v1], /ylog, xtitle='ell', ytitle='dDl [K^2]', title='Randcut_90 Spectra', thick=2,charthick=2,chars=1.8
for ii=0, nfiles-1 do begin
    oplot, l[v1], (dls[ii,*])[v1], color=color_arr[ii], thick=2
endfor
legend, ['full_set', 'azrms_90', 'rand_1', 'rand_2', 'rand_3'], linestyle=[0,0,0,0,0], thick=2,charthick=2,$
  colors=[!white, !red, !skyblue, !lavender, !green], pos=[2000, 5.e-9]
;err = tvread(/png,/nodialog,filename=fdir+'azrms_ps_0617')

window, 2
plot, l[v1], ((dls[0,*] - dls[1,*])/dls[0,*])[v1], thick=2,charthick=2,chars=1.8,$
  xtitle='ell', ytitle='(dl_full - dl_new) / dl_full', title='Fraction change in dls'
for jj=1, nfiles-1 do begin
    oplot, l[v1], ((dls[0,*] - dls[jj,*])/dls[0,*])[v1], color=color_arr[jj], thick=2
endfor
legend, ['azrms_90', 'rand_1', 'rand_2', 'rand_3'], linestyle=[0,0,0,0], thick=2,charthick=2,$
  colors=[!red, !skyblue, !lavender, !green], pos=[600, 0.025]
;err = tvread(/png,/nodialog,filename=fdir+'azrms_ps_frac_0617')

window, 3
plot, l[v1], ((diags[0,*] - diags[1,*])/diags[0,*])[v1], thick=2,charthick=2,chars=1.8,$
  xtitle='ell', ytitle='(err_full - err_new) / err_full', title='Fraction change in errors'
for jj=1, nfiles-1 do begin
    oplot, l[v1], ((diags[0,*] - diags[jj,*])/diags[0,*])[v1], color=color_arr[jj], thick=2
endfor
legend, ['azrms_90', 'rand_1', 'rand_2', 'rand_3'], linestyle=[0,0,0,0], thick=2,charthick=2,$
  colors=[!red, !skyblue, !lavender, !green], pos=[2200, 0.14]
err = tvread(/png,/nodialog,filename=fdir+'azrms_ps_err_frac_0617')

stop
END



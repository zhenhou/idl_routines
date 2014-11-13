;;;
; NAME: script12_0827
; PURPOSE:
;   General script for today (2012)
;
; NOTES:
; 1) Calculate TF for pipe2, spec1 
; 2) run_09
; 3) beam error in sims
;;;


;;;;;;;;;;;;;;;;;;
;
; Get a calibration factor, given data, theory, and icov
;
;;;;;;;;;;;;;;;;;;
FUNCTION get_theory, l_tot=l_tot, dl_tot=dl_tot
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

nell = 15000
wh_f = where(f_ell eq sim_ell[0])
wh_ksz = where(ksz_ell eq sim_ell[0])

l_tot = sim_ell
dl_tot = sim_dl[0:nell] + f_dl[wh_f:wh_f+nell] + ksz_dl[wh_ksz:wh_ksz+nell]
;stop
RETURN, 1
END


;;;;;;;;;;;;;;;;;;;;;;;
;
; 2) Pipeline test2, all
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2_all
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'

; get spec1
restore, edir+'end_ra3h30dec-60_08p2_spec1_kweight_calib.sav'
calib = 1.
dl_p1 = spectrum*(1d12)^2. * (calib)
cov_p1 = cov

; get spec2
restore, edir+'end_ra3h30dec-60_08p2_spec2_kweight_calib.sav'
calib = 1.
dl_p2 = spectrum*(1d12)^2. * (calib)
cov_p2 = cov

; get final combined
restore, edir+'run_08/combined_spectrum_20120717_174249_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

;----------------------------
;Get the theory spectrum 1
dl_th1 = get_theory(l_th1b, dl_th1b)
wh = where(l_th1b ge 50 and l_th1b lt 3998)
l_th1b = l_th1b[wh]
dl_th1b = dl_th1[wh]
dl_th1b_bin = reform(dl_th1b # wf2use)

;----------------------------
;Get the theory spectrum 2
;readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec2,dl_th2
;dl_th2b = dl_th2[50:3998]
;l_vec2b  = l_vec2[50:3998]
dl_th2b = shift(dl_th1b, -10) ; shift by 10:
dl_th2b_bin = dl_th2b # wf2use



;----------------------------
; Plotting
;----------------------------

vec=indgen(47) + 9

window, 1
plot, l[vec], dl_p1[vec], xr=[500,3000], yr=[10, 5000], ystyle=1, /yl

; spec1
oplot, l[vec], dl_p1[vec], color=!red, linestyle=0,thick=2

; spec2
oplot, l[vec], dl_p2[vec], color=!darkgreen, linestyle=0,thick=2

;oplot, l[vec], dl_all[vec]

oplot, l[vec], dl_th1b[vec], color=!purple, linestyle=2,thick=2
oplot, l[vec], dl_th2b[vec], color=!blue, linestyle=2,thick=2

legend,['spec1','spec2','th1','th2'],color=[!red,!darkgreen,!purple,!blue],linestyle=[0,0,2,2],$
  pos=[1800,4000]

ff = '/home/kstory/public_html/notebook/spt_lps12/pipe2_all_0827'
;err=tvread(/png,/nodialog,filename=ff)
stop
END







;;;;;;;;;;;;;;;;;;;;
; 
; run_09
;
;;;;;;;;;;;;;;;;;;;;

PRO end2end_run09_a
tbegin = systime(0, /seconds)
for idx=0, 4 do begin
    print, 'run_end2end_09, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_09, idx, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END


PRO end2end_run09_b
tbegin = systime(0, /seconds)
for idx=5, 9 do begin
    print, 'run_end2end_09, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_09, idx, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END


PRO end2end_run09_c
tbegin = systime(0, /seconds)
for idx=10, 14 do begin
    print, 'run_end2end_09, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_09, idx, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END


PRO end2end_run09_d
tbegin = systime(0, /seconds)
for idx=15, 19 do begin
    print, 'run_end2end_09, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_09, idx, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END






;;;;;;;;;;;;;;;;;;;;;;;
;
; 3) Beam error in sims
;
;;;;;;;;;;;;;;;;;;;;;;;
PRO frac_plots

ybase = 2008
for ii=0, 3 do begin
    year = ybase+ii
    ystr = strtrim(string(year),2)
    bb = get_lps12_beams(year,l,bl,/s09)
    bbx = get_lps12_beams(year,lx,blx,/sim_lmax8000)

    bl2 = bl^2.
    blx2 = blx^2.

    window, ii
    plot, l, (blx2-bl2)/bl2, title=ystr, ytitle='(bl_x^2 - bl_true^2)/bl_true^2'
endfor


fdir='/home/kstory/public_html/notebook/spt_lps12/'

wset, 0
err=tvread(/png,/nodialog,filename=fdir+'frac_beam_err_2008_0827')
wset, 1
err=tvread(/png,/nodialog,filename=fdir+'frac_beam_err_2009_0827')


stop
END



PRO beam_err
comb_file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
restore, comb_file
beam_err_scalefactor_in = 1.0
beamdir = '/data/sptdat/beams/rev3.1/'
yy=[2008, 2009, 2010, 2011]
beam_cor = get_beam_correlation_matrix(l=l, year=yy, weight_years=wyear, band1=150, band2=150, scalefactor=beam_err_scalefactor, beamdir=beamdir)

; get eigen vectors and values
evecs = double(beam_cor)
trired,evecs,evals,e
triql,evals,e,evecs

delta_l = 50. & min_l = 250. & max_l = 3150. & nl = floor((max_l-min_l)/delta_l) & lhi = dindgen(nl)*delta_l + min_l
banddef = lhi
nbands=n_elements(banddef)

; add calibration uncertainty and beam error
beamcal_cor = beam_cor + (cal_error_power)^2
bccov = beam_cor*0.
for i=0, nbands-1 do begin
    for j=0, nbands-1 do begin
        bccov[i,j] = beamcal_cor[i,j]*dl_all[i]*dl_all[j]
    endfor
endfor

; Same thing, but for beam-only
bcov = beam_cor*0.
for i=0, nbands-1 do begin
    for j=0, nbands-1 do begin
        bcov[i,j] = beam_cor[i,j]*dl_all[i]*dl_all[j]
    endfor
endfor


bccov_d = dblarr(58)
bcov_d = dblarr(58)
for i=0, 57 do begin
    bccov_d[i] = sqrt(bccov[i,i])
    bcov_d[i] = sqrt(bcov[i,i])
endfor

; PLOT
fdir='/home/kstory/public_html/notebook/spt_lps12/'

plot, banddef, bccov_d/dl_all,xtitle='ell',title='Beam+Cal error'
oplot, banddef, bcov_d/dl_all, color=!red
legend,['bccov_d/dl_all','bcov_d/dl_all'],colors=[!black,!red],linestyle=[0,0]

err=tvread(/png,/nodialog,filename=fdir+'beam_cal_err_0827')
stop
END



PRO compare_spec, idx
idx=5
edir='/home/kstory/lps12/end2end/'
f=lps12_fieldstruct()

restore, edir+'end_'+f[idx].name+'_08_kweight.sav'
ps08 = spectrum*1d12

restore, edir+'end_'+f[idx].name+'_09_kweight.sav'
ps09 = spectrum*1d12

vec=indgen(47)+8

fdir='/home/kstory/public_html/notebook/spt_lps12/'
; PLOTS
window, 1
plot, banddef[vec], ps08[vec], /yl, ytitle='Dl [uK^2]', title=f[idx].name
oplot, banddef[vec], ps09[vec], color=!red
legend,['ps run_08','ps run_09'],colors=[!black,!red],linestyle=[0,0]
err=tvread(/png,/nodialog,filename=fdir+'ps_idx5_0827')


window, 2
plot, banddef[vec], ((ps09-ps08)/ps09)[vec], yr=[-0.01,0], ytitle='(ps09-ps08)/ps09', title=f[idx].name
err=tvread(/png,/nodialog,filename=fdir+'psFrac_idx5_0827')

stop
END



PRO compare_spec_15
idx=15
edir='/home/kstory/lps12/end2end/'
f=lps12_fieldstruct()

restore, edir+'end_'+f[idx].name+'_08_kweight.sav'
ps08 = spectrum*1d12

restore, edir+'end_'+f[idx].name+'_09_kweight.sav'
ps09 = spectrum*1d12

vec=indgen(47)+8

fdir='/home/kstory/public_html/notebook/spt_lps12/'
; PLOTS
window, 1
plot, banddef[vec], ps08[vec], /yl, ytitle='Dl [uK^2]', title=f[idx].name
oplot, banddef[vec], ps09[vec], color=!red
legend,['ps run_08','ps run_09'],colors=[!black,!red],linestyle=[0,0]
err=tvread(/png,/nodialog,filename=fdir+'ps_idx15_0827')


window, 2
plot, banddef[vec], ((ps09-ps08)/ps09)[vec], yr=[-0.01,0], ytitle='(ps09-ps08)/ps09', title=f[idx].name
err=tvread(/png,/nodialog,filename=fdir+'psFrac_idx15_0827')

stop
END




;;;
; NAME: script_0507
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) coadd remaining fields (17,18,19)
;
; MODIFICATION HISTORY:
;  05/07/2012: (KTS) Created
;;;

; Process a field
PRO process_field, idx
print, 'coadd maps, idx = ', idx
coadd_maps_0420, idx

print, 'make mask, idx = ', idx
make_mask_lps12, field_idx=idx

;print, 'make noise psd, idx = ', idx
;make_noise_psd_lps12, idx
;make_coupling_kernel, idx
;try_end2end, idx, run='03', /use_kweight, /resume
;make_twod_tfs, idx ??
;make_twod_kweights_lps12, idx ??
END

;;;;;;;;;;;;;;;;;;;;;;;;
; Processing for today
;;;;;;;;;;;;;;;;;;;;;;;;


PRO kweight_3 ; DONE
make_twod_kweights_lps12, 3, /save_files
END

PRO end_1 ; DONE
try_end2end, 1, run='03', /use_kweight, /resume
END

PRO process_2 ; running (9:30pm)
idx = 2
time_0 = systime(0,/seconds)

;print, 'Run try_end2end, idx=', idx
;try_end2end, idx, run='03', /resume
time_1 = systime(0,/seconds)

; print, 'Run tf_prep, idx=', idx
; tf_prep, idx, '03'
time_2 = systime(0,/seconds)

print, 'Run make_twod_tfs, idx=', idx
make_twod_tfs, idx, '03', /dosave
time_3 = systime(0,/seconds)

print, 'Run make_twod_kweights, idx=', idx
make_twod_kweights_lps12, idx, /save_files
time_4 = systime(0,/seconds)

print, 'Re-Run try_end2end, idx=', idx
try_end2end, idx, run='03', /use_kweight, /resume
time_5 = systime(0,/seconds)

print, 'Run-times [min]:'
print, 'try_end2end, 1st go: ', (time_1-time_0) / 60.
print, 'save_mode_coupling_lps12: ', (time_2-time_1) / 60.
print, 'make_twod_tfs: ', (time_3-time_2) / 60.
print, 'make_twod_kweights: ', (time_4-time_3) / 60.
print, 'try_end2end, 2nd go: ', (time_5-time_4) / 60.

END

PRO tf_prep_6789 ; done
tf_prep, 6, '03'
tf_prep, 7, '03'
tf_prep, 8, '03'
tf_prep, 9, '03'
END

PRO tf_prep_11121314 ; done
tf_prep, 11, '03'
tf_prep, 12, '03'
tf_prep, 13, '03'
tf_prep, 14, '03'
END

PRO tf_prep_9101617 ; running
tf_prep, 9, '03'
tf_prep, 10, '03'
tf_prep, 16, '03'
tf_prep, 17, '03'
END

PRO make_tf_6789 ; done
make_twod_tfs, 6, '03', /dosave
make_twod_tfs, 7, '03', /dosave
make_twod_tfs, 8, '03', /dosave
make_twod_tfs, 9, '03', /dosave
END

PRO make_tf_10to14 ; done
for ii=10, 14 do begin
    make_twod_tfs, ii, '03', /dosave
endfor
END

PRO make_tf_1617 ; done
make_twod_tfs, 16, '03', /dosave
make_twod_tfs, 17, '03', /dosave
END

PRO kweight_many ; DONE
for ii=6, 14 do begin
    make_twod_kweights_lps12, ii, /save_files
endfor

make_twod_kweights_lps12, 16, /save_files
make_twod_kweights_lps12, 17, /save_files
END

PRO end_45 ; DONE
try_end2end, 4, run='03', /resume
try_end2end, 5, run='03', /resume
END

;;; end2end, with kweight
PRO endkw_3to8 ; running
for ii=3, 8 do begin
    try_end2end, ii, run='03', /resume, /use_kweight
endfor
END

PRO endkw_9to14 ; running
for ii=9, 14 do begin
    try_end2end, ii, run='03', /resume, /use_kweight
endfor
END

; re-do field 10, since something is funky
PRO fix_10 ; done, still not working
idx = 10
;make_noise_psd_lps12, idx - > already re-ran
try_end2end, idx, run='03', /resume
tf_prep, idx, '03'
make_twod_tfs, idx, '03', /dosave
make_twod_kweights_lps12, idx, /save_files
try_end2end, idx, run='03', /resume, /use_kweight
END

; run field 19
PRO run_19
idx = 19
;make_noise_psd_lps12, idx - > already re-ran
try_end2end, idx, run='03', /resume
tf_prep, idx, '03'
make_twod_tfs, idx, '03', /dosave
make_twod_kweights_lps12, idx, /save_files
try_end2end, idx, run='03', /resume, /use_kweight
END

;;;;;;;;;;;;;;;;;;;;;;;;
; jackknives
;;;;;;;;;;;;;;;;;;;;;;;;
PRO jack_lr_0to14
for ii=0, 14 do begin
    lps12_jack, ii, 'lr'
endfor
END

PRO jack_12_0to17
for ii=0, 14 do begin
    lps12_jack, ii, '12'
endfor
lps12_jack, 16, '12'
lps12_jack, 17, '12'

; finish out lr jack
lps12_jack, 16, 'lr'
lps12_jack, 17, 'lr'
END



;;;;;;;;;;;;;;;;;;;;;;;;
; band powers
;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Inputs to kweight
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO check_kweight
field_idx=3

save_files = 0
save_plots = 0
; Set the FWHM for Gaussian smoothing:
l_fwhm=1000.
calib = 0.76 ; use 2008 value

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

; correct the trasfer function beyond 4500 to avoid crazy effects
cut_idx = floor(4500./lbinsize)-2 ; index
whfixell = where(l gt l[cut_idx])
tf[whfixell] = tf[cut_idx-1]

;----------------------
; calculate the weight_2d
;----------------------
print, 'MAKE_TWOD_KWEIGHTS_LPS12: make weight_2d'

; find where tf is not NaN, and tf ne 0 (can't divide by 0)
whgood = where( (finite(tf) ne 0) and (tf ne 0), nwh) 
if (nwh ne 0) then begin
    ntmp = dblarr(nbig, nbig)
    ntmp[whgood] = noise[whgood] / tf[whgood]
    ;weight_2d = ( gauss_smooth(cl_th,l_fwhm,lbinsize) + gauss_smooth(ntmp,l_fwhm,lbinsize) )^(-2.)
    weight_2d = ( cl_th + gauss_smooth(ntmp,l_fwhm,lbinsize) )^(-2.)
endif else begin
    weight_2d = ( gauss_smooth(cl_th,l_fwhm,lbinsize) + gauss_smooth(noise/tf,l_fwhm,lbinsize) )^(-2.)
endelse

w1 = weight_2d ; Plotting only

; renormalize the weight within each l-ring, (l-bin), i.e. divide by
; the maximum weight within that bin.
for j=0,nbin-2 do begin
    ex=execute('wh=wh'+strtrim(floor(j),2))
    ex=execute('nwh=nwh'+strtrim(floor(j),2))
    if nwh gt 0 then begin
        weight_2d[wh] = weight_2d[wh] / max(weight_2d[wh])
    endif
endfor

; Re-set crazy outliers.  This only applies where l > 5000.
weight_2d[where(weight_2d gt 5)] = 5




;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plots
sdir = '/home/kstory/public_html/notebook/spt_lps12/'
; tv_spt_map, shift(noise, 2160, 2160), /norms, reso=5, winnum=10,min=1e-6, max=1e-4, xra=[-4500,4500],yra=[-4500,4500], title='Noise. '+field_name
; wset, 10 & err=tvread(/png,/nodialog, filename=sdir+'noise_0508')

; xtf = tf
; xtf[where( finite(tf) eq 0)] = 0.
; tv_spt_map, shift(xtf, 2160, 2160), /norms, reso=5, winnum=11, min=0.01, max=0.5,xra=[-4500,4500],yra=[-4500,4500], title='TF. '+field_name
; wset, 11 & err=tvread(/png,/nodialog, filename=sdir+'tf_0508')

; tv_spt_map, shift(ntmp, 2160, 2160), /norms, reso=5, winnum=13,min=1e-6,max=3e-4, xra=[-4500,4500],yra=[-4500,4500], title='Noise / TF. '+field_name
; wset, 13 & err=tvread(/png,/nodialog, filename=sdir+'noise_div_tf_0508')

; tv_spt_map, shift(gauss_smooth(ntmp,1000.,5.), 2160, 2160), /norms, reso=5, winnum=14,min=1e-6,max=3e-4, xra=[-4500,4500],yra=[-4500,4500], title='GS_1000{N/TF}. '+field_name
; wset, 14 & err=tvread(/png,/nodialog, filename=sdir+'gs1000_noise_div_tf_0508')

; tv_spt_map, shift(gauss_smooth(ntmp,300.,5.), 2160, 2160), /norms, reso=5, winnum=15,min=1e-6,max=3e-4, xra=[-4500,4500],yra=[-4500,4500], title='GS_300{N/TF}. '+field_name
; wset, 15 & err=tvread(/png,/nodialog, filename=sdir+'gs300_noise_div_tf_0508')

; stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;----------------------
; view this field's weight    
;----------------------
fig_dir = '/home/kstory/public_html/notebook/spt_lps12/'

;tv_spt_map,shift(weight_2d,nbig/2.,nbig/2.), /norms, min=0,max=1,reso=5,xr=[-3000,3000],yr=[-3000,3000],title='weight_2d=1 / ( GS{cl_th} + GS{noise/tf} )^2, '+field_name
tv_spt_map,shift(weight_2d,nbig/2.,nbig/2.), /norms, min=0,max=1,reso=5,xr=[-3000,3000],yr=[-3000,3000],title='weight_2d=1 / ( cl_th + GS{noise/tf} )^2, '+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
err = tvread(/png, filename=fig_dir+'kweight_noSmoothClth_'+field_name+'_0508', /nodialog)

stop
; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,xr=[-3000,3000],yr=[-3000,3000],title='un-normalized LOG weight_2d'+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
;err = tvread(/png, filename=fig_dir+'kweight_unnorm_'+field_name+'_0508', /nodialog)

; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,title='un-normalized LOG weight_2d'+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
;err = tvread(/png, filename=fig_dir+'kweight_unnorm_zoom_'+field_name+'_0508', /nodialog)

; inputs plot
window, 5, xsize=700, ysize=500
vec = indgen(1000)
plot, l[vec,vec], cl_th[vec,vec], /ylog, xtitle='ell', ytitle='Cl [uK-rad]^2', title=field_name, xra=[0,4000]
oplot, l[vec,vec], noise[vec,vec], color=!red
oplot, l[vec,vec], ntmp[vec,vec], color=!orange
legend_str = ['Cl_th', 'noise_psd', 'noise_TF']
legend, legend_str, linestyle=0, colors=[!white, !red, !orange], position=[2500, 0.9e2]
;err = tvread(/png, filename=fig_dir+'kweight_inputs_0508', /nodialog)

stop
END


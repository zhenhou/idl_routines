;;;
; NAME: script_0510
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) coadd remaining fields (17,18,19)
;
; MODIFICATION HISTORY:
;  05/09/2012: (KTS) Created
;;;

; Process a field
PRO process_field, idx
;coadd_maps_0420, idx
;make_mask_lps12, field_idx=idx
;make_noise_psd_lps12, idx
;make_coupling_kernel, idx
;make_noise_psd_lps12, idx

; 1> tf_prep, idx, '03'
; 2> make_twod_tfs, idx, '03', /dosave
; 3> make_twod_kweights_lps12, idx, /save_files
; 4> try_end2end, idx, run='03', /resume, /use_kweight

;time_0 = systime(0,/seconds)
END

;;;;;;;;;;;;;;;;;;;;;;;;
; Processing for today
;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;
; jackknives
;;;;;;;;;;;;;;;;;;;;;;;;
PRO jack_151819 ; DONE
lps12_jack, 15, 'lr'
lps12_jack, 18, 'lr'
lps12_jack, 19, 'lr'
lps12_jack, 15, '12'
lps12_jack, 18, '12'
lps12_jack, 19, '12'
print, systime()
END

PRO jack_tweight_all ; running
print, systime()
t0 = systime(0, /seconds)
for idx=0, 19 do begin
    lps12_jack, idx, 'tweight'
endfor

t1 = systime(0, /seconds)
t_tot = (t1 - t0) / 60.
print, 'process took ', t_tot, ' minutes.'
END

PRO jack_sim_all ; running
print, systime()
t0 = systime(0, /seconds)
for idx=2, 19 do begin
    lps12_jack, idx, 'lr_sim', /sim
endfor

t1 = systime(0, /seconds)
t_tot = (t1 - t0) / 60.
print, 'process took ', t_tot, ' minutes.'
END

; coadd lr sims
; before running this, run:
; IDL> .com /home/cr/code/spt/lowell/coadd_lr_sims.pro
PRO coadd_sims_all
f = lps12_fieldstruct()
;spawn, 'cd /home/cr/code/spt/lowell/'
for idx=3, 19 do begin
    coadd_lr_sims, f[idx].name
endfor
END

; fail...
PRO sss
RESOLVE_ROUTINE, '/home/cr/code/spt/lowell/coadd_lr_sims.pro', /IS_FUNCTION
help, /source_files
END

;;;
; script for oliver
;;;

PRO read_coaddsim, idx
f = lps12_fieldstruct()
field = f[idx].name

info = get_lps12_fieldinfo(idx) 
nx = info.npix[0] & ny = info.npix[1]
n_sim_coadds = 100
;fname = '/home/kstory/lps12/sims/coaddsim_'+field+'.dat'
fname = '/home/kstory/lps12/sims/jacklr/coaddsim_lr_'+field+'.dat'
maps = read_dat_map(fname, nx, ny, n_sim_coadds)

stop
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Inputs to kweight
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO check_kweight
field_idx=0

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
whgood = where( (finite(tf) ne 0) and (tf gt 0), nwh) 
;whgood = where( (finite(tf) ne 0) and (tf gt 0), nwh, complement=whbad) 
if (nwh ne 0) then begin
    ntmp = dblarr(nbig, nbig)
    ntmp[whgood] = noise[whgood] / tf[whgood]
    weight_2d = ( gauss_smooth(cl_th,l_fwhm,lbinsize) + gauss_smooth(ntmp,l_fwhm,lbinsize) )^(-2.)
    print, 'HERE!!!'
    ;weight_2d = ( cl_th + gauss_smooth(ntmp,l_fwhm,lbinsize) )^(-2.)
endif else begin
    weight_2d = ( gauss_smooth(cl_th,l_fwhm,lbinsize) + gauss_smooth(noise/tf,l_fwhm,lbinsize) )^(-2.)
endelse

; TESTING
;weight_2d[whbad] = 0.

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
tv_spt_map, shift(noise, 2160, 2160), /norms, reso=5, winnum=10,min=1e-6, max=1e-4, xra=[-4500,4500],yra=[-4500,4500], title='Noise. '+field_name
; wset, 10 & err=tvread(/png,/nodialog, filename=sdir+'noise_0510')

xtf = tf
xtf[where( finite(tf) eq 0)] = 0.
tv_spt_map, shift(xtf, 2160, 2160), /norms, reso=5, winnum=11, min=0.01, max=0.5,xra=[-4500,4500],yra=[-4500,4500], title='TF. '+field_name
; wset, 11 & err=tvread(/png,/nodialog, filename=sdir+'tf_0510')

tv_spt_map, shift(ntmp, 2160, 2160), /norms, reso=5, winnum=13,min=1e-6,max=3e-4, xra=[-4500,4500],yra=[-4500,4500], title='Noise / TF. '+field_name
; wset, 13 & err=tvread(/png,/nodialog, filename=sdir+'noise_div_tf_0510')

tv_spt_map, shift(gauss_smooth(ntmp,1000.,5.), 2160, 2160), /norms, reso=5, winnum=14,min=1e-6,max=3e-4, xra=[-4500,4500],yra=[-4500,4500], title='GS_1000{N/TF}. '+field_name
; wset, 14 & err=tvread(/png,/nodialog, filename=sdir+'gs1000_noise_div_tf_0510')

tv_spt_map, shift(gauss_smooth(ntmp,300.,5.), 2160, 2160), /norms, reso=5, winnum=15,min=1e-6,max=3e-4, xra=[-4500,4500],yra=[-4500,4500], title='GS_300{N/TF}. '+field_name
; wset, 15 & err=tvread(/png,/nodialog, filename=sdir+'gs300_noise_div_tf_0510')

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
;err = tvread(/png, filename=fig_dir+'kweight_noSmoothClth_'+field_name+'_0510', /nodialog)

; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,xr=[-3000,3000],yr=[-3000,3000],title='un-normalized LOG weight_2d'+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
;err = tvread(/png, filename=fig_dir+'kweight_unnorm_'+field_name+'_0510', /nodialog)

; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,title='un-normalized LOG weight_2d'+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
;err = tvread(/png, filename=fig_dir+'kweight_unnorm_zoom_'+field_name+'_0510', /nodialog)

; inputs plot
window, 5, xsize=700, ysize=500
vec = indgen(1000)
plot, l[vec,vec], cl_th[vec,vec], /ylog, xtitle='ell', ytitle='Cl [uK-rad]^2', title=field_name, xra=[0,4000]
oplot, l[vec,vec], noise[vec,vec], color=!red
oplot, l[vec,vec], ntmp[vec,vec], color=!orange
legend_str = ['Cl_th', 'noise_psd', 'noise_TF']
legend, legend_str, linestyle=0, colors=[!white, !red, !orange], position=[2500, 0.9e2]
;err = tvread(/png, filename=fig_dir+'kweight_inputs_0510', /nodialog)

stop
END


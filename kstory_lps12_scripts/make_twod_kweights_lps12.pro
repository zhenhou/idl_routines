;;;
; NAME: make_twod_kweights_lps12.pro
; PURPOSE:
;   Make weight_2ds for fields
;   seasons
;
;   weight_2d = 1 / ( GS{cl_th} + GS{noise/tf} )^2
;     where, GS{ fn } is gauss_smooth(fn, fwhm=1000,10.)
;
; CALLING SEQUENCE: 
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   l_fwhm,      Full-width half-maximum for Gaussian Smoothing kernel
;   save_files,  Option to save weight_2ds to .sav files
;   save_plots,  Option to save plots of weight_2ds
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   sav files for a weight_2d for each field from 2008 and 2009
;
; NOTES:
;   1) 
;
; MODIFICATION HISTORY:
;   11/22/2010 (RK) : Copied from Ryan's directory, when it was named make_twodim_weight.pro
;   12/14/2010: (KTS) Created make_weight_2d.pro
;   04/27/2012: (KTS) Copy to make_twod_kweights.pro, modify for lps12
;   06/05/2012: (KTS) Modify for run 05
;   06/28/2012: (KTS) Fix default calib to run_07 value (0.825)
;;;

PRO make_twod_kweights_lps12, field_idx, l_fwhm=l_fwhm, calib=calib, save_files=save_files, save_plots=save_plots, sim_lmax4500=sim_lmax4500

if n_elements(save_files) eq 0 then save_files = 0
if n_elements(save_plots) eq 0 then save_plots = 0
; Set the FWHM for Gaussian smoothing:
if n_elements(l_fwhm) eq 0 then l_fwhm=1000.
if n_elements(calib) eq 0 then calib=0.825

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
psd_uK = psd * 1d6 * calib^2.
noise = ( psd_uK )^2


; get the TF for this field
restore, tf_dir+'tf_'+field_name+'.sav'
tf = TF_W_BEAM

; correct the trasfer function beyond 4500 to avoid crazy effects
; This is not needed for sims_lmax8000
if keyword_set(sim_lmax4500) then begin
    cut_idx = floor(4500./lbinsize)-2 ; index
    whfixell = where(l gt l[cut_idx])
    tf[whfixell] = tf[cut_idx-1]
endif

;----------------------
; calculate the weight_2d
;----------------------
print, 'MAKE_TWOD_KWEIGHTS_LPS12: make weight_2d'

; find where tf is not NaN, and tf ne 0 (can't divide by 0)
whgood = where( (finite(tf) ne 0) and (tf gt 0), nwh) 
if (nwh ne 0) then begin
    ntmp = dblarr(nbig, nbig)
    ntmp[whgood] = noise[whgood] / tf[whgood]
    weight_2d = ( gauss_smooth(cl_th,l_fwhm,lbinsize) + gauss_smooth(ntmp,l_fwhm,lbinsize) )^(-2.)
endif else begin
    weight_2d = ( gauss_smooth(cl_th,l_fwhm,lbinsize) + gauss_smooth(noise/tf,l_fwhm,lbinsize) )^(-2.)
endelse

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
;weight_2d[where(weight_2d gt 5)] = 5
weight_2d[where(weight_2d gt 1)] = 1.

;----------------------
; view this field's weight    
;----------------------
tv_spt_map,shift(weight_2d,nbig/2.,nbig/2.), $
  min=0,max=1,reso=5,xr=[-3000,3000],yr=[-3000,3000],/norms, title='weight_2d=1 / ( GS{cl_th} + GS{noise/tf} )^2, '+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data

; plot the inputs
window, 5, xsize=700, ysize=500
vec = indgen(1000)
plot, l[vec,vec], cl_th[vec,vec], /ylog, xtitle='ell', ytitle='Cl [uK-rad]^2', title=field_name, xra=[0,4000]
oplot, l[vec,vec], noise[vec,vec], color=!red
oplot, l[vec,vec], ntmp[vec,vec], color=!orange
legend_str = ['Cl_th', 'noise_psd', 'noise_TF']
legend, legend_str, linestyle=0, colors=[!white, !red, !orange], position=[2500, 0.9e2]

;----------------------
; Save the weight into a .sav file
;----------------------
if (save_files) then begin
    fname=out_dir+'weight_2d_'+field_name+'.sav'
    print, 'saving file: '+fname
    save,weight_2d,filename=fname
endif

; Save the plots
if (save_plots) then begin
    plotname = 'weight_2d_'+field_name
    err = tvread(/png, filename=plotname, /nodialog)
endif

END


;;; 
; obsolete
;;;

; ; get mask
;     restore,'/data/rkeisler/ps09/mask_'+field+'_20101009_055240.sav'
; ; get bigmask (zero-padded mask)
;     bigmask = fltarr(nbig,nbig)
;     bigmask[0:n_elements(mask[*,0])-1, 0:n_elements(mask[0,*])-1] = mask

; ; get data coadd (this is the MAP (left+right), and DMAP
; ; (left-right)).  we'll use the DMAP to estimate the noise power.
;     restore,'/data/rkeisler/ps09/1.18/coadd_'+field+'.sav'

;     bigmap = fltarr(nbig,nbig)
;     bigmap[0:n_elements(map[*,0])-1, 0:n_elements(map[0,*])-1] = map

;     bigdmap = fltarr(nbig,nbig)
;     bigdmap[0:n_elements(dmap[*,0])-1, 0:n_elements(dmap[0,*])-1] = dmap

; ; get simulated coadd for this field
;     ex=execute('s=s'+strtrim(floor(i),2))
;     if s.field ne field then stop
    
;     big_map_power = abs(fft(bigmap*bigmask))^2.
;     big_dmap_power = abs(fft(bigdmap*bigmask))^2.


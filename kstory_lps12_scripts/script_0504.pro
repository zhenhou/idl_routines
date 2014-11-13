;;;
; NAME: script_0504
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) make_maps
; 2) process fields
; 3) check_kweight, make plots of kweight inputs, copy of kweight code
; 4) 
;
; MODIFICATION HISTORY:
;  05/04/2012: (KTS) Created
;;;

; make the maps
PRO make_maps_lps12, idx
f = lps12_fieldstruct()
field = f[idx].name
map_script = '/data/kstory/projects/lps12/maps/20120420/make_maps_'+field+'_150.txt'
log_file   = '/data/kstory/projects/lps12/maps/20120420/run2_'+field+'.out'
spawn, '/home/kstory/sptcdevel/mapping/mpibatch.x '+map_script+' > & '+log_file
END

; Process a field
PRO process_field, idx
print, 'coadd maps, idx = ', idx
coadd_maps_0420, idx

print, 'make mask, idx = ', idx
make_mask_lps12, field_idx=idx

;print, 'make noise psd, idx = ', idx
;make_noise_psd_lps12, idx
;make_coupling_kernel, idx
;try_end2end, idx ??
;make_twod_tfs, idx ??
;make_twod_kweights_lps12, idx ??
END

;;;;;;;;;;;;;;;;;;;;;;;;
; Processing for today
;;;;;;;;;;;;;;;;;;;;;;;;

; fill-in missing map
PRO maps_9
print, 'Make maps, idx = ', 9
make_maps_lps12, 9
END

PRO process_26
print, 'processing field ', 6
process_field, 2
print, 'processing field ', 2
process_field, 6
END

PRO process_many
process_field, 7
process_field, 8
process_field, 9
process_field, 10
process_field, 11
process_field, 12
process_field, 13
process_field, 14
process_field, 15
process_field, 16
END

PRO coupling_many
for ii=2, 9 do begin
    make_coupling_kernel, ii
endfor

for ii=11, 16 do begin
    make_coupling_kernel, ii
endfor    
END


;;;;;;;;;;;;;;;;;;;;;;;;
; band powers
;;;;;;;;;;;;;;;;;;;;;;;;

; plot band powers from end2end
PRO sss
idx = 0
run = '03'

f = lps12_fieldstruct()
field = f[idx].name
edir = '/home/kstory/lps12/end2end/'

restore, edir+'end_'+field+'_'+run+'_kweight.sav'

stop

END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Inputs to kweight
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO check_kweight

save_files = 0
save_plots = 0
; Set the FWHM for Gaussian smoothing:
l_fwhm=1000.
calib = 0.76 ; use 2008 value
field_idx=3

;------------------------
; Setup
;------------------------
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


lbin = findgen(nbin)*lbinsize + 0.

; precompute 2d indices for the l-bins
for j=0,nbin-2 do begin
    ex=execute('wh'+strtrim(floor(j),2)+'=where(l ge lbin[j] and l lt lbin[j+1],nwh)')
    ex=execute('nwh'+strtrim(floor(j),2)+'=nwh')
endfor

; get the average power spectrum from the simulated, noise-free observations.
;  Relevant variables: s#.power
;restore,'/data/rkeisler/ps09/create_twodim_tfs.sav'
;restore, '/home/kstory/lps12/scripts/make_twod_tfs.sav'

;------------------------
; Get Inputs
;------------------------

; Get Theoretical Cl's
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
cl_th = interpol(cl_uK2, l_vec, l)

; get noise psd
restore, '/data/kstory/projects/lps12/noise_psds/noise_psd_'+field_name+'.sav'
; psd is in units of K-rad
psd_uK = psd * 1d6
noise = ( psd_uK )^2


; get the TF for this field
restore, '/data/kstory/projects/lps12/twod_tfs/tf_'+field_name+'.sav'
tf = TF_W_BEAM

; correct the trasfer function beyond 4500 to avoid crazy effects
cut_idx = floor(4500./5.)-2 ; index
whfixell = where(l gt l[cut_idx])
tf[whfixell] = tf[cut_idx-1]

;----------------------
; calculate the weight_2d
;----------------------
print, 'Gathered all inputs.  Now make the weight_2d'
; find where tf is not NaN, and tf ne 0 (can't divide by 0)
whgood = where( (finite(tf) ne 0) and (tf ne 0), nwh) 
if (nwh ne 0) then begin
    ntmp = dblarr(nbig, nbig)
    ntmp[whgood] = noise[whgood] / tf[whgood]
    weight_2d = ( gauss_smooth(cl_th,l_fwhm,5.) + gauss_smooth(ntmp,l_fwhm,5.) )^(-2.)
    ;weight_2d = ( cl_th + gauss_smooth(ntmp,l_fwhm,5.) )^(-2.)
endif else begin
    weight_2d = ( gauss_smooth(cl_th,l_fwhm,5.) + gauss_smooth(noise/tf,l_fwhm,5.) )^(-2.)
    ;weight_2d = ( cl_th + gauss_smooth(noise/tf,l_fwhm,5.) )^(-2.)
endelse

w1 = weight_2d ; unnormalized

; renormalize the weight within each l-ring, (l-bin), i.e. divide by
; the maximum weight within that bin.
for j=0,nbin-2 do begin
    ex=execute('wh=wh'+strtrim(floor(j),2))
    ex=execute('nwh=nwh'+strtrim(floor(j),2))
    if nwh gt 0 then begin
        weight_2d[wh] = weight_2d[wh] / max(weight_2d[wh])
    endif
endfor

for j=0, nbin-2 do begin & ex=execute('wh=wh'+strtrim(floor(j),2)) & ex=execute('nwh=nwh'+strtrim(floor(j),2)) & if nwh gt 0 then weight_2d[wh] = weight_2d[wh] / max(weight_2d[wh]) & endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plots
sdir = '/home/kstory/public_html/notebook/spt_lps12/'
tv_spt_map, shift(noise, 2160, 2160), /norms, reso=5, winnum=10,min=1e-6, max=1e-4, xra=[-4500,4500],yra=[-4500,4500]
;wset, 10 & err=tvread(/png,/nodialog, filename=sdir+'noise_0504')

xtf = tf
xtf[where( finite(tf) eq 0)] = 0.
tv_spt_map, shift(xtf, 2160, 2160), /norms, reso=5, winnum=11, min=0.01, max=0.5,xra=[-4500,4500],yra=[-4500,4500]
;wset, 11 & err=tvread(/png,/nodialog, filename=sdir+'tf_0504')

tv_spt_map, shift(ntmp, 2160, 2160), /norms, reso=5, winnum=13,min=1e-6,max=3e-4, xra=[-5000,5000],yra=[-5000,5000]
;wset, 13 & err=tvread(/png,/nodialog, filename=sdir+'noise_div_tf_0504')

tv_spt_map, shift(gauss_smooth(ntmp,1000.,5.), 2160, 2160), /norms, reso=5, winnum=14,min=1e-6,max=3e-4, xra=[-5000,5000],yra=[-5000,5000]
;wset, 14 & err=tvread(/png,/nodialog, filename=sdir+'gs1000_noise_div_tf_0504')

tv_spt_map, shift(gauss_smooth(ntmp,300.,5.), 2160, 2160), /norms, reso=5, winnum=15,min=1e-6,max=3e-4, xra=[-5000,5000],yra=[-5000,5000]
;wset, 15 & err=tvread(/png,/nodialog, filename=sdir+'gs300_noise_div_tf_0504')

stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;----------------------
; view this field's weight    
;----------------------
fig_dir = '/home/kstory/public_html/notebook/spt_lps12/'

tv_spt_map,shift(weight_2d,nbig/2.,nbig/2.), /norms, min=0,max=1,reso=5,xr=[-3000,3000],yr=[-3000,3000],title='weight_2d=1 / ( GS{cl_th} + GS{noise/tf} )^2, '+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
;err = tvread(/png, filename=fig_dir+'kweight_'+field_name+'_0504', /nodialog)

; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,xr=[-3000,3000],yr=[-3000,3000],title='un-normalized LOG weight_2d'+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
;err = tvread(/png, filename=fig_dir+'kweight_unnorm_'+field_name+'_0504', /nodialog)

; plot the un-normalized kweight, zoomed in
tv_spt_map,shift(alog(w1),nbig/2.,nbig/2.), /norms,reso=5,title='un-normalized LOG weight_2d'+field_name
; draw cirlces at the edges of the low-ell analysis
tvellipse,650,650,0,0,color=!blue,/data
tvellipse,3000,3000,0,0,color=!blue,/data
;err = tvread(/png, filename=fig_dir+'kweight_unnorm_zoom_'+field_name+'_0504', /nodialog)

; inputs plot
window, 5, xsize=700, ysize=500
vec = indgen(1000)
plot, l[vec,vec], cl_th[vec,vec], /ylog, xtitle='ell', ytitle='Cl [uK-rad]^2', title=field_name, xra=[0,4000]
oplot, l[vec,vec], noise[vec,vec], color=!red
oplot, l[vec,vec], ntmp[vec,vec], color=!orange
legend_str = ['Cl_th', 'noise_psd', 'noise_TF']
legend, legend_str, linestyle=0, colors=[!white, !red, !orange], position=[2500, 0.9e2]
;err = tvread(/png, filename=fig_dir+'kweight_inputs_0504', /nodialog)

stop
END


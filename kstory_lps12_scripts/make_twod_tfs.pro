;;;
; NAME: make_twod_tfs.pro
; PURPOSE:
;   Make two-dimensional Transfer Function for one field
;
; INPUTS:
;   idx,         index in lps12_fieldstruct()
;   stub,        identifier from end2end run, must match what was run in save_mode_coupling_lps12.pro
;
; NOTES:
; 1) expects file at scripts/tf_prep[run].sav
;
; MODIFICATION HISTORY:
;  04/25/2012: (KTS) Created from /home/rkeisler/ps09/create_twodim_tfs.pro
;  05/08/2012: (KTS) Re-written to use tf_prep[run].sav
;  06/05/2012: (KTS) Re-written for run 05, lmax8000 sims
;  06/28/2012: (KTS) Set lmax8000 sims as the default
;;;

PRO make_twod_tfs, idx, run, plotit=plotit, dosave=dosave
compile_opt IDL2, HIDDEN

; retrieve default values for structures going into sav file
print, 'restore, '+'tf_prep'+run+'.sav'
restore,'tf_prep'+run+'.sav'

;------------------------
; Setup
;------------------------
f = lps12_fieldstruct()
info = get_lps12_fieldinfo(0)
nbig = long(info.nbig)
n_sim_coadds = 100

reso = 1.0 ; arcmin per pixel
reso_rad = reso/60.*!dtor

; directories
sim_dir  = '/home/kstory/lps12/sims/'
mask_dir = '/home/kstory/lps12/masks/masks_50mJy/'
tf_dir   = '/home/kstory/lps12/twod_tfs/'

; make some noise
print, 'MAKE_TWOD_TFS: for field ', idx, ', '+f[idx].name

; Theory Dl's
;readcol,'/data/rkeisler/low_ell_sims/input/dl_input_true_20101221_061503.txt',lth,dlth
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2


ex=execute('s=s'+strtrim(idx,2))
field_name = s.field_name
nbig = long(s.nbig)
power = s.power*1d12/s.mask_factor*(nbig*nbig)*(reso_rad^2.)

l = make_fft_grid(reso/60.*!dtor/2./!pi,nbig,nbig,fx=lx,fy=ly)
clth_interp = interpol(cl_uK2,l_vec,l)

; get the beam (get the exact beam used by the sims)
year = f[idx].year
beamfiles = get_lps12_beams(year, l_beam, bl_beam, /sim_lmax8000)
if ((run eq '05')) then beamfiles = get_lps12_beams(year, l_beam, bl_beam, /sim_lmax8000)
if ((run eq '03') or (run eq '04')) then beamfiles = get_lps12_beams(year, l_beam, bl_beam, /sim_lmax4500)
print, 'beam files = ', beamfiles ; DEBUGGING

bl = interpol(bl_beam,l_beam,l)
wh_high = where(l gt max(l_beam), n_high)
; replace the part of the beam that is undefined with the beam value
; at the highest ell that is defined.
last_bl_beam = (bl_beam[n_elements(bl_beam)-1])[0]
if n_high gt 0 then bl[wh_high]=last_bl_beam
    
;------------------------
; get the effective change in the TF from the mode-coupling
; correction.
; mode_couplings, l_mode_coupling is in tf_prep.sav
;   RK: "in hindsight this is probably all kinds of wrong and should
;   be ignored.  just ignore mode-coupling when calculating the 2d TF."
;------------------------
; if 0 then begin
;     mode_coupling = mode_couplings[idx,*]
;     mode_coupling[where(l_mode_coupling gt 5e3)]=1.
;     corr_mc = interpol(mode_coupling,l_mode_coupling,l)
;     wh_high = where(l gt max(l_mode_coupling), n_high)
;     if n_high gt 0 then corr_mc[wh_high] = $
;       mode_coupling[n_elements(mode_coupling)-1]
;     tf = power/clth_interp/bl/bl*corr_mc
;     tf_w_beam = power/clth_interp*corr_mc
; endif

;------------------------
; Calculate the TF
;------------------------
tf = power/clth_interp/bl/bl
tf_w_beam = power/clth_interp

if keyword_set(plotit) then tv_spt_map,shift(tf,nbig/2.,nbig/2.),min=0,max=1,xra=[-4000,4000],yra=[-4000,4000],reso=(l[1]-l[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title=s.field_name

README = 'All TFs are in POWER units, as opposed to TEMPERATURE units.'
if keyword_set(dosave) then begin
    savename = tf_dir+'tf_'+field_name+'.sav'
    print, 'Saving twod_TF: ', savename
    save, $
      l,lx,ly,tf,tf_w_beam,README,beamfiles, cl_uK2, clth_interp, $
      filename=savename
endif

END



